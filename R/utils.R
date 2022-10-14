#' Prepare SCP python virtual environment
#'
#' @param python Python path which is used to create the virtual environment.
#' @param pipy_mirror pipy mirrors. Default is "https://pypi.org/simple/". Options can be "https://pypi.tuna.tsinghua.edu.cn/simple", "http://mirrors.aliyun.com/pypi/simple/", "https://pypi.mirrors.ustc.edu.cn/simple/", etc.
#' @param remove_old Whether to remove the old SCP virtual environment.
#' @export
PrepareVirtualEnv <- function(python = NULL, pipy_mirror = "https://pypi.org/simple/", remove_old = FALSE) {
  if (isTRUE(remove_old) || reticulate::virtualenv_exists("SCP")) {
    if (isTRUE(reticulate:::is_python_initialized())) {
      warning("Package installation may fail because python is initialized in the current R session. You may run SCP::PrepareVirtualEnv() in a new R seesion to create SCP virtual environmenet. If you want to disable virtual environment initialization, you can set options(SCP_virtualenv_init = FALSE) before library(SCP).", immediate. = TRUE)
    }
  }
  if (isTRUE(remove_old)) {
    unlink(reticulate:::virtualenv_path("SCP"), recursive = TRUE)
  }
  if (!is.null(python) && file.exists(python)) {
    pythons <- python
  } else {
    pythons <- unique(c(reticulate::virtualenv_python("SCP"), Sys.getenv("RETICULATE_PYTHON"), Sys.which("python3"), Sys.which("python")))
  }
  python_path <- NULL
  sys_bit <- ifelse(grepl("64", Sys.info()["machine"]), "64bit", "32bit")
  for (py in pythons) {
    if (py != "" && !file.exists(py)) {
      # packageStartupMessage(py, " is not a Python executable file.")
      next
    }
    # py <- file.path(normalizePath(dirname(py)), basename(py))
    py_version <- tryCatch(suppressWarnings(reticulate:::python_version(py)),
      error = identity
    )
    if (inherits(py_version, "error") || length(py_version) == 0) {
      next
    }
    if (py_version < numeric_version("3.7.0") || py_version >= numeric_version("3.10.0")) {
      next
    }
    py_bit <- tryCatch(suppressWarnings(system2(command = py, args = " -c \"import platform; print(platform.architecture()[0])\"", stdout = TRUE)),
      error = identity
    )
    if (inherits(py_bit, "error") || length(py_bit) == 0) {
      next
    }
    if (identical(py_bit, sys_bit)) {
      python_path <- py
      # packageStartupMessage("python path: ", python_path)
      break
    } else {
      # packageStartupMessage("System architecture is ", sys_bit, " but ", py, " is ", py_bit)
      next
    }
  }

  if (!reticulate::virtualenv_exists("SCP")) {
    if (is.null(python_path)) {
      packageStartupMessage("Python(3.7-3.9) is unavailable. Install python(3.8.8) automatically ...")
      git_exist <- suppressWarnings(system("git", ignore.stdout = TRUE, ignore.stderr = TRUE))
      if (git_exist == 127) {
        stop("You need to install git first! (http://git-scm.com/download/) or install python manually.")
      }
      python_path <- reticulate::install_python(version = ifelse(sys_bit == "64bit", "3.8.8", "3.8.8-win32"))
    }
    packageStartupMessage("Create SCP virtual environment. The path is: ", reticulate:::virtualenv_path("SCP"))
    reticulate::virtualenv_create(envname = "SCP", python = python_path, packages = FALSE)
    check_Python(
      pkgs = c("pip", "setuptools", "wheel"),
      envname = "SCP",
      pipy_mirror = pipy_mirror,
      force = TRUE
    )
  }
  python_path <- reticulate::virtualenv_python("SCP")
  Sys.setenv(RETICULATE_PYTHON = python_path)

  version <- tryCatch(suppressWarnings(reticulate:::python_version(python_path)), error = identity)
  if (inherits(version, "error")) {
    stop("SCP need python 3.7-3.9! Please install python and reload the SCP!")
  } else {
    if (version < numeric_version("3.7.0") || version >= numeric_version("3.10.0")) {
      stop(
        "SCP currently only support python version 3.7-3.9! The version of Python currently is ", version,
        "!\nPython related functions may not work. Please install the right version of python and reload the SCP!"
      )
    } else {
      reticulate::use_virtualenv("SCP", required = TRUE)

      if (!exist_pkg("pip")) {
        temp <- tempfile()
        download.file("https://bootstrap.pypa.io/get-pip.py", temp)
        suppressWarnings(system2(command = reticulate::virtualenv_python("SCP"), args = temp, stdout = TRUE))
        unlink(temp)
      }
      # reticulate::py_install("numpy==1.23.2", pip_options = "--no-binary='numpy'", ignore_installed = TRUE, envname = "SCP")
      # check_Python(pkgs = c("numba==0.53.1", "python-igraph==0.9.9", "pandas", "matplotlib", "versioned-hdf5", "scanpy", "scvelo", "palantir"), envname = "SCP")
      check_Python(
        pkgs = c(
          "numpy==1.21.6", "numba==0.55.2", "python-igraph==0.9.11",
          "pandas", "matplotlib", "versioned-hdf5", "scanpy", "scvelo", "palantir"
        ),
        envname = "SCP",
        pipy_mirror = pipy_mirror
      )

      pyinfo <- utils::capture.output(reticulate::py_config())
      pyinfo_mesg <- c(
        "======================== SCP python config ========================",
        pyinfo,
        "==================================================================="
      )
      invisible(lapply(pyinfo_mesg, packageStartupMessage))
      invisible(run_Python(command = "import matplotlib", envir = .GlobalEnv))
      if (!interactive()) {
        invisible(run_Python(command = "matplotlib.use('pdf')", envir = .GlobalEnv))
      }
      invisible(run_Python(command = "import matplotlib.pyplot as plt", envir = .GlobalEnv))
      invisible(run_Python(command = "import scanpy", envir = .GlobalEnv))
    }
  }
}

#' Check and install R packages
#'
#' @param pkgs
#'
#' @param pkg_names
#' @param install_methods
#' @param lib
#' @param force
#'
#' @importFrom rlang %||%
#' @importFrom utils packageVersion
#' @export
check_R <- function(pkgs, pkg_names = NULL, install_methods = c("BiocManager::install", "install.packages", "devtools::install_github"), lib = .libPaths()[1], force = FALSE) {
  if (length(pkg_names) != 0 && length(pkg_names) != length(pkgs)) {
    stop("pkg_names must be NULL or a vector of the same length with pkgs")
  }
  status_list <- list()
  for (n in seq_along(pkgs)) {
    pkg <- pkgs[n]
    pkg_info <- pkg
    if (!grepl("/", pkg_info)) {
      pkg_info <- paste0("/", pkg_info)
    }
    if (!grepl("@", pkg_info)) {
      pkg_info <- paste0(pkg_info, "@")
    }
    git <- grep("/", sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\1", x = pkg_info), value = TRUE)
    git <- gsub("/", "", git)
    pkg_name <- pkg_names[n] %||% sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\2", x = pkg_info)
    version <- grep("@", sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\3", x = pkg_info), value = TRUE)
    version <- gsub("@", "", version)
    if (version != "") {
      force <- isTRUE(packageVersion(pkg_name) < package_version(version))
    }
    if (!suppressPackageStartupMessages(requireNamespace(pkg_name, quietly = TRUE)) || isTRUE(force)) {
      message("Install package: '", pkg_name, "' ...")
      status_list[[pkg]] <- FALSE
      i <- 1
      while (!isTRUE(status_list[[pkg]])) {
        tryCatch(expr = {
          if (grepl("BiocManager", install_methods[i])) {
            if (!require("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager", lib = lib)
            }
            eval(str2lang(paste0(install_methods[i], "('", pkg, "', lib='", lib, "', update = FALSE, upgrade = 'never', ask = FALSE, force = TRUE)")))
          } else if (grepl("devtools", install_methods[i])) {
            if (!require("devtools", quietly = TRUE)) {
              install.packages("devtools", lib = lib)
            }
            if (!require("withr", quietly = TRUE)) {
              install.packages("withr", lib = lib)
            }
            eval(str2lang(paste0("withr::with_libpaths(new = '", lib, "', ", install_methods[i], "('", pkg, "', upgrade = 'never', force = TRUE))")))
          } else {
            eval(str2lang(paste0(install_methods[i], "('", pkg, "', lib='", lib, "', force = TRUE)")))
          }
        }, error = function(e) {
          status_list[[pkg]] <- FALSE
        })
        if (requireNamespace(pkg_name, quietly = TRUE)) {
          status_list[[pkg]] <- TRUE
        }
        i <- i + 1
        if (i > length(install_methods)) {
          break
        }
      }
    } else {
      status_list[[pkg]] <- TRUE
    }
  }
  out <- sapply(status_list, isTRUE)
  out <- out[!out]
  if (length(out) > 0) {
    stop("Failed to install the package(s): ", paste0(names(out), collapse = ","), ". Please install manually.")
  }
}

#' @export
exist_pkg <- function(pkg, envname = "SCP") {
  pkg_info <- strsplit(pkg, split = "==")[[1]]
  pkg_name <- pkg_info[1]
  pkg_version <- pkg_info[2]
  pkg_exist <- suppressWarnings(system2(command = reticulate::virtualenv_python(envname), args = paste0("-m pip show ", pkg_name), stdout = TRUE))
  if (length(pkg_exist) == 0 || isTRUE(attr(pkg_exist, "status") == 1)) {
    return(FALSE)
  } else {
    if (!is.na(pkg_version)) {
      pkg_installed_version <- strsplit(pkg_exist, split = "Version: ")[[2]][2]
      if (pkg_installed_version == pkg_version) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(TRUE)
    }
  }
}

#' Check and install python packages
#'
#' @param pkgs
#' @param pkg_names
#' @param envname
#' @param force
#'
#' @importFrom rlang %||%
#' @export
check_Python <- function(pkgs, pkg_names = NULL, envname = "SCP", pipy_mirror = "https://pypi.org/simple/", force = FALSE) {
  if (!reticulate::virtualenv_exists("SCP")) {
    warning("SCP python virtual environment do not exist. Create it with the PrepareVirtualEnv function...", immediate. = TRUE)
    PrepareVirtualEnv()
  }
  if (length(pkg_names) != 0 && length(pkg_names) != length(pkgs)) {
    stop("pkg_names must be NULL or a vector of the same length with pkgs")
  }
  status_list <- list()
  for (n in seq_along(pkgs)) {
    pkg <- pkgs[n]
    pkg_name <- pkg_names[n] %||% gsub("(.*)(\\[.*\\])", "\\1", pkg)
    exist <- exist_pkg(pkg_name, envname)
    if (!isTRUE(exist) || isTRUE(force)) {
      message("Try to install '", pkg, "' ...")
      tryCatch(expr = {
        reticulate::py_install(pkg, envname = envname, pip_options = paste("-i", pipy_mirror))
      }, error = function(e) {
        warning("Something went wrong when installing the package ", pkg)
      })
    }
  }
  for (n in seq_along(pkgs)) {
    pkg <- pkgs[n]
    pkg_name <- pkg_names[n] %||% gsub("(.*)(\\[.*\\])", "\\1", pkg)
    exist <- exist_pkg(pkg_name, envname)
    if (isTRUE(exist)) {
      status_list[[pkg_name]] <- TRUE
    } else {
      status_list[[pkg_name]] <- FALSE
    }
  }
  out <- sapply(status_list, isTRUE)
  out <- out[!out]
  if (length(out) > 0) {
    stop("Failed to install the module(s): ", paste0(names(out), collapse = ","), " into the environment '", envname, "'. Please install manually.")
  }
}

#' @importFrom rlang %||%
#' @export
list_palette <- function(palette_list, names = NULL) {
  if (is.null(names)) {
    names <- names %||% names(palette_list) %||% seq_along(palette_list)
  }
  par(mar = c(0, 0, 0, 0) + 0.1)

  plot(0, 0,
    type = "n", axes = FALSE, bty = "n", xlab = "", ylab = "",
    xlim = c(0, 1), ylim = c(-length(palette_list) - 1, -1)
  )

  for (i in seq_len(length(palette_list))) {
    colors_len <- length(palette_list[[i]])
    breaks <- seq(from = 0, to = 1, length = colors_len + 1)


    text(0, -i, names[i], pos = 4)
    rect(
      xleft = breaks[1:colors_len], xright = breaks[1:colors_len + 1],
      ytop = -0.15 - i, ybottom = -0.8 - i,
      col = palette_list[[i]], border = NA
    )
  }
}

#' @export
run_Python <- function(command, envir = .GlobalEnv) {
  tryCatch(expr = {
    eval(
      {
        reticulate::py_run_string(command)
      },
      envir = envir
    )
  }, error = function(e) {
    message(e)
    stop("Failed to run '", command, "'. Please check manually.")
  })
}

#' @importFrom utils download.file
kegg_get <- function(url) {
  temp <- tempfile()
  ntry <- 0
  status <- NULL
  while (is.null(status)) {
    status <- tryCatch(expr = {
      download.file(url, destfile = temp, method = "wget", quiet = TRUE)
    }, error = function(e) {
      message("Get errors when connecting with KEGG...\nRetrying...")
      Sys.sleep(1)
      return(NULL)
    })
    ntry <- ntry + 1
    if (is.null(status) && ntry >= 5) {
      stop("Stop connecting...")
    }
  }
  content <- readLines(temp) %>%
    strsplit(., "\t") %>%
    do.call("rbind", .)
  res <- data.frame(from = content[, 1], to = content[, 2])
  unlink(temp)
  return(res)
}

rescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  }
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
  if (length(x) == 1) {
    return(TRUE)
  }
  if (length(x) != 2) {
    stop("x must be length 1 or 2")
  }
  if (any(is.na(x))) {
    return(NA)
  }
  if (x[1] == x[2]) {
    return(TRUE)
  }
  if (all(is.infinite(x))) {
    return(FALSE)
  }
  m <- min(abs(x))
  if (m == 0) {
    return(FALSE)
  }
  abs((x[1] - x[2]) / m) < tol
}

#' @importFrom grDevices col2rgb rgb
col2hex <- function(cname) {
  colMat <- col2rgb(cname)
  rgb(red = colMat[1, ] / 255, green = colMat[2, ] / 255, blue = colMat[3, ] / 255)
}

#' @export
#' @importFrom rlang caller_env is_null is_scalar_character is_character is_function set_names env env_get env_bind syms call2
invoke <- function(.fn, .args = list(), ..., .env = caller_env(), .bury = c(".fn", "")) {
  args <- c(.args, list(...))
  if (is_null(.bury) || !length(args)) {
    if (is_scalar_character(.fn)) {
      .fn <- env_get(.env, .fn, inherit = TRUE)
    }
    call <- call2(.fn, !!!args)
    return(.External2(rlang:::ffi_eval, call, .env))
  }
  if (!is_character(.bury, 2L)) {
    abort("`.bury` must be a character vector of length 2")
  }
  arg_prefix <- .bury[[2]]
  fn_nm <- .bury[[1]]
  buried_nms <- paste0(arg_prefix, seq_along(args))
  buried_args <- set_names(args, buried_nms)
  .env <- env(.env, !!!buried_args)
  args <- set_names(buried_nms, names(args))
  args <- syms(args)
  if (is_function(.fn)) {
    env_bind(.env, `:=`(!!fn_nm, .fn))
    .fn <- fn_nm
  }
  call <- call2(.fn, !!!args)
  .External2(rlang:::ffi_eval, call, .env)
}

#' @export
unnest <- function(data, cols, keep_empty = FALSE) {
  if (nrow(data) == 0 || length(cols) == 0) {
    return(data)
  }
  for (col in cols) {
    col_expand <- unlist(data[[col]])
    expand_times <- sapply(data[[col]], length)
    if (isTRUE(keep_empty)) {
      data[[col]][expand_times == 0] <- NA
      col_expand <- unlist(data[[col]])
      expand_times[expand_times == 0] <- 1
    }
    data <- data[rep(seq_len(nrow(data)), times = expand_times), ]
    data[, col] <- col_expand
  }
  rownames(data) <- NULL
  return(data)
}

#' @useDynLib SCP
#' @export
as_matrix <- function(mat) {
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1])

  tmp <- asMatrix(
    rp = row_pos, cp = col_pos, z = mat@x,
    nrows = mat@Dim[1], ncols = mat@Dim[2]
  )

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
