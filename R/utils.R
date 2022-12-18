#' Prepare SCP python virtual environment
#'
#' @param python_version The version of python to install. Default is \code{3.8}
#' @param conda_binary
#' @param miniconda_repo
#' @param force
#'
#' @export
PrepareEnv <- function(conda_binary = NULL, python_version = "3.8",
                       miniconda_repo = "https://repo.anaconda.com/miniconda",
                       force = FALSE, ...) {
  options(reticulate.conda_binary = conda_binary)
  if (python_version < numeric_version("3.7.0") || python_version >= numeric_version("3.10.0")) {
    stop("SCP currently only support python version 3.7-3.9!")
  }

  conda_exist <- reticulate:::conda_installed()
  if (conda_exist) {
    env_exist <- file.exists(paste0(reticulate:::conda_info()$conda_prefix, "/envs/SCP"))
  } else {
    env_exist <- FALSE
  }
  if (isTRUE(force) && isTRUE(env_exist)) {
    env_exist <- FALSE
  }

  if (isFALSE(env_exist)) {
    if (isFALSE(conda_exist)) {
      options(timeout = 360)
      version <- "3"
      info <- as.list(Sys.info())
      if (info$sysname == "Darwin" && info$machine == "arm64") {
        base <- "https://github.com/conda-forge/miniforge/releases/latest/download"
        name <- "Miniforge3-MacOSX-arm64.sh"
        return(file.path(base, name))
      }
      base <- miniconda_repo
      info <- as.list(Sys.info())
      arch <- reticulate:::miniconda_installer_arch(info)
      version <- as.character(version)
      name <- if (reticulate:::is_windows()) {
        sprintf(
          "Miniconda%s-latest-Windows-%s.exe", version,
          arch
        )
      } else if (reticulate:::is_osx()) {
        sprintf(
          "Miniconda%s-latest-MacOSX-%s.sh", version,
          arch
        )
      } else if (reticulate:::is_linux()) {
        sprintf("Miniconda%s-latest-Linux-%s.sh", version, arch)
      } else {
        stopf("unsupported platform %s", shQuote(Sys.info()[["sysname"]]))
      }
      url <- file.path(base, name)
      options(reticulate.miniconda.url = url)
      reticulate::install_miniconda(path = reticulate::miniconda_path(), force = TRUE, update = FALSE)
    }
    python_path <- reticulate::conda_create(envname = "SCP", python_version = python_version)
  }

  packages <- c(
    "numpy==1.21.6", "numba==0.55.2", "python-igraph==0.10.2",
    "scipy", "pandas", "matplotlib", "versioned-hdf5", "leidenalg", "scanpy", "scvelo", "palantir"
  )
  check_Python(packages = packages, envname = "SCP", pip = TRUE, ...)

  python_path <- reticulate::conda_python("SCP")
  reticulate::use_python(python_path, required = TRUE)

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

#' Check and install R packages
#'
#' @param pkgs Package to be installed. Package source can be CRAN, Bioconductor or Github, e.g. scmap, davidsjoberg/ggsankey.
#' @param pkg_names The name of the package that corresponds to the \code{pkgs} parameter, used to check if the package is already installed.
#' By default, the package name is extracted according to the \code{pkgs} parameter.
#' @param install_methods Functions for installing R packages. The default is to try to install packages from CRAN, Bioconductor and Github.
#' @param lib Character vector giving the library directories where to install the packages.
#' @param force Whether to force package installation. Default is FALSE.
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
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager", lib = lib)
            }
            eval(str2lang(paste0(install_methods[i], "('", pkg, "', lib='", lib, "', update = FALSE, upgrade = 'never', ask = FALSE, force = TRUE)")))
          } else if (grepl("devtools", install_methods[i])) {
            if (!requireNamespace("devtools", quietly = TRUE)) {
              install.packages("devtools", lib = lib)
            }
            if (!requireNamespace("withr", quietly = TRUE)) {
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

#' Check if the python package exists in the environment
#'
#' @param packages
#' @param envname The name of, or path to, a Python virtual environment.
#'
#' @export
exist_Python_pkgs <- function(packages, envname = "SCP") {
  conda_exist <- reticulate:::conda_installed()
  if (conda_exist) {
    env_exist <- file.exists(paste0(reticulate:::conda_info()$conda_prefix, "/envs/", envname))
  } else {
    env_exist <- FALSE
  }

  if (isFALSE(env_exist)) {
    stop("Can not find the conda environment: ", envname)
  }
  all_installed <- reticulate:::conda_list_packages(envname, no_pip = FALSE)
  packages_installed <- NULL
  for (pkg in packages) {
    pkg_info <- strsplit(pkg, split = "==")[[1]]
    pkg_name <- pkg_info[1]
    pkg_version <- pkg_info[2]
    if (pkg_name %in% all_installed$package) {
      if (!is.na(pkg_version)) {
        packages_installed[pkg] <- all_installed$version[all_installed$package == pkg_name] == pkg_version
      } else {
        packages_installed[pkg] <- TRUE
      }
    } else {
      packages_installed[pkg] <- FALSE
    }
  }
  return(packages_installed)
}

#' Check and install python packages
#'
#' @param envname The name of, or path to, a Python virtual environment.
#' @param force Whether to force package installation. Default is FALSE.
#' @param packages
#' @param pip
#' @param ...
#'
#' @importFrom rlang %||%
#' @export
check_Python <- function(packages, envname = "SCP", force = FALSE, pip = TRUE, ...) {
  conda_exist <- reticulate:::conda_installed()
  if (conda_exist) {
    env_exist <- file.exists(paste0(reticulate:::conda_info()$conda_prefix, "/envs/", envname))
  } else {
    env_exist <- FALSE
  }
  if (isFALSE(env_exist)) {
    warning(envname, " python environment do not exist. Create it with the PrepareEnv function...", immediate. = TRUE)
    PrepareEnv()
  }
  if (isTRUE(force)) {
    exist <- setNames(rep(FALSE, length(packages)), packages)
  } else {
    exist <- exist_Python_pkgs(packages = packages, envname = envname)
  }
  if (sum(!exist) > 0) {
    pkgs_to_install <- names(exist)[!exist]
    message("Try to install ", paste0(pkgs_to_install, collapse = ","), " ...")
    if (isTRUE(pip)) {
      pkgs_to_install <- c("pip", pkgs_to_install)
    }
    tryCatch(expr = {
      reticulate::conda_install(pkgs_to_install, envname = envname, pip = pip, ...)
    }, error = identity)
  }

  exist <- exist_Python_pkgs(packages = packages, envname = envname)
  if (sum(!exist) > 0) {
    stop("Failed to install the package(s): ", paste0(names(exist)[!exist], collapse = ","), " into the environment '", envname, "'. Please install manually.")
  } else {
    return(invisible(NULL))
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
  }, error = function(error) {
    message(error)
    stop("Failed to run '", command, "'. Please check manually.")
  })
}

#' @param url
#'
#' @param destfile
#' @param methods
#' @param quiet
#' @param attempts
#' @param return_status
#' @param mode
#'
#' @importFrom utils download.file
#' @export
download <- function(url, destfile, methods = c("auto", "wget", "libcurl", "curl", "wininet", "internal"), quiet = FALSE, attempts = 2, return_status = FALSE, mode = "w") {
  if (missing(url) || missing(destfile)) {
    stop("'url' and 'destfile' must be both provided.")
  }
  attempts <- 0
  status <- NULL
  while (is.null(status)) {
    for (method in methods) {
      status <- tryCatch(expr = {
        suppressWarnings(download.file(url, destfile = destfile, method = method, quiet = quiet, mode = mode))
        status <- 1
      }, error = function(error) {
        message(error)
        message("Cannot download from the url: ", url)
        message("Failed to download using '", method, "'. Retry...")
        Sys.sleep(1)
        return(NULL)
      })
      if (!is.null(status)) {
        break
      }
    }
    attempts <- attempts + 1
    if (is.null(status) && attempts >= attempts) {
      stop("Download failed.")
    }
  }
  if (isTRUE(return_status)) {
    return(status)
  } else {
    return(invisible(NULL))
  }
}

kegg_get <- function(url) {
  temp <- tempfile()
  download(url = url, destfile = temp)
  content <- readLines(temp) %>%
    strsplit(., "\t") %>%
    do.call("rbind.data.frame", .)
  unlink(temp)
  return(content)
}

rescale <- function(x, from = range(x, na.rm = TRUE, finite = TRUE), to = c(0, 1)) {
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  } else {
    return((x - from[1]) / diff(from) * diff(to) + to[1])
  }
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
