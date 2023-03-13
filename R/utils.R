#' Prepare SCP python environment
#'
#' @param python_version The version of python to install. Default is \code{3.8}
#' @param miniconda_repo  Repositories for miniconda. Default is \code{https://repo.anaconda.com/miniconda}
#' @param force Whether to force a new environment to be created. If \code{TRUE}, the existing environment will be recreated. Default is \code{FALSE}
#' @inheritParams check_Python
#'
#' @export
PrepareEnv <- function(conda = "auto", miniconda_repo = "https://repo.anaconda.com/miniconda",
                       python_version = "3.8", envname = NULL, force = FALSE, ...) {
  envname <- get_envname(envname)

  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }

  if (is.null(conda)) {
    env <- FALSE
  } else {
    envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    env <- env_exist(conda = conda, envname = envname, envs_dir = envs_dir)
    if (isTRUE(force) && isTRUE(env)) {
      unlink(paste0(envs_dir, "/", envname), recursive = TRUE)
      env <- FALSE
    }
  }

  if (isTRUE(env)) {
    python_path <- conda_python(conda = conda, envname = envname)
    installed_python_version <- reticulate:::python_version(python_path)
    if (installed_python_version < numeric_version("3.7.0") || installed_python_version >= numeric_version("3.10.0")) {
      stop("The python version in the installed SCP environment does not match the requirements. You need to recreate the SCP environment.")
    }
  } else {
    force <- TRUE
    if (is.null(conda)) {
      message("Conda not found. Installing miniconda...")
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
      if (!is.na(Sys.getenv("USER", unset = NA))) {
        miniconda_path <- gsub(pattern = "\\$USER", replacement = Sys.getenv("USER"), reticulate::miniconda_path())
      } else {
        miniconda_path <- reticulate::miniconda_path()
      }
      unlink(miniconda_path, recursive = TRUE)
      reticulate::install_miniconda(path = miniconda_path, force = TRUE, update = FALSE)
      conda <- reticulate:::miniconda_conda(miniconda_path)
      envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    }
    if (python_version < numeric_version("3.7.0") || python_version >= numeric_version("3.10.0")) {
      stop("SCP currently only support python version 3.7-3.9!")
    }
    python_path <- reticulate::conda_create(conda = conda, envname = envname, python_version = python_version)
    env_path <- paste0(envs_dir, "/", envname)
    env <- file.exists(env_path)
    if (isFALSE(env)) {
      print(reticulate:::conda_info(conda = conda))
      print(reticulate::conda_list(conda = conda))
      stop(
        "Unable to find SCP environment under the expected path: ", env_path, "\n",
        "conda: ", conda, "\n",
        "SCP python: ", python_path, "\n"
      )
    }
  }

  packages <- c(
    "numpy==1.21.6", "numba==0.55.2", "scikit-learn==1.1.2", "pandas==1.3.5", "python-igraph==0.10.2", "matplotlib==3.6.3",
    "scipy", "versioned-hdf5", "leidenalg", "scanpy", "scvelo", "palantir"
  )
  check_Python(packages = packages, envname = envname, conda = conda, force = force, ...)

  Sys.unsetenv("RETICULATE_PYTHON")
  python_path <- conda_python(conda = conda, envname = envname)
  reticulate::use_python(python_path, required = TRUE)

  pyinfo <- utils::capture.output(reticulate::py_config())
  pyinfo_mesg <- c(
    "====================== SCP conda environment ======================",
    paste0("conda: ", conda),
    paste0("environment: ", paste0(envs_dir, "/", get_envname())),
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
#' @param packages Package to be installed. Package source can be CRAN, Bioconductor or Github, e.g. scmap, davidsjoberg/ggsankey.
#' @param package_names The name of the package that corresponds to the \code{packages} parameter, used to check if the package is already installed.
#' By default, the package name is extracted according to the \code{packages} parameter.
#' @param install_methods Functions used to install R packages.
#' @param lib  The location of the library directories where to install the packages.
#' @param force Whether to force the installation of packages. Default is \code{FALSE}.
#'
#' @importFrom rlang %||%
#' @importFrom utils packageVersion
#' @export
check_R <- function(packages, package_names = NULL, install_methods = c("BiocManager::install", "install.packages", "devtools::install_github"), lib = .libPaths()[1], force = FALSE) {
  if (length(package_names) != 0 && length(package_names) != length(packages)) {
    stop("package_names must be NULL or a vector of the same length with packages")
  }
  status_list <- list()
  for (n in seq_along(packages)) {
    pkg <- packages[n]
    pkg_info <- pkg
    if (!grepl("/", pkg_info)) {
      pkg_info <- paste0("/", pkg_info)
    }
    if (!grepl("@", pkg_info)) {
      pkg_info <- paste0(pkg_info, "@")
    }
    git <- grep("/", sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\1", x = pkg_info), value = TRUE)
    git <- gsub("/", "", git)
    pkg_name <- package_names[n] %||% sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\2", x = pkg_info)
    version <- grep("@", sub(pattern = "(.*/)(.*)(@.*)", replacement = "\\3", x = pkg_info), value = TRUE)
    version <- gsub("@", "", version)
    if (version != "") {
      force_update <- isTRUE(packageVersion(pkg_name) < package_version(version)) || isTRUE(force)
    } else {
      force_update <- isTRUE(force)
    }
    if (!suppressPackageStartupMessages(requireNamespace(pkg_name, quietly = TRUE)) || isTRUE(force_update)) {
      message("Install package: \"", pkg_name, "\" ...")
      status_list[[pkg]] <- FALSE
      i <- 1
      while (isFALSE(status_list[[pkg]])) {
        tryCatch(expr = {
          if (grepl("BiocManager", install_methods[i])) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager", lib = lib)
            }
            eval(str2lang(paste0(install_methods[i], "(\"", pkg, "\", lib=\"", lib, "\", update = FALSE, upgrade = \"never\", ask = FALSE, force = TRUE)")))
          } else if (grepl("devtools", install_methods[i])) {
            if (!requireNamespace("devtools", quietly = TRUE)) {
              install.packages("devtools", lib = lib)
            }
            if (!requireNamespace("withr", quietly = TRUE)) {
              install.packages("withr", lib = lib)
            }
            eval(str2lang(paste0("withr::with_libpaths(new = \"", lib, "\", ", install_methods[i], "(\"", pkg, "\", upgrade = \"never\", force = TRUE))")))
          } else {
            eval(str2lang(paste0(install_methods[i], "(\"", pkg, "\", lib=\"", lib, "\", force = TRUE)")))
          }
        }, error = function(e) {
          status_list[[pkg]] <- FALSE
        })
        if (version == "") {
          status_list[[pkg]] <- requireNamespace(pkg_name, quietly = TRUE)
        } else {
          if (requireNamespace(pkg_name, quietly = TRUE)) {
            status_list[[pkg]] <- packageVersion(pkg_name) >= package_version(version)
          } else {
            status_list[[pkg]] <- FALSE
          }
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
#' @inheritParams check_Python
#' @export
exist_Python_pkgs <- function(packages, envname = NULL, conda = "auto") {
  envname <- get_envname(envname)
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    stop("Can not find the conda environment: ", envname)
  }
  all_installed <- reticulate:::conda_list_packages(conda = conda, envname = envname, no_pip = FALSE)
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
#' @param packages A character vector, indicating package names which should be installed or removed. Use \code{⁠<package>==<version>}⁠ to request the installation of a specific version of a package.
#' @param envname The name of a conda environment.
#' @param conda The path to a conda executable. Use \code{"auto"} to allow SCP to automatically find an appropriate conda binary.
#' @param force Whether to force package installation. Default is \code{FALSE}.
#' @param pip Whether to use pip for package installation. By default, packages are installed from the active conda channels.
#' @param pip_options An optional character vector of additional command line arguments to be passed to \code{pip}. Only relevant when \code{pip = TRUE}.
#' @param ... Other arguments passed to \code{\link[reticulate]{conda_install}}
#'
#' @examples
#' check_Python(packages = c("bbknn", "scanorama"))
#' \dontrun{
#' check_Python(packages = "scvi-tools==0.20.0", envname = "SCP_env", pip_options = "-i https://pypi.tuna.tsinghua.edu.cn/simple")
#' }
#' @export
check_Python <- function(packages, envname = NULL, conda = "auto", force = FALSE, pip = TRUE, pip_options = character(), ...) {
  envname <- get_envname(envname)
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    warning(envname, " python environment does not exist. Create it with the PrepareEnv function...", immediate. = TRUE)
    PrepareEnv()
  }
  if (isTRUE(force)) {
    pkg_installed <- setNames(rep(FALSE, length(packages)), packages)
    pip_options <- c(pip_options, "--force-reinstall")
  } else {
    pkg_installed <- exist_Python_pkgs(packages = packages, envname = envname, conda = conda)
  }
  if (sum(!pkg_installed) > 0) {
    pkgs_to_install <- names(pkg_installed)[!pkg_installed]
    message("Try to install ", paste0(pkgs_to_install, collapse = ","), " ...")
    if (isTRUE(pip)) {
      pkgs_to_install <- c("pip", pkgs_to_install)
    }
    tryCatch(expr = {
      conda_install(conda = conda, packages = pkgs_to_install, envname = envname, pip = pip, pip_options = pip_options, ...)
    }, error = identity)
  }

  pkg_installed <- exist_Python_pkgs(packages = packages, envname = envname, conda = conda)
  if (sum(!pkg_installed) > 0) {
    stop("Failed to install the package(s): ", paste0(names(pkg_installed)[!pkg_installed], collapse = ","), " into the environment \"", envname, "\". Please install manually.")
  } else {
    return(invisible(NULL))
  }
}

#' Check if a conda environment exists
#'
#' @param envs_dir Directories in which conda environments are located.
#' @inheritParams check_Python
env_exist <- function(conda = "auto", envname = NULL, envs_dir = NULL) {
  envname <- get_envname(envname)
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  if (!is.null(conda)) {
    if (is.null(envs_dir)) {
      envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    }
    exist <- file.exists(paste0(envs_dir, "/", envname))
  } else {
    exist <- FALSE
  }
  return(exist)
}

get_envname <- function(envname = NULL) {
  if (is.character(envname)) {
    envname <- envname
  } else {
    envname <- getOption("SCP_env_name", default = "SCP_env")
  }
  return(envname)
}

#' Find an appropriate conda binary
#'
#' @export
find_conda <- function() {
  conda <- tryCatch(reticulate::conda_binary(conda = "auto"), error = identity)
  conda_exist <- !inherits(conda, "error")
  if (isFALSE(conda_exist)) {
    if (!is.na(Sys.getenv("USER", unset = NA))) {
      miniconda_path <- gsub(pattern = "\\$USER", replacement = Sys.getenv("USER"), reticulate::miniconda_path())
    } else {
      miniconda_path <- reticulate::miniconda_path()
    }
    conda_exist <- reticulate:::miniconda_exists(miniconda_path) && reticulate:::miniconda_test(miniconda_path)
    if (isTRUE(conda_exist)) {
      conda <- reticulate:::miniconda_conda(miniconda_path)
    } else {
      conda <- NULL
    }
  }
  return(conda)
}

#' Installs a list of packages into a specified conda environment
#'
#' @inheritParams reticulate::conda_install
#' @importFrom rlang %||%
conda_install <- function(envname = NULL, packages, forge = TRUE, channel = character(),
                          pip = FALSE, pip_options = character(), pip_ignore_installed = FALSE,
                          conda = "auto", python_version = NULL, ...) {
  envname <- get_envname(envname)
  reticulate:::check_forbidden_install("Python packages")
  if (missing(packages)) {
    if (!is.null(envname)) {
      fmt <- paste("argument \"packages\" is missing, with no default",
        "- did you mean 'conda_install(<envname>, %1$s)'?",
        "- use 'py_install(%1$s)' to install into the active Python environment",
        sep = "\n"
      )
      stopf(fmt, deparse1(substitute(envname)), call. = FALSE)
    } else {
      packages
    }
  }
  conda <- reticulate::conda_binary(conda)
  envname <- reticulate:::condaenv_resolve(envname)
  python_package <- if (is.null(python_version)) {
    NULL
  } else if (grepl("[><=]", python_version)) {
    paste0("python", python_version)
  } else {
    sprintf("python=%s", python_version)
  }
  python <- tryCatch(conda_python(envname = envname, conda = conda), error = identity)
  if (inherits(python, "error") || !file.exists(python)) {
    reticulate::conda_create(envname = envname, packages = python_package %||% "python", forge = forge, channel = channel, conda = conda)
    python <- conda_python(envname = envname, conda = conda)
  }
  if (!is.null(python_version)) {
    args <- reticulate:::conda_args("install", envname, python_package)
    status <- reticulate:::system2t(conda, shQuote(args))
    if (status != 0L) {
      fmt <- "installation of '%s' into environment '%s' failed [error code %i]"
      msg <- sprintf(fmt, python_package, envname, status)
      stop(msg, call. = FALSE)
    }
  }
  if (pip) {
    # target_dir <- system2(command = python, args = c("-c \"import site; print(site.getsitepackages()[0])\""), stdout = TRUE)
    # pip_options <- c(pip_options, paste("--target", target_dir))
    result <- reticulate:::pip_install(
      python = python, packages = packages,
      pip_options = pip_options, ignore_installed = pip_ignore_installed,
      conda = conda, envname = envname
    )
    return(result)
  }
  args <- reticulate:::conda_args("install", envname)
  channels <- if (length(channel)) {
    channel
  } else if (forge) {
    "conda-forge"
  }
  for (ch in channels) args <- c(args, "-c", ch)
  args <- c(args, python_package, packages)
  result <- reticulate:::system2t(conda, shQuote(args))
  if (result != 0L) {
    fmt <- "one or more Python packages failed to install [error code %i]"
    stopf(fmt, result)
  }
  invisible(packages)
}

#' Find the path to Python associated with a conda environment
#'
#' @inheritParams reticulate::conda_python
conda_python <- function(envname = NULL, conda = "auto", all = FALSE) {
  envname <- get_envname(envname)
  envname <- reticulate:::condaenv_resolve(envname)
  if (grepl("[/\\\\]", envname)) {
    suffix <- if (reticulate:::is_windows()) "python.exe" else "bin/python"
    path <- file.path(envname, suffix)
    if (file.exists(path)) {
      return(path)
    }
    fmt <- "no conda environment exists at path '%s'"
    stop(sprintf(fmt, envname))
  }
  conda_envs <- reticulate::conda_list(conda = conda)
  conda_envs <- conda_envs[grep(normalizePath(reticulate:::conda_info(conda = conda)$envs_dirs[1]), x = normalizePath(conda_envs$python), fixed = TRUE), , drop = FALSE]
  env <- conda_envs[conda_envs$name == envname, , drop = FALSE]
  if (nrow(env) == 0) {
    stop("conda environment \"", envname, "\" not found")
  }
  python <- if (all) env$python else env$python[[1L]]
  return(normalizePath(as.character(python)))
}

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
    stop("Failed to run \"", command, "\". Please check manually.")
  })
}

#' Download File from the Internet
#'
#' @inheritParams utils::download.file
#' @param methods Methods to be used for downloading files. The default is to try different download methods in turn until the download is successfully completed.
#' @param attempts Number of attempts for each download method.
#' @param ... Other arguments passed to \code{\link[utils]{download.file}}
#'
#' @importFrom utils download.file
#' @export
download <- function(url, destfile, methods = c("auto", "wget", "libcurl", "curl", "wininet", "internal"), quiet = FALSE, ..., attempts = 2) {
  if (missing(url) || missing(destfile)) {
    stop("'url' and 'destfile' must be both provided.")
  }
  attempts <- 0
  status <- NULL
  while (is.null(status)) {
    for (method in methods) {
      status <- tryCatch(expr = {
        suppressWarnings(download.file(url, destfile = destfile, method = method, quiet = quiet, ...))
        status <- 1
      }, error = function(error) {
        message(error)
        message("Cannot download from the url: ", url)
        message("Failed to download using \"", method, "\". Retry...")
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
  return(invisible(NULL))
}

kegg_get <- function(url) {
  temp <- tempfile()
  download(url = url, destfile = temp)
  content <- as.data.frame(do.call(rbind, strsplit(readLines(temp), split = "\t")))
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

#' Invoke a function with a list of arguments
#' @param .fn A function, or function name as a string.
#' @param .args A list of arguments.
#' @param Other arguments passed to the function.
#' @param .env Environment in which to evaluate the call. This will be most useful if .fn is a string, or the function has side-effects.
#' @importFrom rlang caller_env is_null is_scalar_character is_character is_function set_names env env_get env_bind syms call2
#' @export
invoke <- function(.fn, .args = list(), ..., .env = caller_env()) {
  args <- c(.args, list(...))
  .bury <- c(".fn", "")
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

#' Implement similar functions to the \code{unnest} function in the tidyr package
#' @param data A data frame.
#' @param cols Columns to unnest.
#' @param keep_empty By default, you get one row of output for each element of the list your unchopping/unnesting. This means that if there's a size-0 element (like \code{NULL} or an empty data frame), that entire row will be dropped from the output. If you want to preserve all rows, use \code{keep_empty = TRUE} to replace size-0 elements with a single row of missing values.
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

#' Attempts to turn a dgCMatrix into a dense matrix
#' @param matrix A dgCMatrix
#' @useDynLib SCP
#' @export
as_matrix <- function(matrix) {
  if (!inherits(matrix, "dgCMatrix")) {
    stop("matrix is not a dgCMatrix.")
  }
  row_pos <- matrix@i
  col_pos <- findInterval(seq(matrix@x) - 1, matrix@p[-1])

  out <- asMatrix(
    rp = row_pos, cp = col_pos, z = matrix@x,
    nrows = matrix@Dim[1], ncols = matrix@Dim[2]
  )

  row.names(out) <- matrix@Dimnames[[1]]
  colnames(out) <- matrix@Dimnames[[2]]
  return(out)
}

#' Capitalizes the characters
#' Making the first letter uppercase
#'
#' @param x A vector of character strings to be capitalized.
#' @param force_tolower Whether to force the remaining letters to be lowercase.
#' @export
capitalize <- function(x, force_tolower = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  if (!inherits(x, "character")) {
    stop("x must be the type of character.")
  }
  if (isTRUE(force_tolower)) {
    paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))), sep = "")
  } else {
    paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep = "")
  }
}

str_wrap <- function(x, width = 80) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  x_wrap <- unlist(lapply(x, function(i) paste0(strwrap(i, width = width), collapse = "\n")))
  return(x_wrap)
}

#' Split a vector into the chunks
#'
#' @param x A vector.
#' @param nchunks Number of chunks.
#' @examples
#' x <- 1:10
#' names(x) <- letters[1:10]
#' tochunks(x, nchunks = 3)
#' @export
tochunks <- function(x, nchunks) {
  split(x, cut(seq_along(x), nchunks, labels = FALSE))
}

#' Generate a iterator along chunks of a vector
#' @param x A vector.
#' @param nchunks Number of chunks.
#' @examples
#' \dontrun{
#' library(BiocParallel)
#' x <- 1:100
#' BPPARAM <- bpparam()
#' bpprogressbar(BPPARAM) <- TRUE
#' bpworkers(BPPARAM) <- 10
#' slow_fun <- function(x) {
#'   out <- NULL
#'   for (i in seq_along(x)) {
#'     Sys.sleep(0.5)
#'     out[[i]] <- x[[i]] + 3
#'   }
#'   return(out)
#' }
#' system.time({
#'   res0 <- lapply(x, FUN = slow_fun)
#' })
#' unlist(res0, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res1 <- bplapply(x, FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res1, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res2 <- bplapply(tochunks(x, nchunks = bpworkers(BPPARAM)), FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res2, recursive = FALSE, use.names = FALSE)[71:73]
#' system.time({
#'   res3 <- bpiterate(ITER = iterchunks(x, nchunks = bpworkers(BPPARAM)), FUN = slow_fun, BPPARAM = BPPARAM)
#' })
#' unlist(res3, recursive = FALSE, use.names = FALSE)[71:73]
#' }
#' @export
iterchunks <- function(x, nchunks) {
  chunks <- tochunks(x, nchunks)
  i <- 0L
  function() {
    if (i >= length(chunks)) {
      return(NULL)
    }
    i <<- i + 1L
    x[chunks[[i]]]
  }
}
