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
#' @export
check_R <- function(pkgs, pkg_names = NULL, install_methods = c("BiocManager::install", "install.packages", "devtools::install_github"), lib = .libPaths()[1], force = FALSE) {
  if (length(pkg_names) != 0 && length(pkg_names) != length(pkgs)) {
    stop("pkg_names must be NULL or a vector of the same length with pkgs")
  }
  status_list <- list()
  for (n in seq_along(pkgs)) {
    pkg <- pkgs[n]
    pkg_info <- strsplit(pkg, split = "/|@")[[1]]
    if (length(pkg_info) > 1) {
      pkg_name <- pkg_names[n] %||% pkg_info[2]
    } else {
      pkg_name <- pkg_names[n] %||% pkg_info
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
            eval(str2lang(paste0(install_methods[i], "('", pkg, "', lib='", lib, "', update = FALSE, ask = FALSE, force = TRUE)")))
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
  if (length(pkg_exist) == 0) {
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
check_Python <- function(pkgs, pkg_names = NULL, envname = "SCP", force = FALSE) {
  if (length(pkg_names) != 0 && length(pkg_names) != length(pkgs)) {
    stop("pkg_names must be NULL or a vector of the same length with pkgs")
  }
  status_list <- list()
  for (n in seq_along(pkgs)) {
    pkg <- pkgs[n]
    pkg_name <- pkg_names[n] %||% gsub("(.*)(\\[.*\\])|(==.*)", "\\1", pkg)
    exist <- exist_pkg(pkg_name, envname)
    if (!isTRUE(exist) || isTRUE(force)) {
      message("Try to install '", pkg, "' ...")
      tryCatch(expr = {
        reticulate::py_install(pkg, envname = envname)
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
