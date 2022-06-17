.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)
  options(expressions = 5e5)
  if (Sys.getenv("SCP_PYTHON") != "") {
    packageStartupMessage("Use the environment variable 'SCP_PYTHON'.")
    Sys.setenv(RETICULATE_PYTHON = Sys.getenv("SCP_PYTHON"))
    pythons <- Sys.getenv("SCP_PYTHON")
  } else {
    packageStartupMessage("Default python will be used by SCP. \nOne can use Sys.setenv(SCP_PYTHON='/path/to/python') to specify the python version before library(SCP).")
    pythons <- unique(c(Sys.getenv("RETICULATE_PYTHON"), Sys.which("python3"), Sys.which("python")))
  }

  python_path <- NULL
  sys_bit <- ifelse(grepl("64", Sys.info()["machine"]), "64bit", "32bit")
  for (py in pythons) {
    if (!file.exists(py)) {
      next
    }
    py_version <- tryCatch(suppressWarnings(reticulate:::python_version(py)),
      error = identity
    )
    if (inherits(py_version, "error")) {
      next
    }
    if (py_version < numeric_version("3.7.0") || py_version > numeric_version("3.10.0")) {
      next
    }
    py_bit <- tryCatch(suppressWarnings(system2(command = py, args = " -c 'import platform; print(platform.architecture()[0])'", stdout = TRUE)),
      error = identity
    )
    if (inherits(py_bit, "error") || length(py_bit) == 0) {
      next
    }
    if (identical(py_bit, sys_bit)) {
      python_path <- py
      packageStartupMessage("python path: ", python_path)
      break
    } else {
      packageStartupMessage(py, "architecture is ", py_bit, " and system is ", sys_bit)
      next
    }
  }

  if (is.null(python_path)) {
    packageStartupMessage("Python is unavailable. Install python automatically ...")
    git_exist <- suppressWarnings(system("git", ignore.stdout = TRUE, ignore.stderr = TRUE))
    if (git_exist == 127) {
      warning("You need to install git first! (http://git-scm.com/download/)", immediate. = TRUE)
      return(invisible(NULL))
    }
    python_path <- reticulate::install_python(version = ifelse(sys_bit == "64bit", "3.8.7", "3.8.7-win32"))
  }

  if (!reticulate::virtualenv_exists("SCP")) {
    packageStartupMessage("Create SCP virtual environment. The path is: ", reticulate:::virtualenv_path("SCP"))
    reticulate::virtualenv_create(envname = "SCP", python = python_path)
  }
  Sys.setenv(RETICULATE_PYTHON = reticulate::virtualenv_python("SCP"))

  version <- tryCatch(suppressWarnings(reticulate:::python_version(python_path)), error = identity)
  if (inherits(py_version, "error")) {
    warning("SCP need python 3.7-3.9! Please install python and reload the SCP!", immediate. = TRUE)
    return(invisible(NULL))
  } else {
    if (version < numeric_version("3.7.0") || version > numeric_version("3.10.0")) {
      warning("SCP currently only support python version 3.7-3.9! The version of Python currently is ", version,
        "!\nPython related functions may not work. Please install the right version of python and reload the SCP!",
        immediate. = TRUE
      )
      return(invisible(NULL))
    } else {
      reticulate::use_virtualenv("SCP", required = TRUE)
      pyinfo <- utils::capture.output(reticulate::py_config())
      pyinfo_mesg <- c(
        "======================== SCP python config ========================",
        pyinfo,
        "==================================================================="
      )
      invisible(lapply(pyinfo_mesg, packageStartupMessage))
      if (!exist_pkg("pip")) {
        temp <- tempfile()
        download.file("https://bootstrap.pypa.io/get-pip.py", temp)
        suppressWarnings(system2(command = reticulate::virtualenv_python("SCP"), args = temp, stdout = TRUE))
        unlink(temp)
      }
      check_Python(pkgs = "matplotlib", envname = "SCP")
      run_Python(command = "import matplotlib", envir = .GlobalEnv)
      if (!interactive()) {
        run_Python(command = "matplotlib.use('pdf')", envir = .GlobalEnv)
      }
      run_Python(command = "import matplotlib.pyplot as plt", envir = .GlobalEnv)
      check_Python(pkgs = "versioned-hdf5", envname = "SCP")
      check_Python(pkgs = c("numba==0.53.1", "scanpy"), envname = "SCP")
      run_Python(command = "import scanpy", envir = .GlobalEnv)

      SCP_analysis <- reticulate::import_from_path("SCP_analysis", system.file("python", package = utils::packageName(), mustWork = TRUE))
      lapply(names(SCP_analysis), function(name) assign(name, SCP_analysis[[name]], envir = as.environment("package:SCP")))
    }
  }
}
