.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)
  options(expressions = 5e5)
  conda_exist <- reticulate:::conda_installed()
  if (conda_exist) {
    env_exist <- file.exists(paste0(reticulate:::conda_info()$conda_prefix, "/envs/SCP"))
  } else {
    env_exist <- FALSE
  }
  if (env_exist && (is.null(getOption("SCP_env_init")) || getOption("SCP_env_init") != FALSE)) {
    try({
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
      packageStartupMessage("SCP python environment can be disabled with the command 'options(SCP_env_init = FALSE)' before loading the package")
    })
  }
}
