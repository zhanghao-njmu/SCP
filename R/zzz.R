.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)
  options(expressions = 5e5)
  env_exist <- isTRUE(tryCatch("SCP" %in% reticulate::conda_list()$name[grep(reticulate:::conda_info()$conda_prefix, reticulate::conda_list()$python)], error = identity))
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
