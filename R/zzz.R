.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)
  options(expressions = 5e5)
  conda <- find_conda()
  if (!is.null(conda)) {
    env_exist <- file.exists(paste0(reticulate:::conda_info(conda = conda)$conda_prefix, "/envs/SCP"))
  } else {
    env_exist <- FALSE
    packageStartupMessage("Conda and SCP environment not found.\nIf you have already created an SCP python environment using conda, you can specify the conda path by setting options(reticulate.conda_binary = '/path/to/conda') before loading the package.")
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
