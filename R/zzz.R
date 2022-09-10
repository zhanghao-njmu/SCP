.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)
  options(expressions = 5e5)
  if (reticulate::virtualenv_exists("SCP")) {
    Sys.setenv(RETICULATE_PYTHON = reticulate::virtualenv_python("SCP"))
    reticulate::use_virtualenv("SCP", required = TRUE)
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
