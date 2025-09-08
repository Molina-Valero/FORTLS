

.onLoad <- function(libname, pkgname) {
  # Attempt to use pre-defined environment
  if (reticulate::virtualenv_exists("fortls_env")) {
    reticulate::use_virtualenv("fortls_env", required = FALSE)
  } else if (reticulate::condaenv_exists("fortls_env")) {
    reticulate::use_condaenv("fortls_env", required = FALSE)
  }
}
