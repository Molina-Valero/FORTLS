#' Install Python dependencies required by FORTLS
#'
#' @param envname Name of the virtual environment or conda environment.
#' @param method Either "virtualenv", "conda", or "auto".
#' @export
install_fortls_python_deps <- function(envname = "fortls_env", method = "auto") {
  # Install Miniconda if using conda and it's not available
  if (method == "conda") {
    reticulate::install_miniconda()
  }

  reticulate::py_install(
    packages = c("pandas", "numpy", "jakteristics"),  # Add more as needed
    envname = envname,
    method = method,
    pip = TRUE
  )

  message("Python environment and packages for FORTLS successfully installed.")
}
