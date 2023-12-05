

library(reticulate)

.add_numbers_python <- function(a, b) {

  reticulate::source_python(system.file("python/my_python_script.py", package = "FORTLS"))

  reticulate::py$add_numbers(a, b)

}
