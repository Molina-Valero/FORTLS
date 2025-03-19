
# Downsampling

voxel_grid_downsampling <- function(points, voxel_size) {

  reticulate::source_python(system.file("python", "voxel_grid_downsampling.py", package = "FORTLS"))

  points <- reticulate::np_array(points)

  reticulate::py$voxel_grid_downsampling(points, voxel_size)

}


# Species classification

random.forest.fit <- function(dir.data, file_name = NULL) {

  reticulate::source_python(system.file("python", "random_forest_fit.py", package = "FORTLS"))

  data_path <- paste(dir.data, "/", sep = "")

  if(is.null(file_name)){file_name = "SpeciesClassification.csv"}

  reticulate::py$random_forest_fit(data_path, file_name)

}


random.forest.sp <- function(file_name) {

  reticulate::source_python(system.file("python", "random_forest_sp.py", package = "FORTLS"))

  # data_path <- paste(dir.data, "/", sep = "")

  # if(is.null(file_name)){file_name = "SpeciesClassification_tree.csv"}

  reticulate::py$random_forest_sp(file_name)

}


