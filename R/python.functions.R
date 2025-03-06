# Define the path to the folder containing the Python scripts relative to the package directory
python_folder <- system.file("python", package = "FORTLS")


voxel_grid_downsampling <- function(points, voxel_size) {

  reticulate::source_python(file.path(python_folder, "voxel_grid_downsampling.py"))

  points <- reticulate::np_array(points)

  reticulate::py$voxel_grid_downsampling(points, voxel_size)

}

# Species classification

random.forest.fit <- function(dir.data, file_name = NULL) {

  reticulate::source_python(file.path(python_folder, "random_forest_fit.py"))

  data_path <- paste(dir.data, "/", sep = "")

  if(is.null(file_name)){file_name = "SpeciesClassification.csv"}

  reticulate::py$random_forest_fit(data_path, file_name)

}


random.forest.sp <- function(file_name) {

  reticulate::source_python(file.path(python_folder, "random_forest_sp.py"))

  # data_path <- paste(dir.data, "/", sep = "")

  # if(is.null(file_name)){file_name = "SpeciesClassification_tree.csv"}

  reticulate::py$random_forest_sp(file_name)

}


