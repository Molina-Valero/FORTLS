


# Downsampling

voxel_grid_downsampling <- function(points, voxel_size) {

  reticulate::source_python(system.file("python", "voxel_grid_downsampling.py", package = "FORTLS"))

  points <- reticulate::np_array(points)

  reticulate::py$voxel_grid_downsampling(points, voxel_size)

}



# Geometric features

geometric_features_py <- function(data, dist = 0.1, geom_features_input = NULL, threads = 1) {

  reticulate::source_python(system.file("python", "geometric_features.py", package = "FORTLS"))

  if(is.null(geom_features_input)){
    geom_features_input <- c("verticality", "surface_variation", "planarity")}

  reticulate::py$compute_geom_features(input_df = data,
                                       dist = dist,
                                       features = geom_features_input,
                                       threads = as.integer(threads))

}




