
normalize_optimized <- function(las, normalized = NULL,
                      x.center = NULL, y.center = NULL,
                      x.side = NULL, y.side = NULL,
                      xpoly = NULL, ypoly = NULL,
                      max.dist = NULL, min.height = NULL, max.height = 50,
                      algorithm.dtm = "knnidw", res.dtm = 0.2,
                      csf = list(cloth_resolution = 0.5),
                      intensity = NULL, RGB = NULL,
                      scan.approach = "single",
                      voxel_size = NULL,
                      threads = 1,
                      id = NULL, file = NULL, plot = TRUE,
                      dir.data = NULL, save.result = TRUE, dir.result = NULL,
                      save.las = NULL){


  threads <- max(1, threads)
  lidR::set_lidr_threads(threads)


  if(is.null(normalized)){
  pb <- progress::progress_bar$new(total = 12)} else {
    pb <- progress::progress_bar$new(total = 6)}

  pb$tick()


  # Obtaining working directory for loading files
  if(is.null(dir.data))
    dir.data <- getwd()

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()

  # Reading input (LAS file)

  if(inherits(las, "LAS")){

    las <- las

    } else {

    # Initialize the select option with a default value
    select_option <- "xyz"

    # Modify select_option based on RGB and intensity presence
    if (is.null(RGB) & is.null(intensity)) {

      select_option <- "xyz"

    } else if (is.null(RGB) & !is.null(intensity)) {

      select_option <- "xyzIntensity"

    } else if (!is.null(RGB) & is.null(intensity)) {

      select_option <- "xyzRGB"

    } else {

      select_option <- "xyzIntensityRGB"

    }

    # If normalized = TRUE, append "Classification"
    if (!is.null(normalized)) {
      select_option <- paste(select_option, "Classification")
    }

    # Read LAS file with the selected attributes
    las <- suppressWarnings(suppressMessages(
      lidR::readLAS(file.path(dir.data, las), select = select_option)
    ))

  }

  # Filter duplicates and clean up immediately
  las <- lidR::filter_duplicates(las)


  pb$tick()

  # Establishing center

  if(is.null(x.center)) {

    x.center <- mean(c(las@header@PHB$`Max X`, las@header@PHB$`Min X`))

  }

  if(is.null(y.center)) {

    y.center <- mean(c(las@header@PHB$`Max Y`, las@header@PHB$`Min Y`))

  }


  # Giving the same scale factor to all coordinates

  las@header@PHB[["X scale factor"]] <- 0.001
  las@header@PHB[["Y scale factor"]] <- 0.001
  las@header@PHB[["Z scale factor"]] <- 0.001

  pb$tick()

  # Data filtering at horizontal distances larger than max_dist m in the horizontal plane

  if(!is.null(normalized)) {

    if (!is.null(max.dist)) {

      data <- lidR::clip_circle(las, x.center, y.center, max.dist)

    } else if (!is.null(x.side) && !is.null(y.side)) {

      data <- lidR::clip_rectangle(
        las,
        x.center - (x.side / 2), y.center - (y.side / 2),
        x.center + (x.side / 2), y.center + (y.side / 2)
      )

    } else if (!is.null(xpoly) && !is.null(ypoly)) {

      data <- lidR::clip_polygon(las, xpoly, ypoly)

    } else {

      data <- las

    }

    # Plot if requested
    if (!is.null(plot)) lidR::plot(data)


    # Save LAZ file if required
    if (!is.null(save.las)) {
      lidR::writeLAS(data, paste(dir.result, "/", id, ".laz", sep = ""))
    }

    data <- data.frame(data@data)
    data$slope <- 0
    pb$tick()

  } else {

  # Normalize

  data <- suppressWarnings(suppressMessages(lidR::classify_ground(las, algorithm = lidR::csf(cloth_resolution = csf$cloth_resolution), last_returns = FALSE)))

  pb$tick()


  # Generaion of Digital Terrain Model (DTM)


  if(algorithm.dtm == "knnidw")
    dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(data, res = res.dtm, algorithm = lidR::knnidw())))


  if(algorithm.dtm == "tin")
    dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(data, res = res.dtm, algorithm = lidR::tin())))


  dtm[dtm < min(data@data$Z)] <- NA

  pb$tick()

  # Estimating slope

  if(is.na(sf::st_crs(las))){

    raster::crs(dtm) <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83"

  } else {

    raster::crs(dtm) <- sf::st_crs(las)

  }

  slope <- raster::terrain(dtm, opt=c('slope'), unit='radians')

  pb$tick()

  if(mean(slope@data@values, na.rm = TRUE) > 0.25){


    # Normalize

    data <- suppressWarnings(suppressMessages(lidR::classify_ground(las, algorithm = lidR::csf(sloop_smooth = TRUE), last_returns = FALSE)))


    # Generaion of Digital Terrain Model (DTM)

    if(algorithm.dtm == "knnidw")
      dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(data, res = res.dtm, algorithm = lidR::knnidw())))


    if(algorithm.dtm == "tin")
      dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(data, res = res.dtm, algorithm = lidR::tin())))


    dtm[dtm < min(data@data$Z)] <- NA


  }


  rm(las)

  pb$tick()


  # Normalization of cloud data

  if(algorithm.dtm == "knnidw"){

      # data <- suppressWarnings(suppressMessages(lidR::normalize_height(data, dtm, add_lasattribute = FALSE, na.rm = TRUE)))
      data <- suppressWarnings(suppressMessages(lidR::normalize_height(data, algorithm = lidR::knnidw(), use_class = c(2L, 9L), dtm = NULL)))

  }

  if(algorithm.dtm == "tin"){

    # data <- suppressWarnings(suppressMessages(lidR::normalize_height(data, dtm, add_lasattribute = FALSE, na.rm = TRUE)))
    data <- suppressWarnings(suppressMessages(lidR::normalize_height(data, algorithm = lidR::tin(), use_class = c(2L, 9L), dtm = NULL)))

  }


  data <- lidR::filter_poi(data, Z >= 0)

  rm(dtm)

  pb$tick()



  # Data filtering at horizontal distances larger than max_dist m in the horizontal plane

  if(!is.null(max.dist)){

    data <- lidR::clip_circle(data, x.center, y.center, max.dist)

    } else if (!is.null(x.side) | !is.null(y.side)){

    data <- lidR::clip_rectangle(data, x.center - (x.side / 2), y.center - (y.side / 2),
                                  x.center + (x.side / 2), y.center + (y.side / 2))

    } else if (!is.null(xpoly) | !is.null(ypoly)){

      data <- lidR::clip_polygon(las, xpoly, ypoly)

  }

  pb$tick()


  # Plot

  if (!is.null(plot)) lidR::plot(data)



  # Saving laz file

  if(!is.null(save.las))
    lidR::writeLAS(data, paste(dir.result, "/", id, ".laz", sep = ""))



  # Assigning slope to point cloud

  data <- lidR::merge_spatial(data, slope, "slope")

  rm(slope)

  pb$tick()

  # Test for getting a smaller file - data frame

  data <- data.frame(data@data)
  data <- data.table::setDT(data)



  }

  # Removing points classified as ground - OPTIMIZED with data.table syntax
  data <- data[Classification == 1]



  # Extracting coordinates values - OPTIMIZED to reduce redundancy
  cols <- c("X", "Y", "Z", "slope")
  col_names <- c("x", "y", "z", "slope")

  if (!is.null(intensity)) {
    cols <- c(cols, "Intensity")
    col_names <- c(col_names, "intensity")
  }

  if (!is.null(RGB)) {
    cols <- c(cols, "R", "G", "B")
    col_names <- c(col_names, "R", "G", "B")
  }

  data <- data[, cols, with = FALSE]
  data.table::setnames(data, old = cols, new = col_names)


  # Low and high point filtering - OPTIMIZED with data.table chaining
  if(!is.null(min.height)) {
    data <- data[z > min.height]
  }

  if(!is.null(max.height)) {
    data <- data[z < max.height]
  }

  # Transformation to other coordinate systems - OPTIMIZED
  # Cylindrical coordinate system
  # Pre-compute differences to avoid redundant calculations
  dx <- data$x - x.center
  dy <- data$y - y.center

  data[, rho := sqrt(dx^2 + dy^2)]
  data[, phi := atan2(dy, dx)]
  data[phi < 0, phi := phi + 2 * pi]

  # Spherical coordinates system
  data[, r := sqrt(z^2 + rho^2)]
  data[, theta := atan2(z, rho)]

  data[, point := .I]


  # Green Leaf Algorithm (GLA) (Louhaichi et al., (2001))
  if(!is.null(RGB))
    data[, GLA := (2 * G - R - B) / (2 * G + R + B)]


  # Point crooping process (TLS single-scan)

  if(scan.approach == "single"){

    data[, prob := (r / max(r))^2]
    data[, prob.random := stats::runif(.N)]
    data[, prob.selec := as.integer(prob > prob.random)]
  }


  # For the rest of situations, point cloud is downsampled by voxelization

  if(scan.approach == "multi" & !is.null(voxel_size)){

    # Determine columns for voxelization
    voxel_cols <- c("x", "y", "z", "slope")
    voxel_names <- c("x", "y", "z", "slope")

    if (!is.null(intensity)) {
      voxel_cols <- c(voxel_cols, "intensity")
      voxel_names <- c(voxel_names, "intensity")
    }

    if (!is.null(RGB)) {
      voxel_cols <- c(voxel_cols, "R", "G", "B")
      voxel_names <- c(voxel_names, "R", "G", "B")
    }

    data <- as.data.frame(voxel_grid_downsampling(as.matrix(data[, voxel_cols, with = FALSE]), voxel_size))
    colnames(data) <- voxel_names
    data <- data.table::setDT(data)

    # Recalculate coordinates - OPTIMIZED (fixed bug: swapped arguments)
    dx <- data$x - x.center
    dy <- data$y - y.center

    data[, rho := sqrt(dx^2 + dy^2)]
    data[, phi := atan2(dy, dx)]  # BUGFIX: was atan2(dx, dy)
    data[phi < 0, phi := phi + 2 * pi]
    data[, r := sqrt(z^2 + rho^2)]
    data[, theta := atan2(z, rho)]

    data[, point := .I]
  }


  if(scan.approach == "multi"){
    data[, prob := stats::runif(.N)]
    data[, prob.selec := as.integer(prob > 0.5)]
  }


  # Assign id

  if(!is.null(id)){
    data[, id := id]
  } else {
    data[, id := 1L]
  }


  # File name

  if(!is.null(file)){
    data[, file := file]
  } else {
    data[, file := paste0(id[1], ".txt")]
  }

  # Reorder columns - OPTIMIZED to reduce redundancy
  final_cols <- c("id", "file", "point", "x", "y", "z", "rho", "phi", "r", "theta", "slope")

  if (!is.null(intensity)) {
    final_cols <- c(final_cols, "intensity")
  }

  if (!is.null(RGB)) {
    final_cols <- c(final_cols, "R", "G", "B", "GLA")
  }

  final_cols <- c(final_cols, "prob", "prob.selec")

  data <- data[, final_cols, with = FALSE]

  pb$tick()


  # Adding the point of the plot center - OPTIMIZED direct creation

  new_row <- data.table::data.table(
    id = data$id[1],
    file = data$file[1],
    point = 0L,
    x = x.center,
    y = y.center,
    z = 0,
    rho = 0,
    phi = 0,
    r = 0,
    theta = 0,
    slope = NA_real_,
    prob = 1,
    prob.selec = 1L
  )

  # Add optional columns
  if (!is.null(intensity)) {
    new_row[, intensity := NA_real_]
  }

  if (!is.null(RGB)) {
    new_row[, c("R", "G", "B", "GLA") := list(NA_real_, NA_real_, NA_real_, NA_real_)]
  }

  # Reorder to match data columns
  data.table::setcolorder(new_row, names(data))

  # Bind it to the top of the existing data
  data <- rbind(new_row, data)


  # Saving data

  if(isTRUE(save.result)){

    data.red <- data[prob.selec == 1]

    vroom::vroom_write(data.red, file = file.path(dir.result, data.red$file[1]), delim = ",", progress = FALSE)

  }

  pb$tick()

  return(data)

}
