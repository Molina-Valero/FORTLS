
normalize <- function(las, normalized = NULL,
                      x.center = NULL, y.center = NULL,
                      x.side = NULL, y.side = NULL,
                      xpoly = NULL, ypoly = NULL,
                      max.dist = NULL, min.height = NULL, max.height = 50,
                      algorithm.dtm = "knnidw", res.dtm = 0.2,
                      csf = list(cloth_resolution = 0.5),
                      intensity = NULL, RGB = NULL,
                      scan.approach = "single",
                      voxel_size = NULL,
                      id = NULL, file = NULL, plot = TRUE,
                      dir.data = NULL, save.result = TRUE, dir.result = NULL,
                      save.las = NULL){


  set.seed(123)

  if(is.null(normalized)){
  .pb <- progress::progress_bar$new(total = 12)} else {
    .pb <- progress::progress_bar$new(total = 6)}

  .pb$tick()


  # Obtaining working directory for loading files
  if(is.null(dir.data))
    dir.data <- getwd()

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()

  # Reading input (LAS file)

  if(class(las)[1]=="LAS"){

    .las <- las

    } else {

  # Loading input (LAS file)

  # if(is.null(RGB) & is.null(intensity)){
  # .las <- suppressWarnings(suppressMessages(lidR::readLAS(file.path(dir.data, las), select = "xyz")))}
  #   else if (is.null(RGB) & !is.null(intensity)){
  #     .las <- suppressWarnings(suppressMessages(lidR::readLAS(file.path(dir.data, las), select = "xyzIntensity")))}
  #       else if (!is.null(RGB) & is.null(intensity)){
  #         .las <- suppressWarnings(suppressMessages(lidR::readLAS(file.path(dir.data, las), select = "xyzRGB")))}
  #           else{
  #             .las <- suppressWarnings(suppressMessages(lidR::readLAS(file.path(dir.data, las), select = "xyzIntensityRGB")))}
  #   }

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
    .las <- suppressWarnings(suppressMessages(
      lidR::readLAS(file.path(dir.data, las), select = select_option)
    ))

  }


  .las <- lidR::filter_duplicates(.las)


  .pb$tick()

  # Establishing center

  if(is.null(x.center)) {

    x.center <- mean(c(.las@header@PHB$`Max X`, .las@header@PHB$`Min X`))

  }

  if(is.null(y.center)) {

    y.center <- mean(c(.las@header@PHB$`Max Y`, .las@header@PHB$`Min Y`))

  }


  # Giving the same scale factor to all coordinates

  .las@header@PHB[["X scale factor"]] <- 0.001
  .las@header@PHB[["Y scale factor"]] <- 0.001
  .las@header@PHB[["Z scale factor"]] <- 0.001

  .pb$tick()

  # Data filtering at horizontal distances larger than max_dist m in the horizontal plane

  if(!is.null(normalized)) {

    if (!is.null(max.dist)) {

      .data <- lidR::clip_circle(.las, x.center, y.center, max.dist)

    } else if (!is.null(x.side) && !is.null(y.side)) {

      .data <- lidR::clip_rectangle(
        .las,
        x.center - (x.side / 2), y.center - (y.side / 2),
        x.center + (x.side / 2), y.center + (y.side / 2)
      )

    } else if (!is.null(xpoly) && !is.null(ypoly)) {

      .data <- lidR::clip_polygon(.las, xpoly, ypoly)

    } else {

      .data <- .las

    }

    # Plot if requested
    if (!is.null(plot)) lidR::plot(.data)

    # Save LAZ file if required
    if (!is.null(save.las)) {
      lidR::writeLAS(.data, paste(dir.result, "/", id, ".laz", sep = ""))
    }

    .data <- data.frame(.data@data)
    .data$slope <- 0
    .pb$tick()

  } else {

  # Normalize

  .data <- suppressWarnings(suppressMessages(lidR::classify_ground(.las, algorithm = lidR::csf(cloth_resolution = csf$cloth_resolution), last_returns = FALSE)))

  .pb$tick()


  # Generaion of Digital Terrain Model (DTM)


  if(algorithm.dtm == "knnidw")
    .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::knnidw())))


  if(algorithm.dtm == "tin")
    .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::tin())))


  .dtm[.dtm < min(.data@data$Z)] <- NA

  .pb$tick()

  # Estimating slope

  raster::crs(.dtm) <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83"
  .slope <- raster::terrain(.dtm, opt=c('slope'), unit='radians')

  .pb$tick()

  if(mean(.slope@data@values, na.rm = TRUE) > 0.5){


    # Normalize

    .data <- suppressWarnings(suppressMessages(lidR::classify_ground(.las, algorithm = lidR::csf(sloop_smooth = TRUE), last_returns = FALSE)))


    # Generaion of Digital Terrain Model (DTM)

    if(algorithm.dtm == "knnidw")
      .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::knnidw())))


    if(algorithm.dtm == "tin")
      .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::tin())))


    .dtm[.dtm < min(.data@data$Z)] <- NA


  }


  rm(.las)

  .pb$tick()


  # Normalization of cloud data

  .data <- suppressWarnings(suppressMessages(lidR::normalize_height(.data, .dtm, add_lasattribute = FALSE, na.rm = TRUE)))

  rm(.dtm)

  .pb$tick()



  # Data filtering at horizontal distances larger than max_dist m in the horizontal plane

  if(!is.null(max.dist)){

    .data <- lidR::clip_circle(.data, x.center, y.center, max.dist)

    } else if (!is.null(x.side) | !is.null(y.side)){

    .data <- lidR::clip_rectangle(.data, x.center - (x.side / 2), y.center - (y.side / 2),
                                  x.center + (x.side / 2), y.center + (y.side / 2))

    } else if (!is.null(xpoly) | !is.null(ypoly)){

      .data <- lidR::clip_polygon(.las, xpoly, ypoly)

  }

  .pb$tick()


  # Plot

  if(!is.null(plot))
    lidR::plot(.data)



  # Saving laz file

  if(!is.null(save.las))
    lidR::writeLAS(.data, paste(dir.result, "/", id, ".laz", sep = ""))



  # Assigning slope to point cloud

  .data <- lidR::merge_spatial(.data, .slope, "slope")

  rm(.slope)

  .pb$tick()

  # Test for getting a smaller file - data frame

  .data <- data.frame(.data@data)
  .data <- data.table::setDT(.data)



  }

  # Removing points classified as ground

  .data <- subset(.data, .data$Classification == 1)



  # Extracting coordinates values

  if(is.null(RGB) & is.null(intensity)){
    .data <- .data[, c("X", "Y", "Z", "slope"), drop = FALSE]
    colnames(.data) <- c("x", "y", "z", "slope")}
  else if (is.null(RGB) & !is.null(intensity)){
    .data <- .data[, c("X", "Y", "Z", "slope", "Intensity"), drop = FALSE]
    colnames(.data) <- c("x", "y", "z", "slope", "intensity")}
  else if (!is.null(RGB) & is.null(intensity)){
    .data <- .data[, c("X", "Y", "Z", "slope", "R", "G", "B"), drop = FALSE]
    colnames(.data) <- c("x", "y", "z", "slope", "R", "G", "B")}
  else{
    .data <- .data[, c("X", "Y", "Z", "slope", "Intensity", "R", "G", "B"), drop = FALSE]
    colnames(.data) <- c("x", "y", "z", "slope", "intensity", "R", "G", "B")}


  # Low point filtering

  if(!is.null(min.height)) {

    .data <- .data[which(.data$z > min.height), , drop = FALSE]

  }


  # High point filtering

  if(!is.null(max.height)) {

    .data <- .data[which(.data$z < max.height), , drop = FALSE]

  }

  # Transformation to other coordinate systems

  # Cylindrical coordinate system (https://en.wikipedia.org/wiki/Cylindrical_coordinate_system)
  # rho, axial distance or radial distance (euclidean distance from the z-axis to the point P)
  # phi, azimuth is the angle between the reference direction on the chosen plane and the line from the origin to the projection of P on the plane
  # z, axial coordinate or height z is the signed distance from the chosen plane to the point P
  .data$rho <- sqrt((.data$x - x.center) ^ 2 + (.data$y - y.center) ^ 2)
  .data$phi <- atan2(.data$y - y.center, .data$x - x.center)
  .data$phi <- ifelse(.data$phi < 0, .data$phi + (2 * pi), .data$phi)

  # Spherical coordinates system (https://en.wikipedia.org/wiki/Spherical_coordinate_system)
  # r, radius or radial distance is the Euclidean distance from the origin O to P
  # theta, inclination (or polar angle) is the angle between the zenith direction and the line segment OP
  # phi, azimuth is the angle between the reference direction on the chosen plane and the line from the origin to the projection of P on the plane

  .data$r <- sqrt(.data$z ^ 2 + .data$rho ^ 2)
  .data$theta <- atan2(.data$z, .data$rho)

  .data$point <- as.integer((1:nrow(.data)))


  # Green Leaf Algorithm (GLA) (Louhaichi et al., (2001))
  if(!is.null(RGB))
    .data$GLA <- (2 * .data$G - .data$R - .data$B) / (2 * .data$G + .data$R + .data$B)


  # Point crooping process (TLS single-scan)
  # This is a previous step to obtain a homogeneous density of points in the space
  # This is based on the principle that closer objects (with the same size and shape)
  # have more probability to receive points

  if(scan.approach == "single"){

    .data$prob <- (.data$r / max(.data$r)) ^ 2
    .data$prob.random <- stats::runif(nrow(.data))
    .data$prob.selec <- as.integer(ifelse(.data$prob > .data$prob.random, 1, 0))}


  # For the rest of situations, point cloud is downsampled by voxelization

  # if(scan.approach == "multi" & is.null(voxel_size))
  #   voxel_size <- 0.01


  if(scan.approach == "multi" & !is.null(voxel_size)){

    if(is.null(RGB) & is.null(intensity)){
      .data <- as.data.frame(voxel_grid_downsampling(as.matrix(.data[, c("x", "y", "z", "slope")]), voxel_size))
      colnames(.data) <- c("x", "y", "z", "slope")}

    else if (is.null(RGB) & !is.null(intensity)){
      .data <- as.data.frame(voxel_grid_downsampling(as.matrix(.data[, c("x", "y", "z", "slope", "intensity")]), voxel_size))
      colnames(.data) <- c("x", "y", "z", "slope", "intensity")}

    else if (!is.null(RGB) & is.null(intensity)){
      .data <- as.data.frame(voxel_grid_downsampling(as.matrix(.data[, c("x", "y", "z", "slope", "R", "G", "B", "GLA")]), voxel_size))
      colnames(.data) <- c("x", "y", "z", "slope", "R", "G", "B")}

    else{
      .data <- as.data.frame(voxel_grid_downsampling(as.matrix(.data[, c("x", "y", "z", "slope", "intensity", "R", "G", "B", "GLA")]), voxel_size))
      colnames(.data) <- c("x", "y", "z", "slope", "intensity", "R", "G", "B")}


    .data$rho <- sqrt((.data$x - x.center) ^ 2 + (.data$y - y.center) ^ 2)
    .data$phi <- atan2(.data$x - x.center, .data$y - y.center)
    .data$phi <- ifelse(.data$phi < 0, .data$phi + (2 * pi), .data$phi)
    .data$r <- sqrt(.data$z ^ 2 + .data$rho ^ 2)
    .data$theta <- atan2(.data$z, .data$rho)

    .data$point <- as.integer((1:nrow(.data)))}


    if(scan.approach == "multi"){

    .data$prob <- stats::runif(nrow(.data))
    .data$prob.selec <- as.integer(ifelse(.data$prob > 0.5, 1, 0))}


  # Assign id

  if(!is.null(id)){

    .data$id <- id

  } else {

    .data$id <- as.integer(1)

  }


  # File name

  if(!is.null(file)){

    .data$file <- file

  } else {

    .data$file <- paste(.data$id[1], ".txt", sep = "")

  }

  if(is.null(RGB) & is.null(intensity)){
    .data <- .data[, c("id", "file", "point", "x", "y", "z", "rho", "phi", "r", "theta", "slope", "prob", "prob.selec"), drop = FALSE]}
  else if (is.null(RGB) & !is.null(intensity)){
    .data <- .data[, c("id", "file", "point", "x", "y", "z", "rho", "phi", "r", "theta", "slope", "intensity", "prob", "prob.selec"), drop = FALSE]}
  else if (!is.null(RGB) & is.null(intensity)){
    .data <- .data[, c("id", "file", "point", "x", "y", "z", "rho", "phi", "r", "theta", "slope", "R", "G", "B", "GLA", "prob", "prob.selec"), drop = FALSE]}
  else{.data <- .data[, c("id", "file", "point", "x", "y", "z", "rho", "phi", "r", "theta", "slope", "intensity", "R", "G", "B", "GLA", "prob", "prob.selec"), drop = FALSE]}

  .pb$tick()

  # Saving data

  # Obtaining working directory for saving files

  if(isTRUE(save.result)){

    .data.red <- .data[which(.data$prob.selec == 1), , drop = FALSE]

    vroom::vroom_write(.data.red, path = file.path(dir.result, .data.red$file[1]), delim = ",", progress = FALSE)

  }

  .pb$tick()

  return(.data)

}
