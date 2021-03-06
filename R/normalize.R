
normalize <- function(las,
                      max.dist = NULL, min.height = NULL, max.height = NULL,
                      algorithm.dtm = "tin", res.dtm = 0.2,
                      id = NULL, file=NULL,
                      dir.data = NULL, save.result = TRUE, dir.result = NULL){

  .pb <- progress::progress_bar$new(total = 11)
  .pb$tick()


  # Obtaining working directory for loading files
  if(is.null(dir.data))
    dir.data <- getwd()

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()


  # Loading input (LAS file)

  .las <- suppressWarnings(suppressMessages(lidR::readLAS(file.path(dir.data, las), select = "xyz")))

  .pb$tick()


  # Giving the same scale factor to all coordinates

  .las@header@PHB[["X scale factor"]] <- 0.001
  .las@header@PHB[["Y scale factor"]] <- 0.001
  .las@header@PHB[["Z scale factor"]] <- 0.001


  # Data filtering at horizontal distances larger than max_dist m in the horizontal plane

  if(!is.null(max.dist)) {

    .las <- lidR::clip_circle(.las, 0, 0, max.dist)

  }

  .pb$tick()


  # Normalize

  # .ws  <- seq(3,12, 4)
  # .th  <- seq(0.1, 1.5, length.out = length(.ws))
  # .data <- lidR::classify_ground(.las, algorithm = lidR::pmf(.ws, .th), last_returns = FALSE)

  .data <- suppressWarnings(suppressMessages(lidR::classify_ground(.las, algorithm = lidR::csf(), last_returns = FALSE)))

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

  if(mean(.slope@data@values, na.rm = TRUE) > 0.3){


    # Normalize

    .data <- suppressWarnings(suppressMessages(lidR::classify_ground(.las, algorithm = lidR::csf(sloop_smooth = TRUE), last_returns = FALSE)))


    # Generaion of Digital Terrain Model (DTM)

    if(algorithm.dtm == "knnidw")
      .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::knnidw())))


    if(algorithm.dtm == "tin")
      .dtm <- suppressWarnings(suppressMessages(lidR::grid_terrain(.data, res = res.dtm, algorithm = lidR::tin())))


    .dtm[.dtm < min(.data@data$Z)] <- NA


  }

  .pb$tick()


  # Normalization of cloud data

  .data <- suppressWarnings(suppressMessages(lidR::normalize_height(.data, .dtm, add_lasattribute = FALSE, na.rm = TRUE)))

  .pb$tick()

  # Assigning slope to point cloud

  .data <- lidR::merge_spatial(.data, .slope, "slope")

  .pb$tick()

  # Test for getting a smaller file - data frame

  .data <- data.frame(.data@data)


  # Removing points classified as ground

  .data <- subset(.data, .data$Classification == 1)


  # Extracting coordinates values

  .data <- .data[, c("X", "Y", "Z", "slope"), drop = FALSE]
  colnames(.data) <- c("x", "y", "z", "slope")


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
  .data$rho <- sqrt(.data$x ^ 2 + .data$y ^ 2)
  .data$phi <- atan2(.data$y, .data$x)
  .data$phi <- ifelse(.data$phi < 0, .data$phi + (2 * pi), .data$phi)

  # Spherical coordinates system (https://en.wikipedia.org/wiki/Spherical_coordinate_system)
  # r, radius or radial distance is the Euclidean distance from the origin O to P
  # theta, inclination (or polar angle) is the angle between the zenith direction and the line segment OP
  # phi, azimuth is the angle between the reference direction on the chosen plane and the line from the origin to the projection of P on the plane

  .data$r <- sqrt(.data$z ^ 2 + .data$rho ^ 2)
  .data$theta <- atan2(.data$z, .data$rho)

  .data$point <- (1:nrow(.data))


  # Point crooping process
  # This is a previous step to obtain a homogeneous density of points in the space
  # This is based on the principle that closer objects (with the same size and shape)
  # have more probability to recieve points

  .data$prob <- (.data$r / max(.data$r)) ^ 2
  .data$prob.random <- stats::runif(nrow(.data))
  .data$prob.selec <- ifelse(.data$prob >= .data$prob.random, 1, 0)


  # Assign id

  if(!is.null(id)){

    .data$id <- id

  } else {

    .data$id <- 1

  }


  # File name

  if(!is.null(file)){

    .data$file <- file

  } else {

    .data$file <- paste(.data$id[1], ".txt", sep = "")

  }


  .data <- .data[, c("id", "file", "point", "x", "y", "z", "rho", "phi", "r", "theta", "slope", "prob", "prob.selec"), drop = FALSE]

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
