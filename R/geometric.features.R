

geometric.features <- function(data, dist = 0.05){


  # Select necessary fields from original txt file of point cloud

  data <- as.data.frame(data[, c("point", "x", "y", "z")])

  # Create x and y coordinates for grid

  x <- seq(min(data$x), max(data$x)+1)
  y <- seq(min(data$y), max(data$y)+1)


  if(length(x) < 2 | length(y) < 2){

    dat <- as.matrix(data[, 1:4])

    geo.fea <- geometric_features(dat, dist)

    } else{

      # Empty data frame where coordinates neccesaries for
      # creating grid will be saved

      grid <- data.frame(x = rep(x, each = length(y)),
                         y = rep(y, times = length(x)))

      grid <- sf::st_as_sf(grid, coords = c("x","y"))

      grid <- sf::st_buffer(grid, dist = 0.5 + dist, endCapStyle = "SQUARE")
      grid <- sf::st_cast(grid, "POLYGON")

      grid.2 <- sf::st_intersects(grid, sf::st_as_sf(data, coords = c("x","y")))
      grid.2 <- as.data.frame(grid.2)
      colnames(grid.2) <- c("id", "code")
      data$code <- as.numeric(row.names(data))
      data <- merge(data, grid.2, by = "code")
      data <- data[, 2:ncol(data)]

      rm(grid, grid.2)

      dat <- lapply(split(data[, 1:4], data$id), as.matrix)

      dat <- dat[(lapply(dat, nrow)) > 2]

      cl <- parallel::makeCluster(parallel::detectCores() - 1)

      # geo.fea <- do.call(rbind, lapply(dat, geometric_features, dist))

      geo.fea <- do.call(rbind, parallel::clusterApply(cl, dat, geometric_features, dist))

      parallel::stopCluster(cl)

      }


      geo.fea <- geo.fea[geo.fea$surface_variation > 0 & geo.fea$surface_variation < 9999, ]


      if(is.null(geo.fea)){

        data$first_eigenvalue <- NA
        data$second_eigenvalue <- NA
        data$third_eigenvalue <- NA
        data$sum_eigenvalues <- NA
        data$omnivariance <- NA
        data$eigenentropy <- NA
        data$PCA1 <- NA
        data$PCA2 <- NA
        data$anisotropy <- NA
        data$planarity <- NA
        data$linearity <- NA
        data$surface_variation <- NA
        data$sphericity <- NA
        data$verticality <- NA
        data$number_neighbors <- NA

        data <- data
        return(data)

        }

      data <- merge(data[, c("point", "x", "y", "z")], geo.fea, by = "point", all = TRUE)

      data <- data[!duplicated(data), ]

      data$omnivariance <- (data$first_eigenvalue * data$second_eigenvalue * data$third_eigenvalue) ^ (1 / 3)

      data$eigenentropy <- -(data$first_eigenvalue * log(data$first_eigenvalue) +
                             data$second_eigenvalue * log(data$second_eigenvalue) +
                             data$third_eigenvalue * log(data$third_eigenvalue))

      data$surface_density <- data$number_neighbors / (4 * pi * dist ^ 2)

      data$volume_density <- data$number_neighbors / ((4 / 3) * pi * dist ^ 3)


      data <- data[, c("point", "x", "y", "z",
                       "first_eigenvalue", "second_eigenvalue", "third_eigenvalue",
                       "sum_eigenvalues", "omnivariance", "eigenentropy",
                       "PCA1", "PCA2",
                       "anisotropy",  "planarity", "linearity", "surface_variation",  "sphericity", "verticality",
                       "number_neighbors", "surface_density", "volume_density")]

  return(data)

}
