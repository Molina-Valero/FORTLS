
tree.detection <- function(data, dbh.min = 7.5, dbh.max = 200, ncr.threshold = 0.1, tls.resolution = list(), breaks=c(1.0, 1.3, 1.6), plot.attributes = NULL,
                           save.result = TRUE, dir.result = NULL){

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()

  # Converting arguments to International System of Units
  # Arguments of the forest inventory

  .dbh.min <- dbh.min / 100
  .dbh.max <- dbh.max / 100

  # Arguments of the TLS precision
  # If resolution is defined by points distance at a certain distance from TLS (mm/m):
  .point.dist <- tls.resolution$point.dist / 1000
  .tls.dist <- tls.resolution$tls.dist

  # If resolution is defined by angles (?):
  .vertical.angle <- tls.resolution$vertical.angle * 2 * pi / 360
  .horizontal.angle <- tls.resolution$horizontal.angle * 2 * pi / 360

  # Define angle aperture according to available TLS resolution parameters:

  # Angular resolution:
  if(is.null(.point.dist)){

    .alpha.v <- .vertical.angle
    .alpha.h <- .horizontal.angle

  } else {

    .alpha.v <- atan((.point.dist / 2) / .tls.dist) * 2
    .alpha.h <- .alpha.v

  }

  # Generation of homogenized point cloud

  # .ncr.threshold <- .ncr.threshold.double(data)

  .filteraux <- data.frame(cluster = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(),
                           num.points = as.numeric(), num.points.hom = as.numeric(),
                           phi.left = as.numeric(), phi.right = as.numeric(),
                           arc.circ = as.numeric(), sec = as.numeric())


  for(cuts in breaks){

    message("Computing section: ", cuts, " m")

    .cut <- data[which(data$z > (cuts-0.1) & data$z < (cuts+0.1)), , drop = FALSE]

    .cut <- .ncr.remove.slice.double(.cut)

    .cut <- .cut[which(.cut$ncr < ncr.threshold | is.na(.cut$ncr)), , drop = FALSE]

    # Restrict to slice corresponding to cuts m +/- 5 cm
    .cut <- .cut[which(.cut$z > (cuts-0.05) & .cut$z < (cuts+0.05)), , drop = FALSE]

    # Dbscan parameters
    # .eps <- .dbh.min / 2
    .eps <- (tan(.alpha.h / 2) * (max(.cut$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2)

    # Clustering
    .dbscan <- dbscan::dbscan(.cut[, c("x", "y"), drop = FALSE], eps = .eps)
    .cut$cluster <- .dbscan$cluster
    .cut <- .cut[which(.cut$cluster > 0), , drop = FALSE]

    # Checking if there are clusters
    if(nrow(.cut) < 1)
      next

    .cut$sec <- cuts

    .pb <- progress::progress_bar$new(total = length(unique(.cut$cluster)))

    .filter <- data.frame(cluster = as.numeric(),

                          # Center coordinates
                          center.x = as.numeric(), center.y = as.numeric(),
                          center.phi = as.numeric(), center.rho = as.numeric(),
                          center.r = as.numeric(), center.theta = as.numeric(),

                          # Radius
                          radius = as.numeric(),

                          # Number of points belonging to cluster (craw and after point cropping)
                          num.points = as.numeric(), num.points.hom = as.numeric(),

                          # Phi coordinates of left and right
                          phi.left = as.numeric(), phi.right = as.numeric(),

                          # Circumference arc
                          arc.circ = as.numeric(),

                          # Partial occlusion
                          occlusion = as.numeric())


    for(.i in unique(.cut$cluster)){

      .pb$tick()

      # Select cluster .i
      .dat <- .cut[which(.cut$cluster == .i), , drop = FALSE]

      # First filter
      .n <- (0.1 / (tan(.alpha.v / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2))

      if(nrow(.dat) < .n){next}

      # Second filter
      .h <- (tan(.alpha.h / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2)

      if((max(.dat$phi) - min(.dat$phi)) < pi){

        .dat.2 <- .dat[order(.dat$phi, decreasing = F), ]

      } else {

        .dat.2 <- .dat
        .dat.2$phi <- ifelse(.dat.2$phi < 1, .dat.2$phi + (2 * pi), .dat.2$phi)
        .dat.2 <- .dat.2[order(.dat.2$phi, decreasing = F), ]

      }

      .dist <- sqrt((.dat.2$x[2:nrow(.dat.2)]-.dat.2$x[1:nrow(.dat.2)-1])^2+(.dat.2$y[2:nrow(.dat.2)]-.dat.2$y[1:nrow(.dat.2)-1])^2)
      .dist <- sd(.dist) / .h

      if(.dist > 1){next}

      # Generate mesh

      .x.rang <- max(.dat$x) - min(.dat$x)
      .y.rang <- max(.dat$y) - min(.dat$y)
      # .phi.rang <- max(.dat$phi) - min(.dat$phi)
      # .rho.rang <- max(.dat$rho) - min(.dat$rho)

      # Compute centroids coordinates with respect to TLS
      .x.cent <- (.x.rang / 2) + min(.dat$x)
      .y.cent <- (.y.rang / 2) + min(.dat$y)
      # .phi.cent <- (.phi.rang / 2) + min(.dat$phi)
      # .rho.cent <- (.rho.rang / 2) + min(.dat$rho)

      # Obtain width for mesh to be applied in the cluster
      .ancho.malla <- (max(.x.rang, .y.rang) / 2) * 1.5

      # Maximal and minimal mesh coordinates
      .xmin <- .x.cent - .ancho.malla
      .ymin <- .y.cent - .ancho.malla
      .xmax <- .x.cent + .ancho.malla
      .ymax <- .y.cent + .ancho.malla

      # Obtain second mesh based on cylindrical coordinates (phi and rho), which
      # will be used to compute mean points density by cell
      # .ancho.malla.2 <- (max(.phi.rang, .rho.rang) / 2)
      #
      # .phimin <- .phi.cent - .ancho.malla.2
      # .rhomin <- .rho.cent - .ancho.malla.2
      # .phimax <- .phi.cent + .ancho.malla.2
      # .rhomax <- .rho.cent + .ancho.malla.2

      # Filter
      .h <- 2 * (tan(.alpha.h / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2)

      .x.values <- seq(from = .xmin, to = .xmax, by = .h)
      .y.values <- seq(from = .ymin, to = .ymax, by = .h)

      .h <- .h / 2

      .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

      for(.i in 1:length(.x.values)){
        for(.j in 1:length(.y.values)){

          .den <- .dat[which(.dat$x <= ((.x.values[.i]) + .h) &
                             .dat$x > ((.x.values[.i]) - .h) &
                             .dat$y <= ((.y.values[.j]) + .h) &
                             .dat$y > ((.y.values[.j]) - .h)), , drop = FALSE]

          # Discard cells with less than 2 points for computing mean points
          # density by cell
          .density[.j, .i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

        }

      }


      # Estimate mean density by cell
      .threeshold <- mean(.density, na.rm = T)

      if(is.nan(.threeshold)){next}

      .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))
      .remove <- data.frame(point = as.numeric())

      for(.i in 1:length(.x.values)){
        for(.j in 1:length(.y.values)){

          .den <- .dat[which(.dat$x <= ((.x.values[.i]) + .h) &
                             .dat$x > ((.x.values[.i]) - .h) &
                             .dat$y <= ((.y.values[.j]) + .h) &
                             .dat$y > ((.y.values[.j]) - .h)), , drop = FALSE]

          # Discard cells with less than 2 points for computing mean density by
          # cell
          .density[.j, .i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

          if(nrow(.den) > .threeshold){

            .rem <- data.frame(point = .den$point)
            .remove <- rbind(.remove, .rem)

          }

        }

      }

      .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)

      if(nrow(.dat) < 1){next}

      if(is.nan(mean(.dat$slope, na.rm = TRUE))){

        .n <- (0.1 / (tan(.alpha.v / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2))

      } else {

        .n <- (0.1 / (tan(.alpha.v / 2) * (mean(.dat$r) / cos(mean(.dat$slope, na.rm = TRUE))) * 2))

      }

      # Ratio points

      if(mean(.cut$slope, na.rm = TRUE) > 0.5){.n <- 0.7 * .n}


      # Obtain phi and rho coordinates corresponding to mesh intersections
      .x2.values <- seq(from = min(.dat$phi), to = max(.dat$phi), by = .alpha.h)

      # Define matrix where points number by cell will be stored
      .density <- vector(length = length(.x2.values))

      .remove <- data.frame(point = as.numeric())

      for(.i in 1:length(.x2.values)){

        .den <- .dat[which(.dat$phi <= ((.x2.values[.i]) + (.alpha.h/2)) &
                           .dat$phi >  ((.x2.values[.i]) - (.alpha.h/2))), ]

        # Aquellas celdas con menos de 2 puntos no las tengo en cuenta
        # para luego m?s tarde calcular la densidad media por celda
        .density[.i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))


        if(nrow(.den) > 1){

          .rem <- data.frame(point = .den$point)
          .remove <- rbind(.remove, .rem)

        }

      }


      # Estimate mean density by cell
      .density <- ifelse(is.nan(.density), NA, .density)

      if(is.nan(mean(.density, na.rm = TRUE))){next}

      if(max(.density[which(!is.na(.density))], na.rm = T) < floor(.n)){next}

      # Remove cells containing only 1 point
      .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)

      # If no points remain in .dat after removing, go to next iteration
      if(nrow(.dat) < 1){next}

      # Estimate points number for both the original cloud (.num.points) and the
      # point cloud reduced by the point cropping process (.num.points.hom)
      .num.points <- nrow(.dat)
      .num.points.hom <- nrow(.dat[which(.dat$prob.selec == 1), , drop = FALSE])

      # After this previous filtering, compute cluster centroid

      .x.values <- seq(from = .xmin, to = .xmax, by = 0.01)
      .y.values <- seq(from = .ymin, to = .ymax, by = 0.01)

      # Create an empty matrix where, for each mesh intersection, variance of
      # distances between points and corresponding intersection will be stored
      .matriz <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

      for(.i in 1:length(.x.values)){
        for(.j in 1:length(.y.values)){

          .variance <- stats::var(raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.i], .y.values[.j]), lonlat=FALSE))
          .matriz[.j, .i] <- .variance

        }
      }

      # Consider as section center the intesection where variance is minimal
      .a <- which(.matriz == min(.matriz), arr.ind = TRUE)

      .center.x <- .x.values[.a[2]]
      .center.y <- .y.values[.a[1]]

      .center.phi <- atan2(.center.y, .center.x)
      .center.phi <- ifelse(.center.phi < 0, .center.phi + (2 * pi), .center.phi)
      .center.rho <- sqrt(.center.x ^ 2 + .center.y ^ 2)
      .center.r <- sqrt(.dat$sec[1] ^ 2 + .center.rho ^ 2)
      .center.theta <- atan2(.dat$sec[1], .center.rho)

      # Radius value as the mean distance
      .radio <- mean(raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE))


      # Distances between points and center
      .dat$dist <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE)


      # Center behind tree surface
      if(stats::quantile(.dat$rho, prob = 0.05) > .center.r) {next}

      # At least 95 % of distances should be greater than .radio / 2
      .dat <- .dat[order(.dat$dist, decreasing = FALSE), , drop = FALSE]
      if(stats::quantile(.dat$dist, prob = 0.05) < (.radio / 2)) {next}


      # Compute rho coordinates for section ends

      # if((max(.dat$phi) - min(.dat$phi)) < pi){
      #
      #   .dat.2 <- .dat[order(.dat$phi, decreasing = F), , drop = FALSE]
      #
      # } else {
      #
      #   .dat.2 <- .dat
      #   .dat.2$phi <- ifelse(.dat.2$phi < 1, .dat.2$phi + (2 * pi), .dat.2$phi)
      #   .dat.2 <- .dat.2[order(.dat.2$phi, decreasing = F), , drop = FALSE]
      #
      # }

      # Evaluamos aqui el ratio

      .ratio <- nrow(.dat.2) / (.n * ((max(.dat.2$phi) - min(.dat.2$phi)) / .alpha.h))

      if(.ratio < 0.5){next}


      # Select 1st percentil, if necessary for strange points
      # It remains to be seen what happens if cluster is located in 0 +/- phi
      .pto.left <- stats::quantile(.dat.2$phi, prob = 0.01)
      .rho.left <- mean(.dat.2$rho[which(.dat.2$phi <= .pto.left)])
      .phi.left <- mean(.dat.2$phi[which(.dat.2$phi <= .pto.left)])

      # Select 99th percentil, if necessary for strange points
      .pto.right <- stats::quantile(.dat.2$phi, prob = 0.99)
      .rho.right <- mean(.dat.2$rho[which(.dat.2$phi >= .pto.right)])
      .phi.right <- mean(.dat.2$phi[which(.dat.2$phi >= .pto.right)])

      # For points in section center, select those in half the angle aperture
      # phi +/- TLS aperture .alpha
      .phi.cent <- max(.dat.2$phi) - ((max(.dat.2$phi) - min(.dat.2$phi)) / 2)
      .rho.cent <- mean(.dat.2$rho[which(round(.dat.2$phi, 3) >= round(.phi.cent - .alpha.h, 3) & round(.dat.2$phi, 3) <= round(.phi.cent + .alpha.h, 3))])

      if(is.nan(.rho.cent)){next}

      # Check rho coordinates for ends are greater than center ones
      .arc.circ <- ifelse(.rho.left > .rho.cent & .rho.right > .rho.cent, 1, 0)

      # Convert original coordinates if cluster is located in 0 +/- phi
      .phi.left <- ifelse(.phi.left > (2 * pi), .phi.left - (2 * pi), .phi.left)
      .phi.right <- ifelse(.phi.right > (2 * pi), .phi.right - (2 * pi), .phi.right)

      # If the complete circumference arc is not detected (so previous criterion
      # is not satisfied), check if there is a partial occlusion.
      # If they belong to a tree section, they are found with a systematic
      # regularity so very high correlations between their correlative
      # numeration and phi must be found when they are ordered with respect to
      # phi
      .dat.2$n <- c(1:nrow(.dat.2))
      .cor <- try(stats::cor.test(x = .dat.2$n, y = .dat.2$phi, method = 'pearson'), silent = TRUE) # cor function could be used instead
      # .cor <- try(stats::cor.test(x = .dat2$n, y = .dat2$phi, method = 'spearman'))

      # If error, go to next iteration
      if(class(.cor) == "try-error"){next} else{

        .occlusion <- .cor[[4]]

      }


      # Coefficient of variation for distances among cluster points and the estimated center
      .cv <- stats::sd(raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE)) / .radio

      if(.cv > 0.1){next}

      # Zhang et al., (2019)
      .n.w.ratio <- stats::sd(.dat$z) / sqrt(stats::sd(.dat$x) ^ 2 + stats::sd(.dat$y) ^ 2)

      if(.n.w.ratio > 1){next}

      # Results
      .salida <- data.frame(cluster = .dat$cluster[1],

                            center.x = .center.x, center.y = .center.y,
                            center.phi = .center.phi, center.rho = .center.rho,
                            center.r = .center.r, center.theta = .center.theta,

                            radius = .radio,

                            num.points = .num.points, num.points.hom = .num.points.hom,

                            phi.left = .phi.left, phi.right = .phi.right,

                            arc.circ = .arc.circ, occlusion = .occlusion)

      .filter <- rbind(.filter, .salida)

    }


    # Arch of circumference or partial arch of circumference?
    .filter$tree <- ifelse(.filter$arc.circ == 1, 1,
                           ifelse(.filter$arc.circ == 0 & .filter$occlusion > 0.995, 1, 0))
    .filter <- .filter[which(.filter$tree == 1), , drop = FALSE]

    # Dbh maximum and minimum
    .filter$tree <- ifelse(.filter$radius > (.dbh.max / 2) | .filter$radius < (.dbh.min / 2), 0, 1)
    .filter <- subset(.filter, .filter$tree == 1)


    if(nrow(.filter) < 1){

      .filter1.0 <- data.frame(cluster = as.numeric(),
                               center.x = as.numeric(), center.y = as.numeric(),
                               center.phi = as.numeric(), center.rho = as.numeric(),
                               center.r = as.numeric(), center.theta = as.numeric(),
                               radius = as.numeric(),
                               num.points = as.numeric(), num.points.hom = as.numeric(),
                               phi.left = as.numeric(), phi.right = as.numeric(),
                               arc.cir = as.numeric(), sec = as.numeric())

    } else{

      .filter1.0 <- .filter[, c("cluster",
                                "center.x", "center.y", "center.phi", "center.rho", "center.r", "center.theta",
                                "radius", "num.points", "num.points.hom", "phi.left", "phi.right", "arc.circ"), drop = FALSE]
      .filter1.0$sec <- cuts

    }

    .filteraux<-rbind(.filteraux,.filter1.0)

  }# End of cuts loop

  .filter<-.filteraux

  # Repeat this for the remaining sections!!!
  # Merge of all sections ----

  #.filter <- rbind(.filter1.0, .filter1.3, .filter1.6)

  if(nrow(.filter) < 1) {
    
    # Generate a warning and create empty data.frame to be returned, if no row
    # was included in .filter
    
    if(is.null(data$id)){
      
      # If plot identification (id) is not available
      
      warning("No tree was detected")
      
      .colnames <- c("tree", "x", "y", "phi", "phi.left", "phi.right",
                     "horizontal.distance", "dbh", "num.points",
                     "num.points.hom", "num.points.est", "num.points.hom.est",
                     "partial.occlusion")
    } else {
      
      # If plot identification (id) is available
      
      warning("No tree was detected for plot", data$id[1])
      
      .colnames <- c("id", "file", "tree", "x", "y", "phi", "phi.left",
                     "phi.right", "horizontal.distance", "dbh", "num.points",
                     "num.points.hom", "num.points.est", "num.points.hom.est",
                     "partial.occlusion")
    }
    
    .tree <- data.frame(matrix(nrow = 0, ncol = length(.colnames),
                               dimnames = list(NULL, .colnames)))
  }
  
  else {
    
    # Continue calculations, if any row was included in .filter
    
    .dbscan <- dbscan::dbscan(.filter[, c("center.x", "center.y"), drop = FALSE], eps = mean(.filter$radius), minPts = 1)
    .filter$cluster <- .dbscan$cluster
    .filter <- .filter[order(.filter$cluster, .filter$sec), , drop = FALSE]
    
    # Calculate of taper coefficient as the slope coefficient of linear regression
    
    .taper <- .filter[, c("cluster", "sec", "radius"), drop = FALSE]
    
    
    .lm <- stats::lm(radius ~ sec, data = .taper)
    .slope <- stats::coef(.lm)[2]
    
    
    # If a new column ("dif") is created containing the difference between dbh and
    # section from which radius is estimated, number of used sections (or cuts)
    # will not be important. An estimated radius will always exist
    .filter$dif <- 1.3 - .filter$sec
    .filter$radio.est <- ifelse(.filter$dif == 0, .filter$radius,
                                .filter$radius + .slope * .filter$dif)
    
    # When there are not enough tree to the linear model, and tress were
    # detected at at different sections than 1.3 m, we assume these radius
    # as dbh. This could happen in very few situations.
    .filter$radio.est <- ifelse(is.na(.filter$radio.est), .filter$radius, .filter$radio.est)
    
    .radio.est <- data.frame(radio.est = as.numeric())
    
    for (i in unique(.filter$cluster)) {
      
      .dat <- .filter[which(.filter$cluster == i), ]
      
      if(max(.dat$arc.circ) > 0){
        
        .dat <- .dat[which(.dat$arc.circ > 0), ]
        
      }
      
      .out <- data.frame(radio.est = mean(.dat$radio.est))
      .radio.est <- rbind(.radio.est, .out)
      
    }
    
    
    # Dendrometric variables
    .tree <- data.frame(tree = tapply(.filter$cluster, .filter$cluster, mean, na.rm = TRUE),
                        center.x = tapply(.filter$center.x, .filter$cluster, mean, na.rm = TRUE),
                        center.y = tapply(.filter$center.y, .filter$cluster, mean, na.rm = TRUE),
                        center.phi = tapply(.filter$center.phi, .filter$cluster, mean, na.rm = TRUE),
                        center.rho = tapply(.filter$center.rho, .filter$cluster, mean, na.rm = TRUE),
                        center.r = tapply(.filter$center.r, .filter$cluster, mean, na.rm = TRUE),
                        center.theta = tapply(.filter$center.theta, .filter$cluster, mean,na.rm = TRUE),
                        
                        horizontal.distance = tapply(.filter$center.rho, .filter$cluster, mean, na.rm = TRUE), # repeated line
                        radius = .radio.est$radio.est,
                        
                        phi.left = tapply(.filter$phi.left, .filter$cluster, mean, na.rm = TRUE),
                        phi.right = tapply(.filter$phi.right, .filter$cluster, mean, na.rm = TRUE),
                        
                        partial.occlusion = tapply(.filter$arc.circ, .filter$cluster, min, na.rm = TRUE),
                        
                        num.points = tapply(.filter$num.points, .filter$cluster, mean, na.rm = TRUE),
                        num.points.hom = tapply(.filter$num.points.hom, .filter$cluster, mean, na.rm = TRUE))
    
    # Indicate trees with partial occlusions, those for which none of the sections
    # was identified as circumference arch (ArcCirc)
    .tree$partial.occlusion <- ifelse(.tree$partial.occlusion == 0, 1, 0)
    
    # Compute dbh (cm)
    .tree$dbh <- .tree$radius * 200
    
    # Calculate points belonging to radius unit
    # Since it will be an estimation, select sections completely visible
    # (ArcCirc == 1) in section corresponding to 1.3 m (where dbh is estimated)
    .filter$filter <- ifelse(.filter$sec == 1.3 & .filter$arc.circ == 1, 1, 0)
    .filter2 <- subset(.filter, .filter$filter == 1)
    
    if(nrow(.filter2) < 1)
      .filter2 <- .filter
    
    # Estimate number of points by cluster, with and without point cropping
    # process, corresponding to radius 1 m
    .filter2$points.radio <- .filter2$num.points / .filter2$radio
    .filter2$points.radio.hom <- .filter2$num.points.hom / .filter2$radio
    
    # Average points after point cropping by m of radius
    .tree$points.m <- mean(.filter2$points.radio)
    .tree$points.m.hom <- mean(.filter2$points.radio.hom)
    
    # Finally, compute number of points estimated for each tree according to
    # radius
    .tree$num.points.est <- .tree$points.m * .tree$radius
    .tree$num.points.hom.est <- .tree$points.m.hom * .tree$radius
    
    # If plot identification (id) is not available
    if(is.null(data$id)){
      
      .tree <- .tree[, c("tree", "center.x", "center.y", "center.phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "num.points", "num.points.hom", "num.points.est", "num.points.hom.est", "partial.occlusion"), drop = FALSE]
      colnames(.tree) <- c("tree", "x", "y", "phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "num.points", "num.points.hom", "num.points.est", "num.points.hom.est", "partial.occlusion")
      
    } else{
      
      # If plot identification (id) is available
      
      .tree$id <- data$id[1]
      .tree$file <- data$file[1]
      
      .tree <- .tree[, c("id", "file", "tree", "center.x", "center.y", "center.phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "num.points", "num.points.hom", "num.points.est", "num.points.hom.est", "partial.occlusion"), drop = FALSE]
      colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "num.points", "num.points.hom", "num.points.est", "num.points.hom.est", "partial.occlusion")
      
    }
    
  }

  # .tree$id <- as.integer(.tree$id)

  # Lastly, aggregate attributes table
  if(!is.null(plot.attributes))
    .tree <- merge(.tree, plot.attributes, by = "id", all = FALSE)


  if(isTRUE(save.result)){

    utils::write.csv(.tree,
                     file = file.path(dir.result, "tree.list.tls.csv"),
                     row.names = FALSE)
  }


  #####
  return(.tree)

}
