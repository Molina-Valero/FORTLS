
.sections.single.scan.2 <- function(.dat, .cut, .alpha.v, .alpha.h, .dbh.min, .dbh.max, slice, bark.roughness, x.center, y.center){

  .filter <- data.frame(cluster = as.numeric(),

                        # Center coordinates
                        center.x = as.numeric(), center.y = as.numeric(),
                        center.phi = as.numeric(), center.rho = as.numeric(),
                        center.r = as.numeric(), center.theta = as.numeric(),

                        # Radius
                        radius = as.numeric(),

                        # Number of points belonging to cluster (raw and after point cropping)
                        n.pts = as.numeric(), n.pts.red = as.numeric(),

                        # Phi coordinates of left and right
                        phi.left = as.numeric(), phi.right = as.numeric(),

                        # Circumference arc
                        arc.circ = as.numeric(),

                        # Partial occlusion
                        sec = as.numeric(), stringsAsFactors = FALSE)


  # Estimated number of points per vertical pass

  .n <- (slice / (tan(.alpha.v / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2))

  if(nrow(.dat) < .n){return(.filter)}



  # Compute rho coordinates for section ends

  if((max(.dat$phi) - min(.dat$phi)) < pi){

    .dat.2 <- .dat[order(.dat$phi, decreasing = F), ]

  } else {

    .dat.2 <- .dat
    .dat.2$phi <- ifelse(.dat.2$phi < 1, .dat.2$phi + (2 * pi), .dat.2$phi)
    .dat.2 <- .dat.2[order(.dat.2$phi, decreasing = F), ]

  }

  .dat.3 <- .dat[order(.dat$theta, decreasing = F), ]

  .dist <- diff(.dat.2$phi)
  .dist.2 <- diff(.dat.3$theta)

  .dist <- sd(.dist) / .alpha.h
  .dist.2 <- sd(.dist.2) / .alpha.v
  .dist <- mean(c(.dist, .dist.2), na.rm = TRUE)


  if(.dist > 1){return(.filter)}

  rm(.dist, .dist.2, .dat.3)


  # Generate mesh

  .x.rang <- max(.dat$x) - min(.dat$x)
  .y.rang <- max(.dat$y) - min(.dat$y)


  # Compute centroids coordinates with respect to TLS

  .x.cent <- (.x.rang / 2) + min(.dat$x)
  .y.cent <- (.y.rang / 2) + min(.dat$y)


  # Obtain width for mesh to be applied in the cluster

  .ancho.malla <- (max(.x.rang, .y.rang) / 2) * 1.5


  # Maximal and minimal mesh coordinates

  .xmin <- .x.cent - .ancho.malla
  .ymin <- .y.cent - .ancho.malla
  .xmax <- .x.cent + .ancho.malla
  .ymax <- .y.cent + .ancho.malla


  # Expected mesh resolution according to scan parameters and distance from TLS

  .h <- 2 * (tan(.alpha.h / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2)

  # Define X and Y Grid Values

  .x.values <- seq(from = .xmin, to = .xmax, by = .h)
  .y.values <- seq(from = .ymin, to = .ymax, by = .h)

  .h <- .h / 2


  .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

  for(i in 1:length(.x.values)){
    for(j in 1:length(.y.values)){

      .den <- .dat[.dat$x <= .x.values[i] + .h &
                   .dat$x >  .x.values[i] - .h &
                   .dat$y <= .y.values[j] + .h &
                   .dat$y >  .y.values[j] - .h, , drop = FALSE]

      # Discard cells with less than 2 points for computing mean points density by cell

      .density[j, i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

    }

  }


  .threeshold <- stats::quantile(.density, prob = 0.1, na.rm = TRUE)

  if(is.nan(.threeshold) | is.na(.threeshold)){return(.filter)}



  # Initialize density matrix and remove data frame

  .density <- matrix(NA, ncol = length(.x.values), nrow = length(.y.values))
  .remove <- data.table::data.table(point = numeric())

  # Assign each point in .dat to a grid cell based on its x and y coordinates
  .dat[, grid_x := findInterval(x, .x.values)]
  .dat[, grid_y := findInterval(y, .y.values)]

  # Count the number of points in each grid cell
  point_counts <- .dat[, .N, by = .(grid_x, grid_y)]

  # Fill the density matrix based on point counts
  .density[cbind(point_counts$grid_y, point_counts$grid_x)] <- point_counts$N

  # Identify and remove points exceeding the threshold
  .points_to_remove <- .dat[.N > .threeshold, .(point), by = .(grid_x, grid_y)]
  .remove <- rbind(.remove, .points_to_remove[, .(point)])


  .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)

  .dat <- .dat[!duplicated(.dat$point), ]



  # Filter

  if(nrow(.dat) < .n){return(.filter)}


  # Obtain phi and rho coordinates corresponding to mesh intersections

  .x2.values <- seq(from = min(.dat$phi), to = max(.dat$phi), by = .alpha.h)


  # # Define matrix where points number by cell will be stored
  #
  # .density <- vector(length = length(.x2.values))
  #
  # .remove <- data.frame(point = as.numeric(), stringsAsFactors = FALSE)
  #
  # for(i in 1:length(.x2.values)){
  #
  #   .den <- .dat[.dat$phi <= .x2.values[i] + (.alpha.h/2) &
  #                .dat$phi >  .x2.values[i] - (.alpha.h/2), ]
  #
  #   # Aquellas celdas con menos de 2 puntos no las tengo en cuenta
  #   # para luego mas tarde calcular la densidad media por celda
  #   .density[i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))
  #
  #
  #   if(nrow(.den) > 1){
  #
  #     .remove <- rbind(.remove, data.frame(point = .den$point, stringsAsFactors = FALSE))
  #
  #   }
  #
  # }



  # Pre-allocate the density vector
  .density <- rep(NA, length(.x2.values))

  # Initialize a list to accumulate rows to be removed
  .remove_list <- vector("list", length(.x2.values))

  for (i in seq_along(.x2.values)) {

    # Subset data efficiently using data.table syntax
    .den <- .dat[phi <= .x2.values[i] + (.alpha.h/2) &
                   phi >  .x2.values[i] - (.alpha.h/2)]

    # Calculate density
    .density[i] <- nrow(.den)

    if (nrow(.den) > 1) {
      # Store rows in the list
      .remove_list[[i]] <- .den[, .(point)]
    }
  }

  # Combine all points to remove at once
  .remove <- do.call(rbind, lapply(.remove_list, function(x) if (!is.null(x)) data.frame(point = x, stringsAsFactors = FALSE) else NULL))

  # Remove the NULL entries from the list if there were no points to remove
  .remove <- na.omit(.remove)


  # Estimate mean density by cell
  .density <- ifelse(is.nan(.density), NA, .density)

  if(is.nan(mean(.density, na.rm = TRUE))){return(.filter)}

  if(max(.density[!is.na(.density)], na.rm = T) < floor(.n)){return(.filter)}

  # Remove cells containing only 1 point
  .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)

  # If no points remain in .dat after removing, go to next iteration
  if(nrow(.dat) < .n){return(.filter)}

  # Estimate points number for both the original cloud (.n.pts) and the
  # point cloud reduced by the point cropping process (.n.pts.red)
  .n.pts <- nrow(.dat)
  .n.pts.red <- nrow(.dat[.dat$prob.selec == 1, , drop = FALSE])

  # After this previous filtering, compute cluster centroid

  .x.values <- seq(from = .xmin, to = .xmax, by = 0.01)
  .y.values <- seq(from = .ymin, to = .ymax, by = 0.01)

  # Create an empty matrix where, for each mesh intersection, variance of
  # distances between points and corresponding intersection will be stored
  .matriz <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

  for(i in 1:length(.x.values)){
    for(j in 1:length(.y.values)){

      .variance <- stats::var(raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[i], .y.values[j]), lonlat=FALSE))
      .matriz[j, i] <- .variance

    }
  }

  # Consider as section center the intesection where variance is minimal
  .a <- which(.matriz == min(.matriz), arr.ind = TRUE)

  rm(.matriz)

  .center.x <- .x.values[.a[2]]
  .center.y <- .y.values[.a[1]]


  # Distances between points and center
  .dat$dist <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE)


  dat.i <- rep(list(as.matrix(.dat[, c("x", "y")])), 600)
  kk <- try(do.call(rbind, (lapply(dat.i, .RANSAC))), silent = TRUE)

  rm(dat.i)


  quantiles <- stats::quantile(.dat$dist, probs = c(0.05, 0.25, 0.5, 0.75, 0.99), na.rm = TRUE)


  if(class(kk)[1] == "try-error"){

    if(is.null(bark.roughness)){

      .radio <- mean(.dat$dist[.dat$dist > quantiles[1] & .dat$dist < quantiles[5]])
      .cv <- stats::sd(.dat$dist[.dat$dist > quantiles[1] & .dat$dist < quantiles[5]]) / .radio

    } else if(bark.roughness == 1){

      .radio <- mean(.dat$dist[.dat$dist > quantiles[2] & .dat$dist < quantiles[5]])
      .cv <- stats::sd(.dat$dist[.dat$dist > quantiles[2] & .dat$dist < quantiles[5]]) / .radio

    } else if(bark.roughness == 2){

      .radio <- mean(.dat$dist[.dat$dist > quantiles[3] & .dat$dist < quantiles[5]])
      .cv <- stats::sd(.dat$dist[.dat$dist > quantiles[3] & .dat$dist < quantiles[5]]) / .radio

    } else {

      .radio <- mean(.dat$dist[.dat$dist > quantiles[4] & .dat$dist < quantiles[5]])
      .cv <- stats::sd(.dat$dist[.dat$dist > quantiles[4] & .dat$dist < quantiles[5]]) / .radio

    }

    if(is.na(.cv) | .radio <= 0 | is.na(.radio)){return(.filter)}


  } else {

    kk <- kk[kk$n >= max(kk$n), ]

    if(nrow(kk) > 1)
      kk <- kk[kk$mae <= min(kk$mae), ]

    kk <- kk[1, ]

    # .dat$distRANSAC <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.centerRANSAC$X, .centerRANSAC$Y), lonlat = FALSE)
    .dat$distRANSAC <- raster::pointDistance(cbind(.dat$x,.dat$y), c(kk$x, kk$y), lonlat = FALSE)

    # Radius value as the mean distance
    if(is.null(bark.roughness)){
      .radio <- mean(.dat$dist[.dat$dist > quantiles[1] & .dat$dist < quantiles[5]])
      .cv <- stats::sd(.dat$dist[.dat$dist > quantiles[1] & .dat$dist < quantiles[5]]) / .radio
      .radioRANSAC <- kk$radio
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else if(bark.roughness == 1){
      .radio <- mean(.dat$dist[.dat$dist > quantiles[2] & .dat$dist < quantiles[5]])
      .cv <- stats::sd(.dat$dist[.dat$dist > quantiles[2] & .dat$dist < quantiles[5]]) / .radio
      .radioRANSAC <- kk$radio
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else if(bark.roughness == 2){
      .radio <- mean(.dat$dist[.dat$dist > quantiles[3] & .dat$dist < quantiles[5]])
      .cv <- stats::sd(.dat$dist[.dat$dist > quantiles[3] & .dat$dist < quantiles[5]]) / .radio
      .radioRANSAC <- kk$radio
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else {
      .radio <- mean(.dat$dist[.dat$dist > quantiles[4] & .dat$dist < quantiles[5]])
      .cv <- stats::sd(.dat$dist[.dat$dist > quantiles[4] & .dat$dist < quantiles[5]]) / .radio
      .radioRANSAC <- kk$radio
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}
    }

    if(is.na(.cv)){return(.filter)}

    if(1 * .cv >= .cvRANSAC){

      .radio <- .radioRANSAC
      .cv <- .cvRANSAC
      .dat$dist <- .dat$distRANSAC

      .center.x <- kk$x
      .center.y <- kk$y}

    .dat <- .dat[, 1:(ncol(.dat)-1)]

    if(.radio <= 0 | is.na(.radio)){return(.filter)}

  }

  .center.phi <- atan2(.center.y-y.center, .center.x-x.center)
  .center.phi <- ifelse(.center.phi < 0, .center.phi + (2 * pi), .center.phi)
  .center.rho <- sqrt((.center.x-x.center) ^ 2 + (.center.y-y.center) ^ 2)
  .center.r <- sqrt(.dat$sec[1] ^ 2 + .center.rho ^ 2)
  .center.theta <- atan2(.dat$sec[1], .center.rho)

  if(is.na(.cv) | .cv > 0.1 | length(.dat$dist[.dat$dist > quantiles[2]]) < 2) {return(.filter)}

  # Center behind tree surface
  if(stats::quantile(.dat$rho, prob = 0.05, na.rm = T) > .center.r) {return(.filter)}

  # At least 95 % of distances should be greater than .radio / 2
  if(stats::quantile(.dat$dist, prob = 0.05, na.rm = T) < (.radio / 2)) {return(.filter)}


  # Select 1st percentil, if necessary for strange points
  # It remains to be seen what happens if cluster is located in 0 +/- phi
  .pto.left <- stats::quantile(.dat.2$phi, prob = 0.01, na.rm = T)
  .rho.left <- mean(.dat.2$rho[which(.dat.2$phi <= .pto.left)])
  .phi.left <- mean(.dat.2$phi[which(.dat.2$phi <= .pto.left)])

  # Select 99th percentil, if necessary for strange points
  .pto.right <- stats::quantile(.dat.2$phi, prob = 0.99, na.rm = T)
  .rho.right <- mean(.dat.2$rho[which(.dat.2$phi >= .pto.right)])
  .phi.right <- mean(.dat.2$phi[which(.dat.2$phi >= .pto.right)])

  # For points in section center, select those in half the angle aperture
  # phi +/- TLS aperture .alpha
  .phi.cent <- max(.dat.2$phi) - ((max(.dat.2$phi) - min(.dat.2$phi)) / 2)
  .rho.cent <- mean(.dat.2$rho[which(round(.dat.2$phi, 3) >= round(.phi.cent - .alpha.h, 3) & round(.dat.2$phi, 3) <= round(.phi.cent + .alpha.h, 3))])

  if(is.nan(.rho.cent)){return(.filter)}

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
  .cor <- try(stats::cor.test(x = .dat.2$n, y = .dat.2$phi, method = 'pearson'), silent = TRUE)

  # If error, go to next iteration
  if(methods::is(.cor) == "try-error"){return(.filter)} else{

    .occlusion <- .cor[[4]]

  }


  # Zhang et al., (2019)
  .n.w.ratio <- stats::sd(.dat$z) / sqrt(stats::sd(.dat$x) ^ 2 + stats::sd(.dat$y) ^ 2)

  if(.n.w.ratio > 1 | is.nan(.n.w.ratio)){return(.filter)}

  if(nrow(.dat) < .n){return(.filter)}


  # Results
  .filter <- data.frame(cluster = .dat$cluster[1],

                        center.x = .center.x, center.y = .center.y,
                        center.phi = .center.phi, center.rho = .center.rho,
                        center.r = .center.r, center.theta = .center.theta,

                        radius = .radio,

                        n.pts = .n.pts, n.pts.red = .n.pts.red,

                        phi.left = .phi.left, phi.right = .phi.right,

                        arc.circ = .arc.circ, occlusion = .occlusion,

                        stringsAsFactors = FALSE)


  # Arch of circumference or partial arch of circumference?
  .filter$tree <- ifelse(.filter$arc.circ == 1, 1,
                         ifelse(.filter$arc.circ == 0 & .filter$occlusion > 0.975, 1, 0))
  .filter <- .filter[which(.filter$tree == 1), , drop = FALSE]


  if(nrow(.filter) < 1){

    .filter1.0 <- data.frame(cluster = as.numeric(),
                             center.x = as.numeric(), center.y = as.numeric(),
                             center.phi = as.numeric(), center.rho = as.numeric(),
                             center.r = as.numeric(), center.theta = as.numeric(),
                             radius = as.numeric(),
                             n.pts = as.numeric(), n.pts.red = as.numeric(),
                             phi.left = as.numeric(), phi.right = as.numeric(),
                             arc.cir = as.numeric(), sec = as.numeric(),
                             stringsAsFactors = FALSE)

  } else{

    .filter1.0 <- .filter[, c("cluster",
                              "center.x", "center.y", "center.phi", "center.rho", "center.r", "center.theta",
                              "radius", "n.pts", "n.pts.red", "phi.left", "phi.right", "arc.circ"), drop = FALSE]

    .filter1.0$sec <- .dat$sec[1]

  }


  return(.filter1.0)

}

