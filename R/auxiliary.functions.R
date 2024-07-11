

.sections.single.scan <- function(cut, .cut, .alpha.v, .alpha.h, .dbh.min, .dbh.max, slice, bark.roughness, x.center, y.center){

  .filter <- data.frame(cluster = as.numeric(),

                        # Center coordinates
                        center.x = as.numeric(), center.y = as.numeric(),
                        center.phi = as.numeric(), center.rho = as.numeric(),
                        center.r = as.numeric(), center.theta = as.numeric(),

                        # Radius
                        radius = as.numeric(),

                        # Number of points belonging to cluster (craw and after point cropping)
                        n.pts = as.numeric(), n.pts.red = as.numeric(),

                        # Phi coordinates of left and right
                        phi.left = as.numeric(), phi.right = as.numeric(),

                        # Circumference arc
                        arc.circ = as.numeric(),

                        # Partial occlusion
                        sec = as.numeric())


  # Select cluster i
  .dat <- cut

  if(nrow(.dat) < 10){return(.filter)}
  # print(1)

  # First filter

  .n <- (slice / (tan(.alpha.v / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2))

  if(nrow(.dat) < .n){return(.filter)}
  # print(2)

  # Second filter

  .h <- 3 * (tan(.alpha.h / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2)

  # Compute rho coordinates for section ends

  if((max(.dat$phi) - min(.dat$phi)) < pi){

    .dat.2 <- .dat[order(.dat$phi, decreasing = F), ]

  } else {

    .dat.2 <- .dat
    .dat.2$phi <- ifelse(.dat.2$phi < 1, .dat.2$phi + (2 * pi), .dat.2$phi)
    .dat.2 <- .dat.2[order(.dat.2$phi, decreasing = F), ]

  }

  .dat.3 <- .dat[order(.dat$theta, decreasing = F), ]


  # .dist <- sqrt((.dat.2$x[2:nrow(.dat.2)]-.dat.2$x[1:nrow(.dat.2)-1])^2+(.dat.2$y[2:nrow(.dat.2)]-.dat.2$y[1:nrow(.dat.2)-1])^2)
  # .dist <- sd(.dist) / .h
  .dist <- diff(.dat.2$phi)
  .dist.2 <- diff(.dat.3$theta)

  .dist <- sd(.dist[.dist > 0]) / .alpha.h
  .dist.2 <- sd(.dist.2[.dist.2 > 0]) / .alpha.v
  .dist <- mean(c(.dist, .dist.2))

  if(.dist > 1){return(.filter)}
  # print(3)

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


  # Filter

  .h <- 2 * (tan(.alpha.h / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2)

  .x.values <- seq(from = .xmin, to = .xmax, by = .h)
  .y.values <- seq(from = .ymin, to = .ymax, by = .h)

  .h <- .h / 2

  .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

  for(i in 1:length(.x.values)){
    for(j in 1:length(.y.values)){

      .den <- .dat[which(.dat$x <= ((.x.values[i]) + .h) &
                         .dat$x > ((.x.values[i]) - .h) &
                         .dat$y <= ((.y.values[j]) + .h) &
                         .dat$y > ((.y.values[j]) - .h)), , drop = FALSE]

      # Discard cells with less than 2 points for computing mean points density by cell

      .density[j, i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

    }

  }


  # Estimate mean density by cell
  # .threeshold <- mean(.density, na.rm = T)
  .threeshold <- stats::quantile(.density, prob = 0.1, na.rm = T)

  if(is.nan(.threeshold) | is.na(.threeshold)){return(.filter)}

  .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

  .remove <- data.frame(point = as.numeric())

  for(i in 1:length(.x.values)){
    for(j in 1:length(.y.values)){

      .den <- .dat[which(.dat$x <= ((.x.values[i]) + .h) &
                           .dat$x > ((.x.values[i]) - .h) &
                           .dat$y <= ((.y.values[j]) + .h) &
                           .dat$y > ((.y.values[j]) - .h)), , drop = FALSE]

      # Discard cells with less than 2 points for computing mean density by cell

      .density[j, i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

      if(nrow(.den) > .threeshold){

        .rem <- data.frame(point = .den$point)
        .remove <- rbind(.remove, .rem)

      }

    }

  }

  .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)

  .dat <- .dat[!duplicated(.dat$point), ]


  if(nrow(.dat) < 5){return(.filter)}
  # print(4)

  if(is.nan(mean(.dat$slope, na.rm = TRUE))){

    .n <- (slice / (tan(.alpha.v / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2))

  } else {

    .n <- (slice / (tan(.alpha.v / 2) * (mean(.dat$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2))

  }

  # Ratio points

  if(mean(.cut$slope, na.rm = TRUE) > 0.5){.n <- 0.7 * .n}
  # print(5)


  # Obtain phi and rho coordinates corresponding to mesh intersections
  .x2.values <- seq(from = min(.dat$phi), to = max(.dat$phi), by = .alpha.h)

  # Define matrix where points number by cell will be stored
  .density <- vector(length = length(.x2.values))

  .remove <- data.frame(point = as.numeric())

  for(i in 1:length(.x2.values)){

    .den <- .dat[which(.dat$phi <= ((.x2.values[i]) + (.alpha.h/2)) &
                         .dat$phi >  ((.x2.values[i]) - (.alpha.h/2))), ]

    # Aquellas celdas con menos de 2 puntos no las tengo en cuenta
    # para luego m?s tarde calcular la densidad media por celda
    .density[i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))


    if(nrow(.den) > 1){

      .rem <- data.frame(point = .den$point)
      .remove <- rbind(.remove, .rem)

    }

  }


  # Estimate mean density by cell
  .density <- ifelse(is.nan(.density), NA, .density)

  if(is.nan(mean(.density, na.rm = TRUE))){return(.filter)}
  # print(6)

  if(max(.density[which(!is.na(.density))], na.rm = T) < floor(.n)){return(.filter)}
  # print(7)

  # Remove cells containing only 1 point
  .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)

  # If no points remain in .dat after removing, go to next iteration
  if(nrow(.dat) < 5){return(.filter)}

  # Estimate points number for both the original cloud (.n.pts) and the
  # point cloud reduced by the point cropping process (.n.pts.red)
  .n.pts <- nrow(.dat)
  .n.pts.red <- nrow(.dat[which(.dat$prob.selec == 1), , drop = FALSE])

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

  .center.x <- .x.values[.a[2]]
  .center.y <- .y.values[.a[1]]


  # Distances between points and center
  .dat$dist <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE)


  # dat.i <- rep(list(as.matrix(.dat[, c("x", "y")])), 100)
  dat.i <- rep(list(as.matrix(.dat[, c("x", "y")])), 600)
  kk <- try(do.call(rbind, (lapply(dat.i, .RANSAC))), silent = TRUE)
  # print(.dat$cluster[1])

  rm(dat.i)

  # colnames(.datRANSAC) <- c("X", "Y")

  # .centerRANSAC <- suppressWarnings(try(rTLS::circleRANSAC(data.table::setDT(.datRANSAC),
  #                                                          fpoints = 0.2, pconf = 0.95, poutlier = c(0.75, 0.75), max_iterations = 100, plot = FALSE), silent = TRUE))

  if(class(kk)[1] == "try-error"){
    # if(max(kk$n) < 3){

    if(is.null(bark.roughness)){

      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.05, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.05, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio

    } else if(bark.roughness == 1){

      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio

    } else if(bark.roughness == 2){

      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio

    } else {

      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.75, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.75, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio

    }

    if(.radio <= 0 | is.na(.radio)){return(.filter)}

  } else {

    kk <- kk[kk$n >= max(kk$n), ]
    if(nrow(kk) > 1)
      kk <- kk[kk$mae <= min(kk$mae), ]

    kk <- kk[1, ]

    # .dat$distRANSAC <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.centerRANSAC$X, .centerRANSAC$Y), lonlat = FALSE)
    .dat$distRANSAC <- raster::pointDistance(cbind(.dat$x,.dat$y), c(kk$x, kk$y), lonlat = FALSE)

    # Radius value as the mean distance
    if(is.null(bark.roughness)){
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.05, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.05, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio
      # .radioRANSAC <- .centerRANSAC$radius
      .radioRANSAC <- kk$radio
      # .cvRANSAC <- stats::sd(.dat$distRANSAC[.dat$distRANSAC>stats::quantile(.dat$distRANSAC, prob = 0.05, na.rm = T) & .dat$distRANSAC<stats::quantile(.dat$distRANSAC, prob = 0.95, na.rm = T)]) / .radioRANSAC
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else if(bark.roughness == 1){
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio
      # .radioRANSAC <- .centerRANSAC$radius
      .radioRANSAC <- kk$radio
      # .cvRANSAC <- stats::sd(.dat$distRANSAC[.dat$distRANSAC>stats::quantile(.dat$distRANSAC, prob = 0.25, na.rm = T) & .dat$distRANSAC<stats::quantile(.dat$distRANSAC, prob = 0.95, na.rm = T)]) / .radioRANSAC
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else if(bark.roughness == 2){
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio
      # .radioRANSAC <- .centerRANSAC$radius
      .radioRANSAC <- kk$radio
      # .cvRANSAC <- stats::sd(.dat$distRANSAC[.dat$distRANSAC>stats::quantile(.dat$distRANSAC, prob = 0.5, na.rm = T) & .dat$distRANSAC<stats::quantile(.dat$distRANSAC, prob = 0.95, na.rm = T)]) / .radioRANSAC
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else {
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.75, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.75, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio
      # .radioRANSAC <- .centerRANSAC$radius
      .radioRANSAC <- kk$radio
      # .cvRANSAC <- stats::sd(.dat$distRANSAC[.dat$distRANSAC>stats::quantile(.dat$distRANSAC, prob = 0.75, na.rm = T) & .dat$distRANSAC<stats::quantile(.dat$distRANSAC, prob = 0.95, na.rm = T)]) / .radioRANSAC
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

  if(.cv > 0.1 | length(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T)]) < 2) {return(.filter)}
  # print(8)

  # Center behind tree surface
  if(stats::quantile(.dat$rho, prob = 0.05, na.rm = T) > .center.r) {return(.filter)}
  # print(9)

  # At least 95 % of distances should be greater than .radio / 2
  if(stats::quantile(.dat$dist, prob = 0.05, na.rm = T) < (.radio / 2)) {return(.filter)}
  # print(10)


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
  # print(11)

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
  # print(12)

  if(nrow(.dat) < 5){return(.filter)}
  # print(13)


  # Results
  .filter <- data.frame(cluster = .dat$cluster[1],

                        center.x = .center.x, center.y = .center.y,
                        center.phi = .center.phi, center.rho = .center.rho,
                        center.r = .center.r, center.theta = .center.theta,

                        radius = .radio,

                        n.pts = .n.pts, n.pts.red = .n.pts.red,

                        phi.left = .phi.left, phi.right = .phi.right,

                        arc.circ = .arc.circ, occlusion = .occlusion)


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
                             arc.cir = as.numeric(), sec = as.numeric())

  } else{

    .filter1.0 <- .filter[, c("cluster",
                              "center.x", "center.y", "center.phi", "center.rho", "center.r", "center.theta",
                              "radius", "n.pts", "n.pts.red", "phi.left", "phi.right", "arc.circ"), drop = FALSE]
    .filter1.0$sec <- cut$sec[1]

  }

  return(.filter1.0)

}





.sections.multi.scan <- function(cut, tls.precision, .dbh.min, .dbh.max, slice, bark.roughness, x.center, y.center){

  .filter <- data.frame(cluster = as.numeric(),

                        # Center coordinates
                        center.x = as.numeric(), center.y = as.numeric(),
                        center.phi = as.numeric(), center.rho = as.numeric(),
                        center.r = as.numeric(), center.theta = as.numeric(),

                        # Radius
                        radius = as.numeric(),

                        # Number of points belonging to cluster (craw and after point cropping)
                        n.pts = as.numeric(), n.pts.red = as.numeric(),

                        # Circumference
                        circ = as.numeric(),

                        # Circumference arc
                        arc.circ = as.numeric(),

                        # Partial occlusion
                        occlusion = as.numeric())


  .dat <- cut

  if(nrow(.dat) < 10){return(.filter)}

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

  # Filter
  if(is.null(tls.precision)){.h <- 0.03} else {.h <- tls.precision * 2 }

  .x.values <- seq(from = .xmin, to = .xmax, by = .h)
  .y.values <- seq(from = .ymin, to = .ymax, by = .h)

  .h <- .h / 2

  .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

  for(i in 1:length(.x.values)){
    for(j in 1:length(.y.values)){

      .den <- .dat[which(.dat$x <= ((.x.values[i]) + .h) &
                           .dat$x > ((.x.values[i]) - .h) &
                           .dat$y <= ((.y.values[j]) + .h) &
                           .dat$y > ((.y.values[j]) - .h)), , drop = FALSE]

      # Discard cells with less than 2 points for computing mean points
      # density by cell
      .density[j, i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

    }

  }


  # Estimate mean density by cell
  # .threeshold <- stats::median(.density, na.rm = T)
  # .threeshold <- mean(.density, na.rm = T)
  .threeshold <- stats::quantile(.density, prob = 0.1, na.rm = T)

  if(is.nan(.threeshold) | is.na(.threeshold)){return(.filter)}

  .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))
  .remove <- data.frame(point = as.numeric())

  for(i in 1:length(.x.values)){
    for(j in 1:length(.y.values)){

      .den <- .dat[which(.dat$x <= ((.x.values[i]) + .h) &
                         .dat$x > ((.x.values[i]) - .h) &
                         .dat$y <= ((.y.values[j]) + .h) &
                         .dat$y > ((.y.values[j]) - .h)), , drop = FALSE]

      # Discard cells with less than 2 points for computing mean density by
      # cell
      .density[j, i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

      if(nrow(.den) > .threeshold){

        .rem <- data.frame(point = .den$point)
        .remove <- rbind(.remove, .rem)

      }

    }

  }

  .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)
  .dat <- .dat[!duplicated(.dat$point), ]

  if(nrow(.dat) < 10){return(.filter)}


  # Estimate points number for both the original cloud (.n.pts) and the
  # point cloud reduced by the point cropping process (.n.pts.red)
  .n.pts <- nrow(.dat)
  .n.pts.red <- nrow(.dat[which(.dat$prob.selec == 1), , drop = FALSE])

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

  .center.x <- .x.values[.a[2]]
  .center.y <- .y.values[.a[1]]


  # Distances between points and center
  .dat$dist <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE)

  dat.i <- rep(list(as.matrix(.dat[, c("x", "y")])), 600)

  start_time <- Sys.time()
  kk <- try(do.call(rbind, (lapply(dat.i, .RANSAC))), silent = TRUE)
  end_time <- Sys.time()
  end_time - start_time

  # print(.dat$cluster[1])

  rm(dat.i)

  # colnames(.datRANSAC) <- c("X", "Y")

  # .centerRANSAC <- suppressWarnings(try(rTLS::circleRANSAC(data.table::setDT(.datRANSAC),
  #                                                          fpoints = 0.2, pconf = 0.95, poutlier = c(0.75, 0.75), max_iterations = 100, plot = FALSE), silent = TRUE))

  if(class(kk)[1] == "try-error"){
  # if(max(kk$n) < 3){

    if(is.null(bark.roughness)){

      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.05, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.05, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio

    } else if(bark.roughness == 1){

      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio

    } else if(bark.roughness == 2){

      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio

    } else {

      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.75, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.75, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio

    }

    if(.radio <= 0 | is.na(.radio)){return(.filter)}

  } else {

    kk <- kk[kk$n >= max(kk$n), ]
    if(nrow(kk) > 1)
      kk <- kk[kk$mae <= min(kk$mae), ]

    kk <- kk[1, ]

    # .dat$distRANSAC <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.centerRANSAC$X, .centerRANSAC$Y), lonlat = FALSE)
    .dat$distRANSAC <- raster::pointDistance(cbind(.dat$x,.dat$y), c(kk$x, kk$y), lonlat = FALSE)

    # Radius value as the mean distance
    if(is.null(bark.roughness)){
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.05, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.05, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio
      # .radioRANSAC <- .centerRANSAC$radius
      .radioRANSAC <- kk$radio
      # .cvRANSAC <- stats::sd(.dat$distRANSAC[.dat$distRANSAC>stats::quantile(.dat$distRANSAC, prob = 0.05, na.rm = T) & .dat$distRANSAC<stats::quantile(.dat$distRANSAC, prob = 0.95, na.rm = T)]) / .radioRANSAC
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else if(bark.roughness == 1){
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio
      # .radioRANSAC <- .centerRANSAC$radius
      .radioRANSAC <- kk$radio
      # .cvRANSAC <- stats::sd(.dat$distRANSAC[.dat$distRANSAC>stats::quantile(.dat$distRANSAC, prob = 0.25, na.rm = T) & .dat$distRANSAC<stats::quantile(.dat$distRANSAC, prob = 0.95, na.rm = T)]) / .radioRANSAC
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else if(bark.roughness == 2){
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio
      # .radioRANSAC <- .centerRANSAC$radius
      .radioRANSAC <- kk$radio
      # .cvRANSAC <- stats::sd(.dat$distRANSAC[.dat$distRANSAC>stats::quantile(.dat$distRANSAC, prob = 0.5, na.rm = T) & .dat$distRANSAC<stats::quantile(.dat$distRANSAC, prob = 0.95, na.rm = T)]) / .radioRANSAC
      .cvRANSAC <- kk$cv
      if(is.na(.cvRANSAC)){.cvRANSAC <- 9999}

    } else {
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.75, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)])
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.75, na.rm = T) & .dat$dist<stats::quantile(.dat$dist, prob = 0.99, na.rm = T)]) / .radio
      # .radioRANSAC <- .centerRANSAC$radius
      .radioRANSAC <- kk$radio
      # .cvRANSAC <- stats::sd(.dat$distRANSAC[.dat$distRANSAC>stats::quantile(.dat$distRANSAC, prob = 0.75, na.rm = T) & .dat$distRANSAC<stats::quantile(.dat$distRANSAC, prob = 0.95, na.rm = T)]) / .radioRANSAC
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


  if(.cv > 0.1 | length(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25, na.rm = T)]) < 2){return(.filter)}

  # At least 95 % of distances should be greater than .radio / 2
  if(stats::quantile(.dat$dist, prob = 0.05, na.rm = T) < (.radio / 2)){return(.filter)}


  .dat.2 <- .dat[order(.dat$x, .dat$y, decreasing = F), ]

  .dat.2$a <- sqrt((.dat.2$x - .dat.2$x[1])^2+(.dat.2$y - .dat.2$y[1])^2)
  .dat.2$b <- sqrt((.dat.2$x[1] - .center.x)^2+(.dat.2$y[1] - .center.y)^2)
  .dat.2$c <- sqrt((.dat.2$x - .center.x)^2+(.dat.2$y - .center.y)^2)

  .dat.2$alpha <- suppressWarnings(acos(-(.dat.2$a^2-.dat.2$b^2-.dat.2$c^2)/(2*.dat.2$b*.dat.2$c)))

  # Checking if tree section belong to a full circumference

  k1 <- ifelse(nrow(.dat[.dat$x > .center.x & .dat$y > .center.y, ]) > 0, 1, 0)
  k2 <- ifelse(nrow(.dat[.dat$x > .center.x & .dat$y < .center.y, ]) > 0, 1, 0)
  k3 <- ifelse(nrow(.dat[.dat$x < .center.x & .dat$y > .center.y, ]) > 0, 1, 0)
  k4 <- ifelse(nrow(.dat[.dat$x < .center.x & .dat$y < .center.y, ]) > 0, 1, 0)

  .circ <- ifelse(sum(k1, k2, k3, k4) == 4, 1, 0)
  .arc.circ <- ifelse(sum(k1, k2, k3, k4) == 3 | sum(k1, k2, k3, k4) == 2, 1, 0)

  # If the complete circumference arc is not detected (so previous criterion
  # is not satisfied), check if there is a partial occlusion.
  # If they belong to a tree section, they are found with a systematic
  # regularity so very high correlations between their correlative
  # numeration and phi must be found when they are ordered with respect to phi

  .dat.2 <- .dat.2[order(.dat.2$alpha, decreasing = F), ]

  .dat.2$n <- c(1:nrow(.dat.2))
  .cor <- try(stats::cor.test(x = .dat.2$n, y = .dat.2$alpha, method = 'pearson'), silent = TRUE) # cor function could be used instead


  # If error, go to next iteration
  if(methods::is(.cor) == "try-error"){return(.filter)} else{

    .occlusion <- .cor[[4]]

  }


  # Zhang et al., (2019)
  .n.w.ratio <- stats::sd(.dat$z) / sqrt(stats::sd(.dat$x) ^ 2 + stats::sd(.dat$y) ^ 2)

  if(.n.w.ratio > 1 | .n.w.ratio < 0.1 | is.nan(.n.w.ratio)){return(.filter)}

  .densidad_radio <- .n.pts.red / .radio


  if(nrow(.dat) < 10){return(.filter)}


  # Results
  .filter <- data.frame(cluster = .dat$cluster[1],

                        center.x = .center.x, center.y = .center.y,
                        center.phi = .center.phi, center.rho = .center.rho,
                        center.r = .center.r, center.theta = .center.theta,

                        radius = .radio,

                        n.pts = .n.pts, n.pts.red = .n.pts.red,

                        circ = .circ, arc.circ = .arc.circ, occlusion = .occlusion,

                        density.radio = .densidad_radio)

  # if(!exists(".salida")){next} else {.filter <- rbind(.filter, .salida)}


# Arch of circumference or partial arch of circumference?
.Q1 <- stats::quantile(.filter$density.radio, prob = 0.25, na.rm = T)
.Q3 <- stats::quantile(.filter$density.radio, prob = 0.75, na.rm = T)
.outliers <- .Q1 - 1.5 * (.Q3 - .Q1)

# Minimum number of points
# .filter$tree <- ifelse(.filter$n.pts.red > .outliers, 1, 0)
# .filter <- .filter[which(.filter$tree == 1), , drop = FALSE]

.filter$tree <- ifelse(.filter$circ == 1 & .filter$density.radio >= .outliers, 1,
                       ifelse(.filter$arc.circ == 1 & .filter$occlusion >= 0.95 & .filter$density.radio >= .outliers, 1,
                              ifelse(.filter$circ == 0 & .filter$arc.circ == 0 & .filter$occlusion >= 0.975 & .filter$density.radio >= .outliers, 1, 0)))
.filter <- .filter[which(.filter$tree == 1), , drop = FALSE]


if(nrow(.filter) < 1){

  .filter1.0 <- data.frame(cluster = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(),
                           n.pts = as.numeric(), n.pts.red = as.numeric(),
                           circ = as.numeric(), arc.circ = as.numeric(), sec = as.numeric())

  } else{

  .filter1.0 <- .filter[, c("cluster",
                            "center.x", "center.y", "center.phi", "center.rho", "center.r", "center.theta",
                            "radius", "n.pts", "n.pts.red", "circ", "arc.circ"), drop = FALSE]

  .filter1.0$sec <- cut$sec[1]

  }


  return(.filter1.0)

}


.n.w.ratio <- function(stem){

  n.w.ratio <- stats::sd(stem$z) / sqrt(stats::sd(stem$x) ^ 2 + stats::sd(stem$y) ^ 2)
  out <- data.frame(tree = stem$tree[1], n.w.ratio = n.w.ratio)
  return(out)

}



.volume <- function(data, d.top = NULL, id, den.type){

  datos <- data.frame(tree = as.numeric(),
                      x = as.numeric(), y = as.numeric(),
                      hi = as.numeric(), dhi = as.numeric(),
                      h = as.numeric(), dbh = as.numeric())

  for (i in unique(data$tree)) {

    tree <- data[data$tree == i, ]
    tree <- tree[tree$hi < tree$h, ]

    if (nrow(tree) < 2 | suppressWarnings(min(abs(diff(tree$dhi)))) > 10) {

      datos <- rbind(datos, tree[tree$hi == 1.3, c("tree", "x", "y", "hi", "dhi", "h", "dbh")])

    } else {

      tree$dif.sec <- c(abs(diff(tree$hi)), 0)
      tree$dif <- c(diff(tree$dhi),
                    tree$dhi[nrow(tree)] - tree$dhi[nrow(tree) - 1])
      tree$dif <- ifelse(tree$dif.sec > 1, tree$dif / tree$dif.sec, tree$dif)

      while (suppressWarnings(max(abs(tree$dif))) > 10) {

        tree$filter <- ifelse(tree$dif > 10 | tree$dif < -10, 0, 1)
        tree <- tree[tree$filter == 1, ]
        if(nrow(tree) < 2) {next}

        tree$dif.sec <- c(abs(diff(tree$hi)), 0)
        tree$dif <- c(diff(tree$dhi),
                      tree$dhi[nrow(tree)] - tree$dhi[nrow(tree) - 1])
        tree$dif <- ifelse(tree$dif.sec > 1, tree$dif / tree$dif.sec, tree$dif)

      }

    }

    datos <- rbind(datos, tree[, c("tree", "x", "y", "hi", "dhi", "h", "dbh")])

    # Add additional values if the largest break is smaller than total height
    # minus 0.5 m
    if (all(tree$h - .5 > tree$hi)) {

      aux <- tree[which.max(tree$hi)[1], , drop = FALSE]
      datos <- rbind(datos,
                     data.frame(tree = aux$tree, x = aux$x, y = aux$y,
                                hi = aux$h - .5,
                                dhi = .5 * aux$dhi / (aux$h - aux$hi),
                                aux[, c("h", "dbh"), drop = FALSE]))

    }

  }

  # Removing NA's

  datos <- datos[!is.na(datos$hi), ]

  loes <- stats::lowess(datos$hi, datos$dhi)
  loes <- data.frame(hi = loes$x, dhi.mean = loes$y)
  loes <- loes[!duplicated(loes), ]
  datos <- merge(datos, loes, all = FALSE)
  datos$sep <- abs(datos$dhi - datos$dhi.mean)
  datos <- datos[datos$sep < stats::quantile(datos$sep, prob = 0.9, na.rm = T), ]

  datos$id <- id


  if(nrow(datos) < 3){

    if(is.null(d.top)){

      n <- den.type

      v <- pi * (unique(data$h) ^ (n + 1) / (n + 1)) * ((unique(data$dbh) / 200) ^ 2 / (unique(data$h) - 1.3) ^ n)

      volume <- data.frame(tree = unique(data$tree), v = v)}

    else {

      n <- den.type

      v <- pi * (unique(data$h) ^ (n + 1) / (n + 1)) * ((unique(data$dbh) / 200) ^ 2 / (unique(data$h) - 1.3) ^ n)
      h.lim <- (((d.top / 200) ^ 2) / ((unique(data$dbh) / 200) ^ 2 / (unique(data$h) - 1.3) ^ n)) ^ (1 / n)
      v.com <- pi * ((unique(data$h) ^ (n + 1) - h.lim ^ (n + 1)) / (n + 1)) * ((unique(data$dbh) / 200) ^ 2 / (unique(data$h) - 1.3) ^ n)
      v.com <- ifelse(v.com < 0, 0, v.com)

      volume <- data.frame(tree = unique(data$tree), v = v, v.com = v.com, h.com = h.lim)}

  } else {


  ajuste <- stats::nls(dhi ~ dbh * ((h - hi) / (h - 1.3)) ** b1, data = datos,
                       start = c(b1 = 1), max)
  b1 <- stats::coef(ajuste)[1]

  # Height where diameter limit is reached (d.lim)
  h_d_lim <- function(dbh, h, d.lim, b1) {

    return(h - ((d.lim * (h - 1.3) ** b1) / dbh) ** (1 / b1))

  }

  # Volume between two heights
  vol_m3 <- function(dbh, h, hinf, hsup, b1) {

    return(pi * (dbh ** 2) *
             (((h - hinf) ** (2 * b1 + 1)) - ((h - hsup) ** (2 * b1 + 1))) /
             (40000 * ((h - 1.3) ** (2 * b1)) * (2 * b1 + 1)))

  }

  if (is.null(d.top)) {

    volume <- data.frame(tree = as.numeric(), v = as.numeric())

    for (i in unique(data$tree)) {

      tree <- data[data$tree == i, ]
      volume.i <- data.frame(tree = i,
                             v = vol_m3(tree$dbh[1], tree$h[1], 0, tree$h[1],
                                        b1))
      volume <- rbind(volume, volume.i)

    }

  } else {

    volume <- data.frame(tree = as.numeric(), v = as.numeric(),
                         v.com = as.numeric())

    for (i in unique(data$tree)) {

      tree <- data[data$tree == i, ]

      # If diameter at the bottom is smaller than 'd.top' then 'h.lim' is forced
      # to be zero (consequently, 'v.com' will be also zero)
      if (stats::predict(ajuste, data.frame(hi = 0, h = tree$h[1],
                                     dbh = tree$dbh[1])) < d.top) {

        h.lim <- 0

      } else {

        h.lim <- h_d_lim(tree$dbh[1], tree$h[1], d.top, b1)

      }

      volume.i <- data.frame(tree = i,
                             v = vol_m3(tree$dbh[1], tree$h[1], 0, tree$h[1], b1),
                             v.com = vol_m3(tree$dbh[1], tree$h[1], 0, h.lim, b1),
                             h.com = h.lim)

      volume <- rbind(volume, volume.i)

    }

  }}

  return(volume)

}


.stem.axis <- function(data, scan.approach = "single"){

  if(scan.approach == "multi"){

    s <- sample(nrow(data), round(nrow(data)*0.25))
    data <- data[s, ]

    } else {

    data <- data[data$prob < 0.1 | data$prob > 0.9, ]}


  if(nrow(data) < 100 | nrow(data) > 1000000 | min(data$z) > 1.3){

    eje <- data.frame(tree = as.numeric(), sec = as.numeric(), x = as.numeric(), y = as.numeric(), n.w.ratio = as.numeric())

    } else {

  eje <- data.frame(tree = unique(data$tree), sec = seq(0, round(max(data$z), 1), by = 0.1))


  dbscan <- dbscan::dbscan(data[, c("x", "z"), drop = FALSE], eps = 0.25)
  data$cluster <- as.factor(dbscan$cluster)
  # plot(data$z, data$x, col = data$cluster, asp =1)

  cluster <- data.frame(table(data$cluster))
  colnames(cluster) <- c("cluster", "freq")
  cluster <- cluster[cluster$freq == max(cluster$freq), ]

  data <- merge(data, cluster, by = "cluster", all = FALSE)

  # plot(data$z, data$x, col = data$cluster, asp =1)

  mod.x <- stats::lm(data = data, x ~ z)
  mod.y <- stats::lm(data = data, y ~ z)

  eje$x <- stats::coef(mod.x)[1] + stats::coef(mod.x)[2] * eje$sec
  eje$y <- stats::coef(mod.y)[1] + stats::coef(mod.y)[2] * eje$sec

  n.w.ratio <- as.numeric(.n.w.ratio(data)[2])

  eje$n.w.ratio <- n.w.ratio

  }

  return(eje)


}


# Select part of point cloud free of low vegetation and crown

.getStem <- function(data){

  k <- stats::density(data$z)

  den <- data.frame(x = k$x, y = k$y)

  n.ini <- sum(den$y)

  n <- sum(den[den$x < den[den$y == max(den$y), ]$x, ]$y)
  if(n < n.ini / 2){
    den <- den[den$x > den[den$y == max(den$y), ]$x, ]
    den <- den[den$y > stats::quantile(den$y, probs = 0.25, na.rm = T), ]} else {
      den <- den[den$x < den[den$y == max(den$y), ]$x, ]}


  den$dev1 <- c(diff(den$y), 0)

  den$dev1 <- abs(den$dev1)
  den <- den[den$dev1 > stats::quantile(den$dev1, probs = 0.75, na.rm = T), ]

  den$diff <- c(diff(den$x), 0)
  den$cut <- ifelse(den$diff > .getmode(den$diff) + 0.01, 1, 0)
  den <- den[den$cut == 1, ]
  den <- den[den$diff == max(den$diff), ]

  return(den[, c("x", "diff")])

}





# Add quotes and paste a vector

.quot.past <- function(x, quot = "'", sep = "", collapse = ", ") {

  paste(quot, trimws(x), quot, sep = sep, collapse = collapse)

}


# Check object class and dimension, and values validity

.check.class <- function(x, x.class, name, col.of = NULL, n = NULL,
                         val = NULL, order.val = TRUE) {

  # Check object class
  if (x.class == "integer")
    .check <- is.numeric(x) && all(round(x[!is.na(x)], 0) == x[!is.na(x)])
  else
    .check <- get(paste("is", x.class, sep = "."))(x)
  if (!.check)
    stop(ifelse(is.null(col.of), .quot.past(name),
                paste("Column", .quot.past(name), "of", .quot.past(col.of))),
         " must be a ", x.class,
         ifelse(x.class %in% c("integer", "numeric", "character", "logical") &
                  !is.null(n) & n > 1, " vector", ""),
         ".")

  # Check object dimension
  if (!is.null(n)) {

    if (!x.class %in% c("data.frame")) {

      if (length(x) == 0)
        stop(ifelse(is.null(col.of), .quot.past(name),
                    paste("Column", .quot.past(name), "of",
                          .quot.past(col.of))),
             " must have length different from zero.")
      if (n == 1) {

        x <- unique(x)
        if (length(x) > 1) {

          warning("Only first value in ",
                  ifelse(is.null(col.of), .quot.past(name),
                         paste("column", .quot.past(name), "of",
                               .quot.past(col.of))),
                  " was taken into account during the execution.",
                  immediate. = TRUE)
          x <- x[1]

        }

      }

    } else {

      if (nrow(x) == 0) stop(.quot.past(name), " must have at least one row.")
      if (n == 1) {

        x <- unique(x)
        if (nrow(x) > 1) {

          warning("Only first row in ", .quot.past(name), " was taken into ",
                  "account during the execution.", immediate. = TRUE)
          x <- x[1, , drop = FALSE]

        }

      }

    }

  }

  # Check values and order
  if (!is.null(val)) {

    .inval <- x[!x %in% val]
    if (length(.inval) > 0)
      stop("Value(s) ", .quot.past(.inval), " in ",
           ifelse(is.null(col.of), .quot.past(name),
                  paste("column", .quot.past(name), "of",
                        .quot.past(col.of))),
           " is(are) not valid.\nValid value(s) is(are) ", .quot.past(val), ".")
    if (order.val)
      x <- val[val %in% x]


  }

  return(x)

}


# Check mandatory columns and their classes for a data.frame

.check.col <- function(x, x.mand, x.class, name, def = NULL,
                       na.action = c("stop", "warning")[1]) {

  # Check mandatory columns existence, assign values by default (if necessary),
  # and remove non necessary columns
  .miss <- colnames(x.mand)[apply(x.mand, 2, any)]
  .miss <- .miss[!.miss %in% colnames(x)]
  if (length(.miss) > 0)
    stop("According to specified arguments, ", .quot.past(name), " must have ",
         "column(s) named: ", .quot.past(.miss), ".")
  if (!is.null(def)) {

    .miss <- colnames(def)[!colnames(def) %in% colnames(x)]
    if (length(.miss) > 0)
      x <- cbind(x, do.call(cbind, lapply(def[, .miss], rep, times = nrow(x))))

  }
  x <- x[, colnames(x.mand), drop = FALSE]

  # Check mandatory columns classes
  for (.i in colnames(x)[colnames(x) %in% names(x.class)])
    x[, .i] <- .check.class(x = x[, .i], x.class = x.class[.i], name = .i,
                            col.of = name, n = nrow(x))

  # Check all values are no-NA
  .col.na <- colnames(x)[apply(is.na(x), 2, any)]
  if (length(.col.na) > 0) {

    if (na.action %in% "stop")
      stop("Column(s) ", .quot.past(.col.na), " of ", .quot.past(name),
           " has(have) missing values.")
    else if (na.action %in% "warning")
      warning("Column(s) ", .quot.past(.col.na),  "of ", .quot.past(name),
              " has(have) missing values.", immediate. = TRUE)

  }

  return(x)

}


# Compute number of decimals places

.decimals <- function(x){

  x <- format(x, scient = FALSE)
  ifelse(!base::grepl(".", x, fix = TRUE), 0,
         nchar(base::strsplit(x, ".", fix = TRUE)[[1]][2]))

}


# Compute radius, k and BAF, and tree variables according to plot design(s) and
# 'tree.var'. Currently available tree variables: basal area (g) and volume (v)

.tree.calc <- function(tree, plot.design, tree.var,
                       v.calc = c("coeff", "parab")[2]) {

  # Create data.frame where results will be saved
  .col.names <- plot.design
  if ("k.tree" %in% names(plot.design))
    .col.names <- c(.col.names,
                    paste("radius", plot.design["k.tree"], sep = "."))
  .col.names <- c(.col.names, tree.var)
  tree <- cbind(tree, matrix(nrow = nrow(tree), ncol = length(.col.names),
                             dimnames = list(NULL, .col.names)))

  # Compute radius (m) for fixed area plot design
  if (plot.design["fixed.area"] %in% colnames(tree))
    tree[, plot.design["fixed.area"]] <- tree[, "h.dist"]

  # Compute k (trees) and associated radius (m) for k-tree design
  if (plot.design["k.tree"] %in% colnames(tree)) {

    # Order trees by horizontal distance
    .ord <- order(tree[, "h.dist"], decreasing = F)
    # Save order as k
    tree[.ord, plot.design["k.tree"]] <- 1:nrow(tree)
    # Compute associated radius
    .dist <- c(tree[.ord[-1], "h.dist"],
               utils::tail(tree[.ord, "h.dist"], n = 1))
    tree[.ord, paste("radius", plot.design["k.tree"], sep = ".")] <-
      (tree[.ord, "h.dist"] + .dist) / 2

  }

  # Compute BAF (m2/ha) threshold for angle-count design
  if (plot.design["angle.count"] %in% colnames(tree))
    tree[, plot.design["angle.count"]] <- 2500 / (tree[, "h.dist"] /
                                                    tree[, "dbh"]) ^ 2

  # Compute basal area (m^2)
  if ("g" %in% colnames(tree)) tree[, "g"] <- (pi / 4) * tree[, "dbh"] ^ 2

  # Compute volume (m^3)
  if ("v" %in% colnames(tree)) tree[, "v"] <- tree[, "v"]

  # Compute comercial volume (m^3)
  if ("v.com" %in% colnames(tree)) tree[, "v.com"] <- tree[, "v.com"]

  # if ("v" %in% colnames(tree)) {
  #
  #   if (v.calc == "coeff") {
  #
  #     # Coefficient of 0.45
  #     tree[, "v"] <- (pi / 4) * tree[, "dbh"] ^ 2 * tree[, "h"] * 0.45
  #
  #   } else if (v.calc == "parab") {
  #
  #     # Paraboloid
  #     tree[, "v"] <- pi * (tree[, "h"] ^ 2 / 2) *
  #       ((tree[, "dbh"] / 2) ^ 2 / (tree[, "h"] - 1.3))
  #
  #   } else
  #     stop("Argument for tree volume calculation must be 'coeff' or 'parab'.")
  #
  # }

  return(tree)

}


# Custom rounding of numbers

.customCeiling <- function(x, Decimals = 1) {

  ceiling(x * 10 ^ Decimals) / 10 ^ Decimals

}

.customFloor <- function(x, Decimals = 1) {

  floor(x * 10 ^ Decimals) / 10 ^ Decimals

}


# Custom format for numbers

.format.numb <- function(x, dec) {

  format(x, trim = TRUE, nsmall = dec)

}


# Compute several weighted mean functions for a numeric vector

.wmean.calculation <- function(data, w, mean.names) {

  sapply(mean.names,
         function(x, data, w) {
           get(paste("weighted_mean", x, sep = "_"))(data, w)
         },
         data = data, w = w)

}


# Compute expansion factors with occlusion corrections and estimate stand
# variables per ha, and compute mean diameters and heights

.stand.calc <- function(val, plot.design, tree, var.metr, ds.meth,
                        var.field.user, mean.d, mean.h) {

  # Select data according to specified radius/k/BAF value, and define radius
  if (names(plot.design) %in% c("fixed.area", "k.tree")) {

    tree <- tree[tree[, plot.design] <= val, , drop = FALSE]
    .radius <- switch(names(plot.design), fixed.area = val,
                      k.tree = tree[nrow(tree), paste("radius", plot.design,
                                                      sep = ".")])

  } else if (names(plot.design) %in% "angle.count")
    tree <- tree[tree[, plot.design] >= val, , drop = FALSE]

  # Add auxiliary column for density calculations
  .col.names <- c(paste("N", c("tls", names(ds.meth), "sh", "pam"), sep = "."),
                  "N")
  .col.names <- .col.names[.col.names %in% var.metr]
  if (length(.col.names) > 0) tree <- cbind(tree, n = 1)

  # Initialize data.frame where results will be saved: id, radius/k/BAF, and
  # variables and/or metrics according to 'var.metr'
  .col.names <- c("id", plot.design, var.metr)
  if ("stratum" %in% colnames(tree)) .col.names <- c("stratum", .col.names)
  stand <- data.frame(matrix(nrow = length(val), ncol = length(.col.names),
                             dimnames = list(names(val), .col.names)),
                      stringsAsFactors = FALSE)

  # Plot ID and radius/k/BAF
  stand[, "id"] <- unique(tree[, "id"])
  stand[, plot.design] <- val


  # Compute expansion factors to convert units per plot into units per ha, ----
  # with occlusion corrections, and estimate stand variables per ha: density
  # (trees/ha), basal area (m2/ha), volume (m3/ha), and optionally volume
  # (m3/ha) and biomass (Mg/ha) derived from trees attributes provided by the
  # user

  # Identify stand variables to be computed, and associated tree variable and
  # expansion factor
  .col.names <- c(sapply(c("N", "G", "V", "V.com", "h.com"), paste,
                         c("tls", names(ds.meth), "sh", "pam"), sep = "."),
                  "N", "G", "V", names(var.field.user))
  .col.names <- matrix("", nrow = length(.col.names), ncol = 2,
                       dimnames = list(.col.names, c("var", "ef")))
  for (.i in c("N", "G", "V"))
    .col.names[c(sapply(.i, paste, c("tls", names(ds.meth), "sh", "pam"),
                        sep = "."), .i), "var"] <-
    switch(.i, N = "n", G = "g", V = "v")

  .col.names[c("V.com.tls", "V.com.hn", "V.com.hr", "V.com.hn.cov", "V.com.hr.cov", "V.com.sh", "V.com.pam"), "var"] <- "v.com"
  .col.names[c("h.com.tls", "h.com.hn", "h.com.hr", "h.com.hn.cov", "h.com.hr.cov", "h.com.sh", "h.com.pam"), "var"] <- "h.com"

  for (.i in names(var.field.user))
    .col.names[.i, "var"] <- var.field.user[[.i]]
  .col.names[c(paste(c("N", "G", "V"), "tls", sep = "."), "N", "G", "V",
               names(var.field.user)), "ef"] <- "EF"

  .col.names["V.com.tls", "ef"] <- "EF"
  .col.names["h.com.tls", "ef"] <- "EF"

  for (.i in c(names(ds.meth), "sh", "pam"))
    .col.names[paste(c("N", "G", "V", "V.com", "h.com"), .i, sep = "."), "ef"] <- paste("EF", .i,
                                                                      sep =".")
  .col.names <- .col.names[rownames(.col.names) %in% colnames(stand), ,
                           drop = FALSE]
  if (nrow(.col.names) > 0) {

    # Compute expansion factor according to plot design (all plot designs)
    if (names(plot.design) %in% c("fixed.area", "k.tree"))
      .EF <- 10000 / (pi * .radius ^ 2)
    else if (names(plot.design) %in% "angle.count")
      .EF <- val / (pi * (tree[, "dbh"] / 2) ^ 2)
    .EF <- matrix(.EF, nrow = nrow(tree), dimnames = list(NULL, "EF"))

    # Compute expansion factor using distance sampling based correction (fixed
    # area and k-tree plot designs)
    .col.names.2 <- ds.meth[paste("EF", names(ds.meth), sep = ".") %in%
                              .col.names[, "ef"]]
    if (length(.col.names.2) > 0) {

      .P <- tree[nrow(tree), .col.names.2, drop = FALSE]
      colnames(.P) <- paste("EF", names(.col.names.2), sep = ".")
      if (nrow(tree) > 1) .P <- apply(.P, 2, rep, times = nrow(tree))
      .EF <- cbind(.EF, .EF[, "EF"] / .P)

    }

    # Compute expansion factor with correction of the shadowing effect (fixed
    # area and k-tree plot designs)
    if ("EF.sh" %in% .col.names[, "ef"]) {

      # Compute shadow area for each tree
      .sh <- ifelse(tree[, "partial.occlusion"] == 0,

                    # Non-occluded trees
                    ((((pi * .radius ^ 2) - (pi * tree[, "h.dist"] ^ 2)) /
                        (2 * pi)) * atan2(tree[, "dbh"], tree[, "h.dist"])) -
                      ((pi * (tree[, "dbh"] / 2) ^ 2) / 2),

                    # Occluded trees
                    ((((pi * .radius ^ 2) - (pi * tree[, "h.dist"] ^ 2)) /
                        (2 * pi)) * (tree[, "wide"])) -
                      (((pi * (tree[, "dbh"] / 2) ^ 2) * tree[, "wide"]) /
                         atan2(tree[, "dbh"], tree[, "h.dist"])))

      # Correction: sh = 0 if tree is not completely inside the plot
      .sh[tree[, "h.dist"] + tree[, "dbh"] / 2 > .radius] <- 0

      # Compute expansion factor
      .sh <- 1 + (sum(.sh, na.rm = TRUE) / (pi * .radius ^ 2))
      .EF <- cbind(.EF, EF.sh = .EF[, "EF"] * .sh)

    }

    # Compute expansion factor with gap probability attenuation correction
    # (angle-count plot design)
    if ("EF.pam" %in% .col.names[, "ef"]) {

      .De <- mean(tree[, "dbh"]) *
        (1 + (stats::sd(tree[, "dbh"], na.rm = TRUE) /
                mean(tree[, "dbh"], na.rm = TRUE)) ^ 2) ^ 0.5
      .t <- (((val / (pi * (tree[, "dbh"] / 2) ^ 2)) / 10000) * .De *
               tree[, "dbh"]) / (2 * sqrt(val / 10000))
      .Ft <- (2 / .t ^ 2) * (1 - exp(-.t) * (1 + .t))
      .EF <- cbind(.EF, EF.pam = .EF[, "EF"] / .Ft)

    }

    # Estimate variables per ha
    stand[, rownames(.col.names)] <-
      apply(tree[, .col.names[, "var"], drop = FALSE] *
              .EF[, .col.names[, "ef"], drop = FALSE], 2, sum)

  }


  # Compute mean diameters and heights ----

  # Compute mean diameters
  .col.names <- sapply(names(mean.d), paste, c("", ".tls"), sep = "")
  .col.names <- .col.names[apply(apply(.col.names, 1:2,  "%in%",
                                       colnames(stand)), 1, any), ]
  if (length(.col.names) > 0)
    stand[, .col.names] <- .wmean.calculation(data = tree[, "dbh"],
                                              w = rep(1, nrow(tree)),
                                              mean.names = mean.d)

  # Compute mean heights
  .col.names <- sapply(names(mean.h), paste, c("", ".tls"), sep = "")
  .col.names <- .col.names[apply(apply(.col.names, 1:2,  "%in%",
                                       colnames(stand)), 1, any), ]
  if (length(.col.names) > 0)
    stand[, .col.names] <- .wmean.calculation(data = tree[, "h"],
                                              w = rep(1, nrow(tree)),
                                              mean.names = mean.h)


  # Compute number of points
  .col.names <- paste("n.pts", c("", ".est", ".red", ".red.est"), sep = "")
  .col.names <- .col.names[.col.names %in% colnames(stand)]
  if (length(.col.names) > 0)
    stand[, .col.names] <- apply(tree[, .col.names, drop = FALSE], 2, sum)

  return(stand)

}


# LiDAR metrics: 'mean', 'max', 'min', 'sd', 'var', 'mode', 'kurtosis',
# 'skewness', 'p.a.mean', 'p.a.mode', 'weibull_b', 'weibull_c'

.getmode <- function(v) {

  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]

}

.c_function <- function(c, media, varianza){

  varianza - (media^2) * (gamma(1 + 2 / c) -
                            (gamma(1 + 1 / c))^2) / (gamma(1 + 1 / c))^2

}

.points.metrics <- function(rho_seq, data, metr) {

  # Restrict data
  data <- data[data[, "z"] > 0.1, , drop = FALSE]
  data <- data[data[, "rho"] > 0.1, , drop = FALSE]
  data <- data[data[, "r"] > 0.1, , drop = FALSE]

  # Compute metrics for coordinate z
  .metr <- lapply(rho_seq,
                  function(rho, data, metr) {

                    .sub <- data[data[, "rho"] <= rho, "z"]

                    .metr <- rep(NA, length(metr))
                    names(.metr) <- metr

                    if (any(c("mean.z", "p.a.mean.z", "weibull_c.z", "weibull_b.z") %in% names(.metr)))
                      .metr["mean.z"] <- mean(.sub)
                    if ("mean.q.z" %in% names(.metr))
                      .metr["mean.q.z"] <- sqrt(mean(.sub ^ 2))
                    if ("mean.g.z" %in% names(.metr))
                      .metr["mean.g.z"] <- exp(mean(log(.sub[.sub>0])))
                    if ("mean.h.z" %in% names(.metr))
                      .metr["mean.h.z"] <- length(.sub) / sum(1 / .sub)
                    if ("median.z" %in% names(.metr))
                      .metr["median.z"] <- stats::median(.sub)
                    if (any(c("mode.z", "p.a.mode.z") %in% names(.metr)))
                      .metr["mode.z"] <- .getmode(.sub)
                    if (any(c("max.z", "weibull_c.z", "weibull_b.z") %in% names(.metr)))
                      .metr["max.z"] <- max(.sub)
                    if (any(c("min.z", "weibull_c.z", "weibull_b.z") %in% names(.metr)))
                      .metr["min.z"] <- min(.sub)
                    if (any(c("var.z", "weibull_c.z", "weibull_b.z") %in% names(.metr)))
                      .metr["var.z"] <- var(.sub)
                    if ("sd.z" %in% names(.metr))
                      .metr["sd.z"] <- sd(.sub)
                    if ("CV.z" %in% names(.metr))
                      .metr["CV.z"] <- .metr["sd.z"] / .metr["mean.z"]
                    if ("D.z" %in% names(.metr))
                      .metr["D.z"] <- .metr["max.z"] - .metr["min.z"]
                    if ("ID.z" %in% names(.metr))
                      .metr["ID.z"] <- stats::quantile(.sub, prob = 0.75, na.rm = T) - stats::quantile(.sub, prob = 0.25, na.rm = T)

                    if ("kurtosis.z" %in% names(.metr))
                      .metr["kurtosis.z"] <- moments::kurtosis(.sub)
                    if ("skewness.z" %in% names(.metr))
                      .metr["skewness.z"] <- moments::skewness(.sub)

                    if ("p.a.mean.z" %in% names(.metr))
                      .metr["p.a.mean.z"] <- mean(.sub > .metr["mean.z"]) * 100
                    if ("p.a.mode.z" %in% names(.metr))
                      .metr["p.a.mode.z"] <- mean(.sub > .metr["mode.z"]) * 100
                    if ("p.a.2m.z" %in% names(.metr))
                      .metr["p.a.2m.z"] <- mean(.sub > 2) * 100

                    if ("p.b.mean.z" %in% names(.metr))
                      .metr["p.b.mean.z"] <- mean(.sub < .metr["mean.z"]) * 100
                    if ("p.b.mode.z" %in% names(.metr))
                      .metr["p.b.mode.z"] <- mean(.sub < .metr["mode.z"]) * 100
                    if ("p.b.2m.z" %in% names(.metr))
                      .metr["p.b.2m.z"] <- mean(.sub < 2) * 100

                    if ("CRR.z" %in% names(.metr))
                      .metr["CRR.z"] <- .metr["mean.z"] / .metr["max.z"]

                    if ("L2.z" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 1) - choose(length(.sub) - .sub, 1)) * .sub
                      .metr["L2.z"] <- (1 / (2 * choose(length(.sub), 2))) * sum(sumatorio)
                    if ("L3.z" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 2) - 2 * choose(.sub - 1, 1) * choose(length(.sub) - .sub, 1) + choose(length(.sub) - .sub, 2)) * .sub
                      .metr["L3.z"] <- (1 / (3 * choose(length(.sub), 3))) * sum(sumatorio)
                    if ("L4.z" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 3) - 3 * choose(.sub - 1, 2) * choose(length(.sub) - .sub, 1) +
                                      3 * choose(.sub - 1, 1) * choose(length(.sub) - .sub, 2) - choose(length(.sub) - .sub, 3)) * .sub
                      .metr["L4.z"] <- (1 / (4 * choose(length(.sub), 4))) * sum(sumatorio)


                    if ("L.CV.z" %in% names(.metr))
                      .metr["L.CV.z"] <- .metr["L2.z"] / .metr["mean.z"]
                    if ("L.skewness.z" %in% names(.metr))
                      .metr["L.skewness.z"] <- .metr["L3.z"] / .metr["L2.z"]
                    if ("L.kurtosis.z" %in% names(.metr))
                      .metr["L.kurtosis.z"] <- .metr["L4.z"] / .metr["L2.z"]


                    if ("median.a.d.z" %in% names(.metr))
                      .metr["median.a.d.z"] <- stats::median(abs(.sub - .metr["mean.z"]))
                    if ("mode.a.d.z" %in% names(.metr))
                      .metr["mode.a.d.z"] <- .getmode(abs(.sub - .metr["mean.z"]))


                    if (any(c("weibull_c.z", "weibull_b.z") %in% names(.metr))) {

                      .error <- try(stats::uniroot(.c_function, media = .metr["mean.z"],
                                                   varianza = .metr["var.z"],
                                                   interval = c(.metr["min.z"],
                                                                .metr["max.z"]))$root, silent = TRUE)

                      if(class(.error)[1] == "try-error"){

                        .metr["weibull_c.z"] <- NA
                        .metr["weibull_b.z"] <- NA

                      } else {

                        .metr["weibull_c.z"] <-
                          stats::uniroot(.c_function, media = .metr["mean.z"],
                                         varianza = .metr["var.z"],
                                         interval = c(.metr["min.z"],
                                                      .metr["max.z"]))$root

                      }

                      if ("weibull_b.z" %in% names(.metr)) {

                        .metr["weibull_b.z"] <-
                          .metr["mean.z"] / gamma(1 + 1 / .metr["weibull_c.z"])}

                    }



                    # Compute metrics for coordinate rho

                    .sub <- data[data[, "rho"] <= rho, "rho"]

                    if (any(c("mean.rho", "p.a.mean.rho", "weibull_c.rho", "weibull_b.rho")  %in% names(.metr)))
                      .metr["mean.rho"] <- mean(.sub)
                    if ("mean.q.rho" %in% names(.metr))
                      .metr["mean.q.rho"] <- sqrt(mean(.sub ^ 2))
                    if ("mean.g.rho" %in% names(.metr))
                      .metr["mean.g.rho"] <- exp(mean(log(.sub[.sub>0])))
                    if ("mean.h.rho" %in% names(.metr))
                      .metr["mean.h.rho"] <- length(.sub) / sum(1 / .sub)
                    if ("median.rho" %in% names(.metr))
                      .metr["median.rho"] <- stats::median(.sub)
                    if (any(c("mode.rho", "p.a.mode.rho") %in% names(.metr)))
                      .metr["mode.rho"] <- .getmode(.sub)
                    if (any(c("max.rho", "weibull_c.rho", "weibull_b.rho") %in% names(.metr)))
                      .metr["max.rho"] <- max(.sub)
                    if (any(c("min.rho", "weibull_c.rho", "weibull_b.rho") %in% names(.metr)))
                      .metr["min.rho"] <- min(.sub)
                    if (any(c("var.rho", "weibull_c.rho", "weibull_b.rho") %in% names(.metr)))
                      .metr["var.rho"] <- var(.sub)
                    if ("sd.rho" %in% names(.metr))
                      .metr["sd.rho"] <- sd(.sub)
                    if ("CV.rho" %in% names(.metr))
                      .metr["CV.rho"] <- .metr["sd.rho"] / .metr["mean.rho"]
                    if ("D.rho" %in% names(.metr))
                      .metr["D.rho"] <- .metr["max.rho"] - .metr["min.rho"]
                    if ("ID.rho" %in% names(.metr))
                      .metr["ID.rho"] <- stats::quantile(.sub, prob = 0.75, na.rm = T) - stats::quantile(.sub, prob = 0.25, na.rm = T)

                    if ("kurtosis.rho" %in% names(.metr))
                      .metr["kurtosis.rho"] <- moments::kurtosis(.sub)
                    if ("skewness.rho" %in% names(.metr))
                      .metr["skewness.rho"] <- moments::skewness(.sub)

                    if ("p.a.mean.rho" %in% names(.metr))
                      .metr["p.a.mean.rho"] <- mean(.sub > .metr["mean.rho"]) * 100
                    if ("p.a.mode.rho" %in% names(.metr))
                      .metr["p.a.mode.rho"] <- mean(.sub > .metr["mode.rho"]) * 100
                    if ("p.b.mean.rho" %in% names(.metr))
                      .metr["p.b.mean.rho"] <- mean(.sub < .metr["mean.rho"]) * 100
                    if ("p.b.mode.rho" %in% names(.metr))
                      .metr["p.b.mode.rho"] <- mean(.sub < .metr["mode.rho"]) * 100

                    if ("CRR.rho" %in% names(.metr))
                      .metr["CRR.rho"] <- .metr["mean.rho"] / .metr["max.rho"]

                    if ("L2.rho" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 1) - choose(length(.sub) - .sub, 1)) * .sub
                    .metr["L2.rho"] <- (1 / (2 * choose(length(.sub), 2))) * sum(sumatorio)
                    if ("L3.rho" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 2) - 2 * choose(.sub - 1, 1) * choose(length(.sub) - .sub, 1) + choose(length(.sub) - .sub, 2)) * .sub
                    .metr["L3.rho"] <- (1 / (3 * choose(length(.sub), 3))) * sum(sumatorio)
                    if ("L4.rho" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 3) - 3 * choose(.sub - 1, 2) * choose(length(.sub) - .sub, 1) +
                                      3 * choose(.sub - 1, 1) * choose(length(.sub) - .sub, 2) - choose(length(.sub) - .sub, 3)) * .sub
                    .metr["L4.rho"] <- (1 / (4 * choose(length(.sub), 4))) * sum(sumatorio)


                    if ("L.CV.rho" %in% names(.metr))
                      .metr["L.CV.rho"] <- .metr["L2.rho"] / .metr["mean.rho"]
                    if ("L.skewness.rho" %in% names(.metr))
                      .metr["L.skewness.rho"] <- .metr["L3.rho"] / .metr["L2.rho"]
                    if ("L.kurtosis.rho" %in% names(.metr))
                      .metr["L.kurtosis.rho"] <- .metr["L4.rho"] / .metr["L2.rho"]

                    if ("median.a.d.rho" %in% names(.metr))
                      .metr["median.a.d.rho"] <- stats::median(abs(.sub - .metr["mean.rho"]))
                    if ("mode.a.d.rho" %in% names(.metr))
                      .metr["mode.a.d.rho"] <- .getmode(abs(.sub - .metr["mean.rho"]))



                    if (any(c("weibull_c.rho", "weibull_b.rho") %in% names(.metr))) {

                      .error <- try(stats::uniroot(.c_function, media = .metr["mean.rho"],
                                                   varianza = .metr["var.rho"],
                                                   interval = c(.metr["min.rho"],
                                                                .metr["max.rho"]))$root, silent = TRUE)

                      if(class(.error)[1] == "try-error"){

                        .metr["weibull_c.rho"] <- NA
                        .metr["weibull_b.rho"] <- NA

                      } else {

                      .metr["weibull_c.rho"] <-
                        stats::uniroot(.c_function, media = .metr["mean.rho"],
                                       varianza = .metr["var.rho"],
                                       interval = c(.metr["min.rho"],
                                                    .metr["max.rho"]))$root

                      }

                    if ("weibull_b.rho" %in% names(.metr)) {

                      .metr["weibull_b.rho"] <-
                        .metr["mean.rho"] / gamma(1 + 1 / .metr["weibull_c.rho"])}

                    }


                    # Compute metrics for coordinate r

                    .sub <- data[data[, "rho"] <= rho, "r"]

                    if (any(c("mean.r", "p.a.mean.r", "weibull_c.r", "weibull_b.r") %in% names(.metr)))
                      .metr["mean.r"] <- mean(.sub)
                    if ("mean.q.r" %in% names(.metr))
                      .metr["mean.q.r"] <- sqrt(mean(.sub ^ 2))
                    if ("mean.g.r" %in% names(.metr))
                      .metr["mean.g.r"] <- exp(mean(log(.sub[.sub>0])))
                    if ("mean.h.r" %in% names(.metr))
                      .metr["mean.h.r"] <- length(.sub) / sum(1 / .sub)
                    if ("median.r" %in% names(.metr))
                      .metr["median.r"] <- stats::median(.sub)
                    if (any(c("mode.r", "p.a.mode.r") %in% names(.metr)))
                      .metr["mode.r"] <- .getmode(.sub)
                    if (any(c("max.r", "weibull_c.r", "weibull_b.r") %in% names(.metr)))
                      .metr["max.r"] <- max(.sub)
                    if (any(c("min.r", "weibull_c.r", "weibull_b.r") %in% names(.metr)))
                      .metr["min.r"] <- min(.sub)
                    if (any(c("var.r", "weibull_c.r", "weibull_b.r") %in% names(.metr)))
                      .metr["var.r"] <- var(.sub)
                    if ("sd.r" %in% names(.metr))
                      .metr["sd.r"] <- sd(.sub)
                    if ("CV.r" %in% names(.metr))
                      .metr["CV.r"] <- .metr["sd.r"] / .metr["mean.r"]
                    if ("D.r" %in% names(.metr))
                      .metr["D.r"] <- .metr["max.r"] - .metr["min.r"]
                    if ("ID.r" %in% names(.metr))
                      .metr["ID.r"] <- stats::quantile(.sub, prob = 0.75, na.rm = T) - stats::quantile(.sub, prob = 0.25, na.rm = T)

                    if ("kurtosis.r" %in% names(.metr))
                      .metr["kurtosis.r"] <- moments::kurtosis(.sub)
                    if ("skewness.r" %in% names(.metr))
                      .metr["skewness.r"] <- moments::skewness(.sub)

                    if ("p.a.mean.r" %in% names(.metr))
                      .metr["p.a.mean.r"] <- mean(.sub > .metr["mean.r"]) * 100
                    if ("p.a.mode.r" %in% names(.metr))
                      .metr["p.a.mode.r"] <- mean(.sub > .metr["mode.r"]) * 100
                    if ("p.b.mean.r" %in% names(.metr))
                      .metr["p.b.mean.r"] <- mean(.sub < .metr["mean.r"]) * 100
                    if ("p.b.mode.r" %in% names(.metr))
                      .metr["p.b.mode.r"] <- mean(.sub < .metr["mode.r"]) * 100

                    if ("CRR.r" %in% names(.metr))
                      .metr["CRR.r"] <- .metr["mean.r"] / .metr["max.r"]

                    if ("L2.r" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 1) - choose(length(.sub) - .sub, 1)) * .sub
                    .metr["L2.r"] <- (1 / (2 * choose(length(.sub), 2))) * sum(sumatorio)
                    if ("L3.r" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 2) - 2 * choose(.sub - 1, 1) * choose(length(.sub) - .sub, 1) + choose(length(.sub) - .sub, 2)) * .sub
                    .metr["L3.r"] <- (1 / (3 * choose(length(.sub), 3))) * sum(sumatorio)
                    if ("L4.r" %in% names(.metr))
                      sumatorio <- (choose(.sub - 1, 3) - 3 * choose(.sub - 1, 2) * choose(length(.sub) - .sub, 1) +
                                      3 * choose(.sub - 1, 1) * choose(length(.sub) - .sub, 2) - choose(length(.sub) - .sub, 3)) * .sub
                    .metr["L4.r"] <- (1 / (4 * choose(length(.sub), 4))) * sum(sumatorio)


                    if ("L.CV.r" %in% names(.metr))
                      .metr["L.CV.r"] <- .metr["L2.r"] / .metr["mean.r"]
                    if ("L.skewness.r" %in% names(.metr))
                      .metr["L.skewness.r"] <- .metr["L3.r"] / .metr["L2.r"]
                    if ("L.kurtosis.r" %in% names(.metr))
                      .metr["L.kurtosis.r"] <- .metr["L4.r"] / .metr["L2.r"]

                    if ("median.a.d.r" %in% names(.metr))
                      .metr["median.a.d.r"] <- stats::median(abs(.sub - .metr["mean.r"]))
                    if ("mode.a.d.r" %in% names(.metr))
                      .metr["mode.a.d.r"] <- .getmode(abs(.sub - .metr["mean.r"]))


                    if (any(c("weibull_c.r", "weibull_b.r") %in% names(.metr))) {

                      .error <- try(stats::uniroot(.c_function, media = .metr["mean.r"],
                                                   varianza = .metr["var.r"],
                                                   interval = c(.metr["min.r"],
                                                                .metr["max.r"]))$root, silent = TRUE)

                      if(class(.error)[1] == "try-error"){

                        .metr["weibull_c.r"] <- NA
                        .metr["weibull_b.r"] <- NA

                      } else {

                        .metr["weibull_c.r"] <-
                          stats::uniroot(.c_function, media = .metr["mean.r"],
                                         varianza = .metr["var.r"],
                                         interval = c(.metr["min.r"],
                                                      .metr["max.r"]))$root

                      }

                      if ("weibull_b.r" %in% names(.metr)) {

                        .metr["weibull_b.r"] <-
                          .metr["mean.r"] / gamma(1 + 1 / .metr["weibull_c.r"])}

                    }

                    return(.metr)

                  },

                  data = data, metr = metr)

  .metr <- do.call(rbind, .metr)

  return(.metr)

}


# Simulate plots and compute stand variables and metrics

.sim.calc <- function(funct, tree.tls, tree.ds, tree.field,
                      plot.design, plot.parameters, scan.approach, var.metr,
                      v.calc, dbh.min, h.min, max.dist, dir.data, save.result,
                      dir.result) {


  # Define available values for arguments and other auxiliary objects ----

  # Define a list containing mandatory trees' database(s) according to 'function
  # case'. Currently available 'function cases': 'sim', 'metr' and 'est'
  if(is.null(tree.field))
    .funct <- list(sim = c(field = "field", tls = "TLS"), metr = c(tls = "TLS"),
                   est = c(tls = "TLS"))

  else {

   .funct <- list(sim = c(field = "field", tls = "TLS"), metr = c(field = "field", tls = "TLS"),
                  est = c(tls = "TLS"))
  }


  # Define a character vector containing index name (radius, k or BAF) for each
  # available plot design. Currently available plot designs: 'fixed.area',
  # 'k.tree' and 'angle.count'
  .plot.design <- c(fixed.area = "radius", k.tree = "k", angle.count = "BAF")

  # Define a list containing mandatory columns in 'tree.tls' argument
  # associated to each available scan approach. Currently available scan
  # approaches: 'single' and 'multi'
  .scan.approach <- list(single = c("phi.left", "phi.right"),
                         multi = character())

  # Define a character vector containing mandatory columns in 'tree.ds'
  # argument associated to available distance sampling methodologies. Currently
  # available methodologies: 'hn' (half normal function), 'hn.cov' (half normal
  # function with dbh as covariate), 'hr' (half rate function) and 'hr.cov'
  # (half rate function with dbh as covariate)
  .ds.meth <- c(hn = "P.hn", hr = "P.hr", hn.cov = "P.hn.cov", hr.cov = "P.hr.cov")

  # Define character vectors containing the mean function to be used for each
  # mean diameter/height computation. Currently available means: arithmetic,
  # quadratic, geometric and harmonic
  .mean.d <- c(d = "arit", dg = "sqrt", dgeom = "geom", dharm = "harm")
  .mean.h <- c(h = "arit", hg = "sqrt", hgeom = "geom", hharm = "harm")

  # Define a numeric vector containing probabilities for height percentiles
  # calculation
  .prob = c(1, 5, 10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90, 95, 99)

  # Define a list containing mandatory columns in 'tree.field' argument
  # associated to each available field variables to be computed from optional
  # trees attributes provided by the user. Currently available field variables:
  # 'V.user' and 'W.user'
  .var.field.user <- list(V.user = "v.user", W.user = "w.user")

  # Define a list containing the available TLS variables and metrics according
  # to 'function case', 'tree.ds' availability, scan approach and plot
  # design(s), and the available field variables according to 'function case'
  # and trees variables provided by user in 'tree.field'.
  # Currently available:
  # - TLS variables and metrics (1): 'N.tls', 'N.hn' (2), 'N.hr' (2),
  #   'N.hn.cov' (2), 'N.hr.cov' (2), 'N.sh' (3), 'N.pam' (4); 'G.tls',
  #   'G.hn' (2), 'G.hr' (2), 'G.hn.cov' (2), 'G.hr.cov' (2), 'G.sh' (3),
  #   'G.pam' (4); 'V.tls', 'V.hn' (2), 'V.hr' (2), 'V.hn.cov' (2),
  #   'V.hr.cov'(2), 'V.sh' (3), 'V.pam' (4); 'd.tls', 'dg.tls', 'dgeom.tls',
  #   'dharm.tls'; 'h.tls', 'hg.tls', 'hgeom.tls', 'hharm.tls'; 'd.0.tls',
  #   'dg.0.tls', 'dgeom.0.tls', 'dharm.0.tls'; 'h.0.tls', 'hg.0.tls',
  #   'hgeom.0.tls', 'hharm.0.tls'; 'n.pts', 'n.pts.est',
  #   'n.pts.red', 'n.pts.red.est'; 'P01', 'P05', 'P10', 'P20', 'P25',
  #   'P30', 'P40', 'P50', 'P60', 'P70', 'P75', 'P80', 'P90', 'P95', 'P99';
  #   'mean', 'max', 'min', 'sd', 'var', 'mode', 'kurtosis', 'skewness',
  #   'p.a.mode', 'p.a.mean', 'weibull_c', 'weibull_b'.
  # - Field variables (5): 'N'; 'G'; 'V', 'V.user' (6); 'W.user' (6); 'd', 'dg',
  #   'dgeom', 'dharm'; 'h', 'hg', 'hgeom', 'hharm'; 'd.0', 'dg.0', 'dgeom.0',
  #   'dharm.0'; 'h.0', 'hg.0', 'hgeom.0', 'hharm.0'.
  # Remarks. (1) For function case 'est' only 'N.tls' and 'G.tls' are computed.
  # (2) Only for circular fixed area and k-tree plot designs if scan approach is
  # 'simple' and 'tree.ds' is not NULL. (3) Only for circular fixed area
  # and k-tree plot design if scan approach is 'simple'. (4) Only for
  # angle-count plot design if scan approach is 'simple'. (5) Only for function
  # case 'sim'. (6) Only if corresponding trees variable ('v.user' for 'V.user',
  # 'w.user' for 'W.user') are provided by user in 'tree.field' argument.
  .var.metr <- list(tls = NULL, field = NULL)
  .var.metr$tls <- c(sapply(if(length(tree.tls$v.com) > 0){
                   c("N", "G", "V", "V.com", "h.com")} else{
                   c("N", "G", "V")}, paste,
                            c("tls", names(.ds.meth), "sh", "pam"), sep = "."),
                     t(sapply(names(c(.mean.d, .mean.h)),
                              paste, c(".tls", ".0.tls"), sep = "")),
                     paste("n.pts", c("", ".est", ".red", ".red.est"),
                           sep = ""),
                     sprintf("P%02i", .prob),

                     # Z coordinate
                     "mean.z", "mean.q.z", "mean.g.z", "mean.h.z", "median.z", "mode.z",
                     "max.z", "min.z", "var.z", "sd.z", "CV.z", "D.z", "ID.z",
                     "kurtosis.z", "skewness.z",
                     "p.a.mean.z", "p.a.mode.z", "p.a.2m.z",
                     "p.b.mean.z", "p.b.mode.z", "p.b.2m.z", "CRR.z",
                     "L2.z", "L3.z", "L4.z", "L.CV.z", "L.skewness.z", "L.kurtosis.z",
                     "median.a.d.z", "mode.a.d.z",
                     "weibull_c.z", "weibull_b.z",

                     # Rho coordinate
                     "mean.rho", "mean.q.rho", "mean.g.rho", "mean.h.rho", "median.rho", "mode.rho",
                     "max.rho", "min.rho", "var.rho", "sd.rho", "CV.rho", "D.rho", "ID.rho",
                     "kurtosis.rho", "skewness.rho",
                     "p.a.mean.rho", "p.a.mode.rho",
                     "p.b.mean.rho", "p.b.mode.rho", "CRR.rho",
                     "L2.rho", "L3.rho", "L4.rho", "L.skewness.rho", "L.kurtosis.rho",
                     "L.CV.rho",
                     "median.a.d.rho", "mode.a.d.rho",
                     "weibull_c.rho", "weibull_b.rho",

                     # R coordinate
                     "mean.r", "mean.q.r", "mean.g.r", "mean.h.r", "median.r", "mode.r",
                     "max.r", "min.r", "var.r", "sd.r", "CV.r", "D.r", "ID.r",
                     "kurtosis.r", "skewness.r",
                     "p.a.mean.r", "p.a.mode.r",
                     "p.b.mean.r", "p.b.mode.r", "CRR.r",
                     "L2.r", "L3.r", "L4.r", "L.CV.r", "L.skewness.r", "L.kurtosis.r",
                     "median.a.d.r", "mode.a.d.r",
                     "weibull_c.r", "weibull_b.r")

  .var.metr$field <- c("N", "G", names(.var.field.user),
                       t(sapply(names(c(.mean.d, .mean.h)), paste, c("", ".0"),
                                sep = "")))

  # Define available calculations for trees volume. Currently available
  # calculations: 'coeff' and 'parab'
  .v.calc <- c("coeff", "parab")

  # Define values by default for certain columns in 'plot.parameters'
  .plot.parameters <- data.frame(radius.incr = 0.1, k.incr = 1, BAF.incr = 0.1,
                                 num.trees = 100)


  # Check arguments 'funct', 'plot.design', 'scan.approach', 'v.calc', ----
  # 'dbh.min', 'h.min', 'max.dist', 'save.result', 'dir.data', 'dir.result'

  # 'funct' must be a character string containing the 'function  case'
  funct <- .check.class(x = funct, x.class = "character", name = "funct",
                        n = 1, val = names(.funct))

  # 'plot.design' must be a character vector containing the plot designs to be
  # computed
  plot.design <- .check.class(x = plot.design, x.class = "character",
                              name = "plot.design", n = Inf,
                              val = names(.plot.design))

  # 'scan.approach' must be a character string indicating the scan approach used
  # for obtaining LAS files
  scan.approach <- .check.class(x = scan.approach, x.class = "character",
                                name = "scan.approach", n = 1,
                                val = names(.scan.approach))

  # 'v.calc' must be a character string indicating how to calculate trees volume
  v.calc <- .check.class(x = v.calc, x.class = "character", name = "v.calc",
                         n = 1, val = .v.calc)

  # 'dbh.min', 'h.min' and 'max.dist' must be positive numeric values
  for (.i in c("dbh.min", "h.min", "max.dist")) {

    # Check object class and dimension
    assign(.i, .check.class(x = get(.i), x.class = "numeric", name = .i, n = 1))

    # Check value is  no-NA and strictly positive
    if (is.na(get(.i)) || get(.i) <= 0)
      stop(.quot.past(.i), " argument must be strictly positive.")

  }

  # 'save.result' must be a logical
  save.result <- .check.class(x = save.result, x.class = "logical",
                              name = "save.result", n = 1)

  # 'dir.data' and 'dir.result' must be NULL or a character string containing
  # the absolute path to a existing directory
  # Remark. If 'dir.data' is NULL, working directory is assigned by default. If
  # 'save.result' is FALSE, 'dir.result' is forced to be NULL; otherwise, if
  # 'dir.result' is NULL, working directory is assigned by default
  for (.i in c("dir.data", "dir.result")) {

    if (.i %in% "dir.result" & !save.result) assign(.i, NULL)
    else {

      if (is.null(get(.i))) assign(.i, getwd())
      else {

        # Check object class and dimension
        assign(.i, .check.class(x = get(.i), x.class = "character", name = .i,
                                n = 1))


        # Check directory existence
        if (!dir.exists(get(.i))) stop(.quot.past(.i), " directory must exist.")

      }

    }

  }


  # Check argument 'var.metr' ----

  # Restrict available TLS variables and metrics according to 'function case',
  # 'tree.ds', scan approach and plot design(s)
  if (!"tls" %in% names(.funct[[funct]])) .var.metr$tls <- NULL
  else {

    if (funct %in% "est")
      .var.metr$tls <- .var.metr$tls[.var.metr$tls %in%
                                       paste(c("N", "G"), "tls", sep = ".")]

    if (is.null(tree.ds))
      .var.metr$tls <- .var.metr$tls[!.var.metr$tls %in%
                                       sapply(c("N", "G", "V"), paste,
                                              names(.ds.meth), sep = ".")]

    if (is.null(tree.ds) & length(tree.tls$v.com) > 0)
      .var.metr$tls <- .var.metr$tls[!.var.metr$tls %in%
                                       sapply(c("N", "G", "V", "V.com", "h.com"), paste,
                                              names(.ds.meth), sep = ".")]

    if (!scan.approach %in% "single")
      .var.metr$tls <- .var.metr$tls[!.var.metr$tls %in%
                                       sapply(c("N", "G", "V"), paste,
                                              c(names(.ds.meth), "sh", "pam"),
                                              sep = ".")]

    if (!scan.approach %in% "single" & length(tree.tls$v.com) > 0)
      .var.metr$tls <- .var.metr$tls[!.var.metr$tls %in%
                                       sapply(c("N", "G", "V", "V.com", "h.com"), paste,
                                              c(names(.ds.meth), "sh", "pam"),
                                              sep = ".")]

    .var.metr$tls <- matrix(TRUE, nrow = length(plot.design),
                            ncol = length(.var.metr$tls),
                            dimnames = list(plot.design, .var.metr$tls))

    if (length(tree.tls$v.com) > 0)
    .var.metr$tls[!rownames(.var.metr$tls) %in% c("fixed.area", "k.tree"),
                  colnames(.var.metr$tls) %in% sapply(c("N", "G", "V", "V.com", "h.com"), paste,
                                                      c(names(.ds.meth), "sh"),
                                                      sep = ".")] <- FALSE

    if (length(tree.tls$v.com) < 1)
    .var.metr$tls[!rownames(.var.metr$tls) %in% c("fixed.area", "k.tree"),
                  colnames(.var.metr$tls) %in% sapply(c("N", "G", "V"), paste,
                                                      c(names(.ds.meth), "sh"),
                                                      sep = ".")] <- FALSE

    if (length(tree.tls$v.com) > 0)
    .var.metr$tls[!rownames(.var.metr$tls) %in% "angle.count",
                  colnames(.var.metr$tls) %in% paste(c("N", "G", "V", "V.com", "h.com"), "pam",
                                                     sep = ".")] <- FALSE

    if (length(tree.tls$v.com) < 1)
    .var.metr$tls[!rownames(.var.metr$tls) %in% "angle.count",
                  colnames(.var.metr$tls) %in% paste(c("N", "G", "V"), "pam",
                                                     sep = ".")] <- FALSE

    .var.metr$tls <- .var.metr$tls[, apply(.var.metr$tls, 2, any), drop = FALSE]

  }

  # Restrict field variables according to 'function case' and trees variables
  # provided by user in 'tree.field'
  if (!"field" %in% names(.funct[[funct]])) .var.metr$field <- NULL
  else {

    if (!is.data.frame(tree.field) ||
        any(!.var.field.user %in% colnames(tree.field)))
      .var.metr$field <-
        .var.metr$field[!.var.metr$field %in%
                          names(.var.field.user)[!.var.field.user %in%
                                                   colnames(tree.field)]]
    .var.metr$field <- matrix(TRUE, nrow = length(plot.design),
                              ncol = length(.var.metr$field),
                              dimnames = list(plot.design, .var.metr$field))

  }

  # 'var.metr' must be a list with one or two elements named 'tls' and 'field'
  var.metr <- .check.class(x = var.metr, x.class = class(.var.metr),
                           name = 'var.metr', n = Inf)
  # if (length(var.metr) != length(.var.metr))
  #   stop("'var.metr' must have length equal to ", length(.var.metr), ".")
  .miss <- names(.var.metr)[! names(.var.metr) %in% names(var.metr)]
  if (length(.miss) > 0)
    stop("Element(s) named ", .quot.past(.miss), " is(are) missing in ",
         "'var.metr' argument.")

  # 'var.metr$tls' must be NULL or a character vector containing the TLS
  # variables and metrics to be computed
  # Remark. If '.var.metr$tls' is NULL, 'var.metr$tls' is forced to be NULL;
  # otherwise, if 'var.metr$tls', all possible TLS variables and metrics
  # are assigned by default
  for (.i in names(var.metr)) {

    if (!.i %in% names(.var.metr)) var.metr[[.i]] <- NULL
    else {

      if (is.null(var.metr[[.i]]))
        var.metr[[.i]] <- colnames(.var.metr[[.i]])
      else
        var.metr[[.i]] <- .check.class(x = var.metr[[.i]],
                                       x.class = "character",
                                       name = paste("var.metr", .i, sep ="$"),
                                       n = Inf, val = colnames(.var.metr[[.i]]))

    }

  }

  # Remove plot design(s) with no metric/variable to be computed
  .inval <- do.call(cbind, .var.metr)[, do.call(c, var.metr), drop = FALSE]
  .inval <- rownames(.inval)[apply(!.inval, 1, all)]
  if (length(.inval) > 0) {

    warning("Plot design(s) ", .quot.past(.inval), " in 'plot.design' ",
            "argument was(were) discarded, since no associated metric or ",
            "variable needs to be computed.", immediate. = TRUE)
    plot.design <- plot.design[!plot.design %in% .inval]

  }


  # Check 'plot.parameters' argument ----

  # Define a logical value indicating if results will be computed by stratum
  # (only for function case 'metr' when 'stratum' column is included in
  # 'plot.parameters' argument)
  .by.stratum <- funct %in% "metr" &
    (is.data.frame(plot.parameters) & "stratum" %in% colnames(plot.parameters))

  # Define a logical vector containing all possible mandatory columns in
  # 'plot.parameters' according to 'function case', '.by.stratum' and plot
  # design(s)
  if (funct %in% c("sim", "est")) {

    .col.mand <- sapply(.plot.design[plot.design], paste, c("max", "incr"),
                        sep = ".")
    if (funct %in% "sim") .col.mand <- c(.col.mand, "num.trees")
    .col.mand <- matrix(FALSE, nrow = 1, ncol = length(.col.mand),
                        dimnames = list(NULL, .col.mand))
    .col.mand[, paste(.plot.design[plot.design], "max", sep = ".")] <- TRUE

  } else {

    .col.mand <- c(.plot.design[plot.design], "num.trees")
    .col.mand <- matrix(FALSE, nrow = 1, ncol = length(.col.mand),
                        dimnames = list(NULL, .col.mand))
    .col.mand[, .plot.design[plot.design]] <- TRUE

  }

  if (.by.stratum) .col.mand <- cbind("stratum" = TRUE, .col.mand)

  # Define a character vector containing (if necessary) the class of all
  # possible mandatory columns in 'plot.parameters'
  .col.class <- rep("numeric", ncol(.col.mand))
  names(.col.class) <- colnames(.col.mand)
  .col.class[names(.col.class) %in% "stratum"] <- NA
  .col.class[names(.col.class) %in%
               c(paste("k", c("max", "incr"), sep = "."), "num.trees")] <-
    "integer"
  .col.class <- .col.class[!is.na(.col.class)]

  # 'plot.parameters' must be a data.frame with at least one row
  plot.parameters <- .check.class(x = plot.parameters, x.class = "data.frame",
                                  name = "plot.parameters",
                                  n = ifelse(funct %in% c("sim", "est"), 1,
                                             Inf))

  # Check mandatory columns existence, assign values by default (if necessary),
  # remove non necessary columns, check mandatory columns classes and check all
  # values are no-NA
  plot.parameters <- .check.col(x = plot.parameters, x.mand = .col.mand,
                                x.class = .col.class, name = "plot.parameters",
                                def = .plot.parameters)

  # Check there are no duplicated strata
  if (.by.stratum) {

    # Define stratum code
    .stratum <- plot.parameters[, "stratum", drop = FALSE]
    .stratum <- cbind(code = apply(.stratum, 1, .quot.past), .stratum)

    # Check duplicated
    .dupl <- unique(.stratum[duplicated(.stratum[, "code"]), "code"])
    if (length(.dupl) > 0)
      stop("'plot.parameters' argument has several different rows for the ",
           "stratum(strata) ", .quot.past(.dupl, quot = ""), ".")

  }

  # Check 'radius', 'radius.max', 'k', 'k.max', 'BAF', 'BAF.max' and/or
  # 'num.trees' are strictly positive
  .inval <- c(sapply(.plot.design, paste, c("", ".max"), sep = ""), "num.trees")
  .inval <- colnames(plot.parameters)[colnames(plot.parameters) %in% .inval]
  .inval <- .inval[!apply(plot.parameters[, .inval, drop = FALSE] > 0, 2, all)]
  if (length(.inval) > 0)
    stop("All values in column(s) ", .quot.past(.inval), " of ",
         "'plot.parameters' argument must be strictly positive.")

  # Check 'radius.incr', 'k.incr" and 'BAF.incr' are positive
  .inval <- paste(.plot.design, "incr", sep = ".")
  .inval <- colnames(plot.parameters)[colnames(plot.parameters) %in% .inval]
  .inval <- .inval[!apply(plot.parameters[, .inval, drop = FALSE] >= 0, 2, all)]
  if (length(.inval) > 0)
    stop("All values in column(s) ", .quot.past(.inval), " of ",
         "'plot.parameters' argument must be positive.")

  # Add number of decimals places to be considered
  for (.i in .plot.design[plot.design]) {

    .dec <- sapply(plot.parameters[, ifelse(funct %in% c("sim", "est"),
                                            paste(.i, "incr", sep = "."), .i)],
                   .decimals)
    .dec <- matrix(.dec, ncol = 1,
                   dimnames = list(NULL, paste(.i, "dec", sep = ".")))
    plot.parameters <- cbind(plot.parameters, .dec)

  }


  # Check 'tree.tls' argument ----

  if (!"tls" %in% names(.funct[[funct]])) tree.tls <- NULL
  else {

    # Define a logical value indicating if results will be computed by stratum
    # (only for function case 'metr' when 'stratum' column is included in
    # 'plot.parameters' argument) or if strata will be differentiated in charts
    # with different colours (only for function case 'est' when 'stratum' column
    # is included in 'tree.tls' argument)
    .by.stratum.2 <- .by.stratum |
      (funct %in% "est" & (is.data.frame(tree.tls) &
                             "stratum" %in% colnames(tree.tls)))

    # Define a logical matrix containing all possible mandatory columns in
    # 'tree.tls' according to 'function case', '.by.stratum.2' and 'var.metr'
    if(length(tree.tls$v.com) < 1)
    .col.mand <- c("id", "file", "tree", unlist(.scan.approach), "h.dist",
                   "dbh", "h", "v",
                   paste("n.pts", c("", ".est", ".red", ".red.est"),
                         sep = ""),
                   "partial.occlusion")

    if(length(tree.tls$v.com) > 0)
      .col.mand <- c("id", "file", "tree", unlist(.scan.approach), "h.dist",
                     "dbh", "h", "h.com", "v", "v.com",
                     paste("n.pts", c("", ".est", ".red", ".red.est"),
                           sep = ""),
                     "partial.occlusion")

    if (.by.stratum.2) .col.mand <- c("stratum", .col.mand)
    .col.mand <- matrix(FALSE, nrow = length(var.metr$tls),
                        ncol = length(.col.mand),
                        dimnames = list(var.metr$tls, .col.mand))

    if(length(tree.tls$v.com) < 1)
    .col.mand[, colnames(.col.mand) %in%
                c("stratum", "id", "tree", "h.dist", "dbh", "h", "v")] <- TRUE

    if(length(tree.tls$v.com) > 0)
    .col.mand[, colnames(.col.mand) %in%
                c("stratum", "id", "tree", "h.dist", "dbh", "h", "h.com", "v", "v.com")] <- TRUE

    .col.mand[rownames(.col.mand) %in%
                c(sprintf("P%02i", .prob),
                  # Z coordinate
                  "mean.z", "mean.q.z", "mean.g.z", "mean.h.z", "median.z", "mode.z",
                  "max.z", "min.z", "var.z", "sd.z", "CV.z", "D.z", "ID.z",
                  "kurtosis.z", "skewness.z",
                  "p.a.mean.z", "p.a.mode.z", "p.a.2m",
                  "p.b.mean.z", "p.b.mode.z", "p.b.2m.z", "CRR.z",
                  "L2.z", "L3.z", "L4.z", "L.CV.z", "L.skewness.z", "L.kurtosis.z",
                  "median.a.d.z", "mode.a.d.z",
                  "weibull_c.z", "weibull_b.z",

                  # Rho coordinate
                  "mean.rho", "mean.q.rho", "mean.g.rho", "mean.h.rho", "median.rho", "mode.rho",
                  "max.rho", "min.rho", "var.rho", "sd.rho", "CV.rho", "D.rho", "ID.rho",
                  "kurtosis.rho", "skewness.rho",
                  "p.a.mean.rho", "p.a.mode.rho",
                  "p.b.mean.rho", "p.b.mode.rho", "CRR.rho",
                  "L2.rho", "L3.rho", "L4.rho", "L.skewness.rho", "L.kurtosis.rho",
                  "L.CV.rho",
                  "median.a.d.rho", "mode.a.d.rho",
                  "weibull_c.rho", "weibull_b.rho",

                  # R coordinate
                  "mean.r", "mean.q.r", "mean.g.r", "mean.h.r", "median.r", "mode.r",
                  "max.r", "min.r", "var.r", "sd.r", "CV.r", "D.r", "ID.r",
                  "kurtosis.r", "skewness.r",
                  "p.a.mean.r", "p.a.mode.r",
                  "p.b.mean.r", "p.b.mode.r", "CRR.r",
                  "L2.r", "L3.r", "L4.r", "L.CV.r", "L.skewness.r", "L.kurtosis.r",
                  "median.a.d.z", "mode.a.d.z",
                  "weibull_c.r", "weibull_b.r"),

              colnames(.col.mand) %in% "file"] <- TRUE

    if(length(tree.tls$v.com) < 1)
    .col.mand[rownames(.col.mand) %in% paste(c("N", "G", "V"), "sh", sep = "."),
              colnames(.col.mand) %in%
                c(unlist(.scan.approach), "partial.occlusion")] <- TRUE

    if(length(tree.tls$v.com) > 0)
    .col.mand[rownames(.col.mand) %in% paste(c("N", "G", "V", "V.com", "h.com"), "sh", sep = "."),
              colnames(.col.mand) %in%
                c(unlist(.scan.approach), "partial.occlusion")] <- TRUE

    for (.i in paste("n.pts", c("", ".est", ".red", ".red.est"), sep = ""))
      .col.mand[rownames(.col.mand) %in% .i,
                colnames(.col.mand) %in% .i] <- TRUE

    .col.mand <- .col.mand[, apply(.col.mand, 2, any), drop = FALSE]

    # Define a character vector containing (if necessary) the class of all
    # possible mandatory columns in 'tree.tls'
    .col.class <- rep("numeric", ncol(.col.mand))
    names(.col.class) <- colnames(.col.mand)
    .col.class[names(.col.class) %in% c("stratum", "id", "file", "tree")] <- NA
    .col.class <- .col.class[!is.na(.col.class)]

    # 'tree.tls' must be a data.frame with at least one row
    tree.tls <- .check.class(x = tree.tls, x.class = "data.frame",
                                  name = "tree.tls", n = Inf)

    # Check mandatory columns existence, remove non necessary columns, check
    # mandatory columns classes and check all values are no-NA
    tree.tls <- .check.col(x = tree.tls, x.mand = .col.mand,
                                x.class = .col.class, name = "tree.tls")

    # Define code for plot IDs
    .id <- unique(tree.tls[, "id", drop = FALSE])
    .id <- cbind(code = apply(.id, 1, .quot.past), .id)

    # Check strata
    if (.by.stratum.2) {

      # Check strata are included in 'plot.parameters' (only for function case
      # 'metr' when 'stratum' column is included in 'plot.parameters' argument)
      if (.by.stratum) {

        .miss <- apply(unique(tree.tls[, "stratum", drop = FALSE]), 1,
                       .quot.past)
        .miss <- .miss[!.miss %in%
                         .stratum[.stratum[, "code"] %in% .miss, "code"]]
        if (length(.miss) > 0)
          stop("Stratum(strata) ", .quot.past(.miss, quot = ""), " in ",
               "'tree.tls' argument is(are) missing in 'plot.parameters'",
               " argument.")

      }

      # Define code for unique pairs ('stratum', 'id')
      .stratum.id <- unique(tree.tls[, c("stratum", "id"), drop = FALSE])
      .stratum.id <- cbind(code = paste("(", apply(.stratum.id, 1, .quot.past),
                                        ")", sep = ""),
                           .stratum.id)

      # Check each plot ID has an unique associated stratum
      .dupl <- unique(.stratum.id[duplicated(.stratum.id[, "id"]), "id"])
      if (length(.dupl) > 0)
        stop("Plot ID(s) ", .quot.past(.dupl), " has(have) several different ",
             "strata associated in 'tree.tls' argument.")

      # Save code for pairs ('stratum', 'id'), and remove from 'tree.tls'
      .stratum.id <- .stratum.id[match(.id[, "id"], .stratum.id[, "id"]), ,
                                 drop = FALSE]
      .id <- cbind(.id, stratum.id.code = .stratum.id[, "code"],
                   stratum = .stratum.id[, "stratum"])
      tree.tls <- tree.tls[, !colnames(tree.tls) %in% "stratum",
                                     drop = FALSE]

    }

    # Check TXT files
    if ("file" %in% colnames(tree.tls)) {

      # Define code for unique pairs ('id', 'file')
      .id.file <- unique(tree.tls[, c("id", "file"), drop = FALSE])
      .id.file <- cbind(code = paste("(", apply(.id.file, 1, .quot.past), ")",
                                     sep = ""),
                        .id.file)

      # Check each plot ID has an unique associated TXT file
      .dupl <- unique(.id.file[duplicated(.id.file[, "id"]), "id"])
      if (length(.dupl) > 0)
        stop("Plot ID(s) ", .quot.past(.dupl), " has(have) several different ",
             "TXT files associated in 'tree.tls' argument.")

      # Check each TXT file has an unique associated plot ID
      .dupl <- unique(.id.file[duplicated(.id.file[, "file"]), "file"])
      if (length(.dupl) > 0)
        stop("TXT file(s) ", .quot.past(.dupl), " has(have) several different ",
             "plot IDs associated in 'tree.tls' argument.")

      # Check TXT files existence
      .miss <- .id.file[, "file"]
      .miss <- .miss[!file.exists(file.path(dir.data, .miss))]
      if (length(.miss) > 0)
        stop("TXT file(s) ", .quot.past(.miss), " in 'tree.tls' argument ",
             "does(do) not exists in the directory specified in 'dir.data' ",
             "argument.")

      # Save code for pairs ('id', 'file'), and remove file from 'tree.tls'
      .id.file <- .id.file[match(.id[, "id"], .id.file[, "id"]), , drop = FALSE]
      .id <- cbind(.id, id.file.code = .id.file[, "code"],
                   file = .id.file[, "file"])
      tree.tls <- tree.tls[, !colnames(tree.tls) %in% "file",
                                     drop = FALSE]

    }

    # Define code for pairs ('id', 'tree') and check they are not duplicated
    .id.tree <- tree.tls[, c("id", "tree"), drop = FALSE]
    .id.tree <- cbind(code = paste("(", apply(.id.tree, 1, .quot.past), ")",
                                   sep = ""),
                      .id.tree)
    .dupl <- unique(.id.tree[duplicated(.id.tree[, "code"]), "code"])
    if (length(.dupl) > 0)
      stop("'tree.tls' argument has several different rows for the ",
           "pair(s) of plot ID and tree ", .quot.past(.dupl, quot = ""), ".")

    # Check 'dbh, 'h', 'n.pts', 'n.pts.red', 'n.pts.est' and/or
    # 'n.pts.red.est' are strictly positive
    .inval <- c("dbh", "h", paste("n.pts",
                                  c("", ".est", ".red", ".red.est"), sep = ""))
    .inval <- colnames(tree.tls)[colnames(tree.tls) %in% .inval]
    .inval <- .inval[!apply(tree.tls[, .inval, drop = FALSE] > 0, 2, all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.tls' argument must be strictly positive.")

    # Check 'n.pts.red'/'n.pts.red.est' is greater or equal to
    # 'n.pts'/'n.pts.est'
    .inval <- t(sapply(paste("n.pts", c("", ".red"), sep = ""), paste,
                       c("", ".est"), sep = ""))
    .inval <- .inval[, apply(apply(.inval, 1:2, "%in%",
                                   colnames(tree.tls)), 2, all),
                     drop = FALSE]
    if (ncol(.inval) > 0) {

      .inval <- .inval[, apply(.inval, 2,
                               function(dat.par, dat) {
                                 any(dat[, dat.par[1]] < dat[, dat.par[2]])
                               },
                               dat = tree.tls),
                       drop = FALSE]
      if (ncol(.inval) > 0)
        stop(paste("All values in column(s) ", .quot.past(.inval[1, ]), " of ",
                   "'tree.tls' argument must be greater or equal to the ",
                   "corresponding ones in column(s) ", .quot.past(.inval[2, ]),
                   ". ", sep = ""))

    }

    # Check 'h.dist' is positive
    .inval <- "h.dist"
    .inval <- colnames(tree.tls)[colnames(tree.tls) %in% .inval]
    .inval <- .inval[!apply(tree.tls[, .inval, drop = FALSE] >= 0, 2, all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.tls' argument must be positive.")

    # Check 'phi.left' and 'phi.right' belong to the interval [0, 2 * pi]
    .inval <- unlist(.scan.approach)
    .inval <- colnames(tree.tls)[colnames(tree.tls) %in% .inval]
    .inval <- .inval[!apply(tree.tls[, .inval, drop = FALSE] >= 0 &
                              tree.tls[, .inval, drop = FALSE] <= 2 * pi,
                            2, all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.tls' argument must be in the interval [0, 2 * pi].")

    # Check 'partial.occlusion' is 0 or 1
    .inval <- "partial.occlusion"
    .inval <- colnames(tree.tls)[colnames(tree.tls) %in% .inval]
    .inval <- .inval[!apply(tree.tls[, .inval, drop = FALSE] == 0 |
                              tree.tls[, .inval, drop = FALSE] == 1,
                            2, all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.tls' argument must be 0 or 1.")

  }


  # Check 'tree.ds' argument ----

  if (!"tls" %in% names(.funct[[funct]])) tree.ds <- NULL
  else if (!is.null(tree.ds)) {

    # Define a logical matrix containing all possible mandatory columns in
    # 'tree.ds' according to 'var.metr'
    .col.mand <- c("id", "tree", .ds.meth)
    .col.mand <- matrix(FALSE, nrow = length(var.metr$tls),
                        ncol = length(.col.mand),
                        dimnames = list(var.metr$tls, .col.mand))
    for (.i in names(.ds.meth))
      .col.mand[rownames(.col.mand) %in% paste(c("N", "G", "V", "V.com", "h.com"), .i, sep = "."),
                colnames(.col.mand) %in%
                  c("id", "tree", paste("P", .i, sep = "."))] <- TRUE
    .col.mand <- .col.mand[, apply(.col.mand, 2, any), drop = FALSE]

    # Define a character vector containing (if necessary) the class of all
    # possible mandatory columns in 'tree.ds'
    .col.class <- rep("numeric", ncol(.col.mand))
    names(.col.class) <- colnames(.col.mand)
    .col.class[names(.col.class) %in% c("id", "tree")] <- NA
    .col.class <- .col.class[!is.na(.col.class)]


    # If there are no mandatory columns, 'tree.ds' is forced to be NULL;
    # otherwise, checking process is continued
    if (ncol(.col.mand) == 0) tree.ds <- NULL
    else {

      # 'tree.ds' must be a list with at least an element named 'tree'
      tree.ds <- .check.class(x = tree.ds, x.class = "list",
                                   name = "tree.ds", n = Inf)
      .miss <- c("tree")[!c("tree") %in% names(tree.ds)]
      if (length(.miss) > 0)
        stop("Element(s) named ", .quot.past(.miss), " is(are) missing in ",
             "'tree.ds' argument.")
      tree.ds <- tree.ds$tree

      # 'tree.ds$tree' must be a data.frame with at least one row
      tree.ds <- .check.class(x = tree.ds, x.class = "data.frame",
                                   name = "tree.ds$tree", n = Inf)

      # Check mandatory columns existence, remove non necessary columns, check
      # mandatory columns classes and check all values are no-NA
      tree.ds <- .check.col(x = tree.ds, x.mand = .col.mand,
                                 x.class = .col.class,
                                 name = "tree.ds$tree")

      # Discard rows corresponding to plot IDs not included in 'tree.tls'
      .miss <- unique(tree.ds[!tree.ds[, "id"] %in% .id[, "id"],
                                   "id"])
      if (length(.miss) > 0) {

        warning("Row(s) of 'tree.ds$tree' argument associated to plot ",
                "ID(s) ", .quot.past(.miss, quot = ""), " was(were) discarded ",
                "because this(these) plot(s) is(are) missing in ",
                "'tree.tls' argument.", immediate. = TRUE)
        tree.ds <- tree.ds[!tree.ds[, "id"] %in% .miss, ,
                                     drop = FALSE]

      }

      # Check plot IDs included in 'tree.tls' are no missing
      .miss <- .id[, "id"]
      .miss <- .miss[!.miss %in% unique(tree.ds[, "id"])]
      if (length(.miss) > 0)
        stop("Plot ID(s) ", .quot.past(.miss), " in 'tree.tls' ",
             "argument is(are) missing in 'tree.ds$tree' argument.")

      # Check there are no duplicated pairs ('id', 'tree')
      .dupl <- tree.ds[, c("id", "tree"), drop = FALSE]
      .dupl <- paste("(", apply(.dupl, 1, .quot.past), ")", sep = "")
      .dupl <- unique(.dupl[duplicated(.dupl)])
      if (length(.dupl) > 0)
        stop("'tree.ds$tree' argument has several different rows for the ",
             "pair(s) of plot ID and tree ", .quot.past(.dupl, quot = ""), ".")

      # Check pairs ('id', 'tree') in 'tree.tls' are included in
      # 'tree.ds'
      .miss <- tree.ds[, c("id", "tree"), drop = FALSE]
      .miss <- paste("(", apply(.miss, 1, .quot.past), ")", sep = "")
      .miss <- .id.tree[, "code"][!.id.tree[, "code"] %in% .miss]
      if (length(.miss) > 0)
        stop("Pair(s) of plot ID and tree ", .quot.past(.miss, quot = ""),
             " in 'tree.tls' argument is(are) missing in ",
             "'tree.ds$tree' argument.")

      # Check 'P.hn', 'P.hr', 'P.hn.cov' and 'P.hr.cov' belong to the interval
      # [0, 1]
      .inval <- .ds.meth
      .inval <- colnames(tree.ds)[colnames(tree.ds) %in% .inval]
      .inval <- .inval[!apply(tree.ds[, .inval, drop = FALSE] >= 0 &
                                tree.ds[, .inval, drop = FALSE] <= 1,
                              2, all)]
      if (length(.inval) > 0)
        stop("All values in column(s) ", .quot.past(.inval), " of ",
             "'tree.ds$tree' argument must be in the interval [0, 1].")

    }

  }


  # Check 'tree.field' argument ----

  if (!"field" %in% names(.funct[[funct]])) tree.field <- NULL
  else {

    # Define a logical matrix containing all possible mandatory columns in
    # 'tree.field' according to 'var.metr'
    .col.mand <- c("id", "tree", "h.dist", "dbh", "h", unlist(.var.field.user))
    .col.mand <- matrix(FALSE, nrow = length(var.metr$field),
                        ncol = length(.col.mand),
                        dimnames = list(var.metr$field, .col.mand))
    .col.mand[, colnames(.col.mand) %in%
                c("id", "tree", "h.dist", "dbh", "h")] <- TRUE
    for (.i in names(.var.field.user))
      .col.mand[rownames(.col.mand) %in% .i,
                colnames(.col.mand) %in% .var.field.user[[.i]]] <- TRUE
    .col.mand <- .col.mand[, apply(.col.mand, 2, any), drop = FALSE]

    # Define a character vector containing (if necessary) the class of all
    # possible mandatory columns in 'tree.field'
    .col.class <- rep("numeric", ncol(.col.mand))
    names(.col.class) <- colnames(.col.mand)
    .col.class[names(.col.class) %in% c("id", "tree")] <- NA
    .col.class <- .col.class[!is.na(.col.class)]

    # 'tree.field' must be a data.frame with at least one row
    tree.field <- .check.class(x = tree.field, x.class = "data.frame",
                                    name = "tree.field", n = Inf)

    # Check mandatory columns existence, remove non necessary columns, check
    # mandatory columns classes and check all values are no-NA
    tree.field <- .check.col(x = tree.field, x.mand = .col.mand,
                                  x.class = .col.class,
                                  name = "tree.field")

    # Discard rows corresponding to plot IDs not included in 'tree.tls'
    .miss <- unique(tree.field[!tree.field[, "id"] %in% .id[, "id"],
                                    "id"])
    if (length(.miss) > 0) {

      warning("Row(s) of 'tree.field' argument associated to plot ID(s) ",
              .quot.past(.miss, quot = ""), " was(were) discarded because ",
              "this(these) plot(s) is(are) missing in 'tree.tls' ",
              "argument.", immediate. = TRUE)
      tree.field <- tree.field[!tree.field[, "id"] %in% .miss, ,
                                         drop = FALSE]

    }

    # Check plot IDs included in 'tree.tls' are no missing
    .miss <- .id[, "id"]
    .miss <- .miss[!.miss %in% unique(tree.field[, "id"])]
    if (length(.miss) > 0)
      stop("Plot ID(s) ", .quot.past(.miss), " in 'tree.tls' ",
           "argument is(are) missing in 'tree.field' argument.")

    # Check there are no duplicated pairs ('id', 'tree')
    .dupl <- tree.field[, c("id", "tree"), drop = FALSE]
    .dupl <- paste("(", apply(.dupl, 1, .quot.past), ")", sep = "")
    .dupl <- unique(.dupl[duplicated(.dupl)])
    if (length(.dupl) > 0)
      stop("'tree.field' argument has several different rows for the ",
           "pair(s) of plot ID and tree ", .quot.past(.dupl, quot = ""), ".")

    # Check 'dbh, 'h', 'v.user' and/or 'w.user' are strictly positive
    .inval <- c("dbh", "h", unlist(.var.field.user))
    .inval <- colnames(tree.field)[colnames(tree.field) %in% .inval]
    .inval <- .inval[!apply(tree.field[, .inval, drop = FALSE] > 0, 2,
                            all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.field' argument must be strictly positive.")

    # Check 'h.dist' is positive
    .inval <- "h.dist"
    .inval <- colnames(tree.field)[colnames(tree.field) %in% .inval]
    .inval <- .inval[!apply(tree.field[, .inval, drop = FALSE] >= 0, 2,
                            all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.field' argument must be positive.")

  }


  # Create a list containing empty data.frames where results will be saved ----
  # and restrict trees' database(s) according to 'dbh.min', 'h.min' and/or
  # 'max.dist'

  # Create a list containing empty data.frames where results will be saved for
  # each plot design
  stand <- vector("list", length(plot.design))
  names(stand) <- plot.design
  for (.i in names(stand)) {

    .col.names <- c("id", .plot.design[.i])
    for (.j in names(.funct[[funct]]))
      .col.names <- c(.col.names,
                      var.metr[[.j]][.var.metr[[.j]][.i, var.metr[[.j]]]])
    if (.by.stratum.2) .col.names <- c("stratum", .col.names)
    stand[[.i]] <- data.frame(matrix(numeric(), ncol = length(.col.names),
                                     dimnames = list(NULL, .col.names)),
                              stringsAsFactors = FALSE)

  }

  # Restrict trees' database(s) according to 'dbh.min', 'h.min' and/or
  # 'max.dist'
  for (.i in names(.funct[[funct]])) {

    .tree <- get(paste("tree", .i, sep ="."))
    .inval <- which(.tree[, "dbh"] < dbh.min | .tree[, "h"] < h.min |
                      .tree[, "h.dist"] > max.dist)
    if (length(.inval) > 0) {

      warning(length(.inval), " tree(s) in 'tree.list.", .i, "' argument ",
              "was(were) discarded, since dbh or height were less than minima ",
              "values indicated in 'dbh.min' and 'h.min' arguments, or ",
              "horizontal distance is greater than 'max.dist' argument.",
              immediate. = TRUE)
      .tree <- .tree[- .inval, , drop = FALSE]
      assign(paste("tree.list", .i, sep ="."), .tree)

      # Remove plot ID(s) without trees for calculations below
      .miss <- .id[!.id[, "id"] %in% unique(.tree[, "id"]), "id"]
      if (length(.miss) > 0) {

        if (length(.miss) < nrow(.id)) {

          warning("Plot ID(s) ", .quot.past(.miss), " in 'tree.list.", .i,
                  "' argument was(were) discarded, since all its(their) trees ",
                  "have dbh or height less than minima values indicated in ",
                  "'dbh.min' and 'h.min' arguments, or horizontal distance ",
                  "greater than 'max.dist' argument.", immediate. = TRUE)
          .id <- .id[!.id[, "id"] %in% .miss, , drop = FALSE]

        } else
          stop("All plot ID(s) in 'tree.list.", .i, "' argument was(were) ",
               "discarded, since all trees have dbh or height less than ",
               "minima values indicated in 'dbh.min' and 'h.min' arguments, ",
               "or horizontal distance greater than 'max.dist' argument.")

      }

    }

  }


  # Loop for each TLS plot ----

  for (.i in 1:nrow(.id)) {

    # Define initial time, and print message with the plot ID
    t0 <- Sys.time()
    message("Computing ",
            ifelse(funct %in% c("sim", "est"), "simulations", "metrics"),
            " for plot: ", .id[.i, "code"])

    # Define list where results associated to the plot will be saved
    .stand <- vector("list", length(.funct[[funct]]))
    names(.stand) <- names(.funct[[funct]])

    # Select plot parameters according to '.by.stratum'
    .str <- ifelse(.by.stratum,
                   match(.id[.i, "stratum"], plot.parameters[, "stratum"]), 1)
    .par <- plot.parameters[.str, , drop = FALSE]


    # Loop for each trees' database (TLS and field data) ----

    for (.j in names(.stand)) {

      # Define list where results associated to the plot will be saved
      .stand[[.j]] <- vector("list", length(plot.design))
      names(.stand[[.j]]) <- plot.design


      # Create trees' database (TLS and field data) ----

      # Select data corresponding to the plot from the complete trees' database
      .tree <- get(paste("tree", .j, sep ="."))
      .tree <- .tree[.tree[, "id"] %in% .id[.i, "id"], , drop = FALSE]
      rownames(.tree) <- NULL

      # Convert dbh (cm) to International System of Units (m)
      .tree[, "dbh"] <- .tree[, "dbh"] / 100

      # Compute radius, k and BAF, and tree variables according to plot
      # design(s) and 'tree.var'. Currently available tree variables: basal area
      # (g) and volume (v)
      .col.mand <- c("g", "v", "v.com", "h.com")

      .col.mand <- matrix(FALSE, nrow = length(var.metr[[.j]]),
                          ncol = length(.col.mand),
                          dimnames = list(var.metr[[.j]], .col.mand))
      .col.mand[rownames(.col.mand) %in%
                  c(paste("G", c("tls", names(.ds.meth), "sh", "pam"),
                          sep = "."), "G"),
                colnames(.col.mand) %in% "g"] <- TRUE
      .col.mand[rownames(.col.mand) %in%
                  c(paste("V", c("tls", names(.ds.meth), "sh", "pam"),
                          sep = "."), "V"),
                colnames(.col.mand) %in% "v"] <- TRUE
      .col.mand[rownames(.col.mand) %in%
                  c(paste("V.com", c("tls", names(.ds.meth), "sh", "pam"),
                          sep = "."), "V.com"),
                colnames(.col.mand) %in% "v.com"] <- TRUE
      .col.mand[rownames(.col.mand) %in%
                  c(paste("h.com", c("tls", names(.ds.meth), "sh", "pam"),
                          sep = "."), "h.com"),
                colnames(.col.mand) %in% "v.com"] <- TRUE
      .col.mand <- colnames(.col.mand)[apply(.col.mand, 2, any)]
      .tree <- .tree.calc(tree = .tree, plot.design = .plot.design[plot.design],
                          tree.var = .col.mand, v.calc = v.calc)

      # Compute angular aperture (TLS data)
      if (.j == "tls" & all(c("phi.left", "phi.right") %in% colnames(.tree))) {

        # Compute angular aperture
        .wide <- .tree[, "phi.right"] - .tree[, "phi.left"]
        .wide <- ifelse(.wide < 0, (2 * pi) + .wide, .wide)

        # Remove angular coordinates, and add angular aperture
        .tree <- .tree[, !colnames(.tree) %in% c("phi.left", "phi.right"),
                       drop = FALSE]
        .tree <- cbind(.tree, wide = .wide)

      }

      # Add detection probability from distance sampling database (TLS data)
      if (.j == "tls" & !is.null(tree.ds)) {

        .ds <- tree.ds[tree.ds[, "id"] %in% .id[.i, "id"], ,
                            drop = FALSE]
        .tree <- cbind(.tree, .ds[match(.tree[, "tree"], .ds[, "tree"]),
                                  colnames(.ds) %in% .ds.meth, drop = FALSE])

      }


      # Create points database(s) (TLS data) ----

      .data.tls <- NULL
      if (.j == "tls" & "file" %in% colnames(.id)) {

        # Read points' database
        .data.tls <-
          suppressMessages(vroom::vroom(file.path(dir.data, .id[.i, "file"]),
                                        col_select = c("x", "y", "z", "rho", "r"),
                                        progress = FALSE))
        .data.tls <- as.matrix(.data.tls)

        # Discard points according to 'max.dist'
        .inval <- which(.data.tls[, "rho"] > max.dist)
        if (length(.inval) > 0) {

          warning(length(.inval), " point(s) in point cloud associated to ",
                  "plot ", .id[.i, "code"], " was(were) discarded, since ",
                  "horizontal distance is greater than 'max.dist' argument.",
                  immediate. = TRUE)
          .data.tls <- .data.tls[- .inval, , drop = FALSE]

        }

        # Order .data.tls by rho, select columns required for calculations
        # below, and convert to matrix
        .data.tls <- .data.tls[order(.data.tls[, "rho"], decreasing = FALSE),
                               c("z", "rho", "r"), drop = FALSE]

      }


      # Loop for each plot design (TLS and field data) ----

      for (.k in plot.design) {

        # Order trees' database by radius/k/BAF in order to simplify
        # calculations below
        if (.k %in% c("fixed.area", "k.tree"))
          .tree <- .tree[order(.tree[, .plot.design[.k]], decreasing = FALSE), ,
                         drop = FALSE]
        else if (.k %in% "angle.count")
          .tree <- .tree[order(.tree[, .plot.design[.k]], decreasing = TRUE), ,
                         drop = FALSE]


        # Create radius/k/BAF sequence according to 'function case', ----
        # 'plot.parameters' argument, and radius/k/BAF values in trees' database

        # Select number of decimals
        .par.dec <- .par[, paste(.plot.design[.k], "dec", sep = ".")]

        # Compute range of values in trees' database
        .funct.rang <- switch(.k, fixed.area = ".customCeiling",
                              k.tree = ".customCeiling",
                              angle.count = ".customFloor")
        .db.rang <- range(get(.funct.rang)(.tree[, .plot.design[.k]],
                                           Decimals = .par.dec))

        # Select minimum value (ensuring at least one tree is always selected
        # for fixed area and k-tree plot designs)
        .par.min <- .db.rang[1]

        # Select maximum value
        if (funct %in% c("sim", "est"))
          .par.max <- .par[, paste(.plot.design[.k], "max", sep = ".")]
        else .par.max <- .par[, .plot.design[.k]]
        # Adjust maximum value according to 'max.dist' for fixed area plot
        # design, number of trees for k-tree plot design, and maximum BAF for
        # angle-count plot design
        .ref.max <- switch(.k, fixed.area = max.dist, k.tree = .db.rang[2],
                           angle.count = .db.rang[2])
        if (.par.max > .ref.max) {

          .par.max <- .ref.max
          warning("For ", .funct[[funct]][.j], " plot ", .id[.i, "code"], ", ",
                  ifelse(funct %in% c("sim", "est"), "maximum ", ""),
                  .quot.past(.plot.design[.k]), " was reduced to ",
                  .par.max,
                  switch(.k,
                         fixed.area = paste(", since this is the maximum",
                                            "horizontal distance indicated in",
                                            "'max.dist' argument"),
                         k.tree = paste(", since this is the number of trees",
                                        "in the plot"),
                         angle.count = paste(" to ensure that at least one",
                                             "tree is always selected")),
                  ".", immediate. = TRUE)

        }

        # Adjust maximum value according to minimum one
        if (.par.min > .par.max) {

          .par.max <- .par.min
          warning("For ", .funct[[funct]][.j], " plot ", .id[.i, "code"], ", ",
                  ifelse(funct %in% c("sim", "est"), "maximum ", ""),
                  .quot.past(.plot.design[.k]), " was increased to ", .par.max,
                  switch(.k,
                         fixed.area = paste(" to ensure that at least one",
                                            "tree is always selected"),
                         k.tree = paste(" to ensure that at least one",
                                        "tree is always selected"),
                         angle.count = paste(", since this ensures all trees",
                                             "in the plot are selected")),
                  ".", immediate. = TRUE)

        }

        # Force coincidence between minimum and maximum value, and select
        # increment according to 'function case'
        if (!funct %in% c("sim", "est")) {

          .par.min <- .par.max
          .par.incr <- 0

        } else .par.incr <- .par[, paste(.plot.design[.k], "incr", sep = ".")]

        # Create sequence
        .par.seq <- seq(from = .par.max, to = .par.min, by = - .par.incr)
        .par.seq <- sort(unique(round(.par.seq, .par.dec)))
        names(.par.seq) <- .format.numb(x = .par.seq, dec = .par.dec)


        # Compute stand variables and/or metrics according to plot design  ----
        # and 'var.metr'

        # Compute expansion factors with occlusion corrections and estimate
        # stand variables per ha, and compute mean diameters and heights
        .col.names <- var.metr[[.j]][.var.metr[[.j]][.k, var.metr[[.j]]]]
        .col.names <-
          .col.names[.col.names %in%
                       c(sapply(c("N", "G", "V", "V.com", "h.com"), paste,
                                c("tls", names(.ds.meth), "sh", "pam"),
                                sep = "."),
                         paste(names(c(.mean.d, .mean.h)), "tls", sep = "."),
                         paste("n.pts", c("", ".est", ".red", ".red.est"),
                               sep = ""),
                         "N", "G", "V", names(.var.field.user),
                         names(c(.mean.d, .mean.h)))]
        .stand[[.j]][[.k]] <- lapply(.par.seq, .stand.calc,
                                     plot.design = .plot.design[.k],
                                     tree = .tree, var.metr = .col.names,
                                     ds.meth = .ds.meth,
                                     var.field.user = .var.field.user,
                                     mean.d = .mean.d, mean.h = .mean.h)
        .stand[[.j]][[.k]] <- do.call(rbind, .stand[[.j]][[.k]])

        # Compute mean dominant diameters and heights
        .col.names <- paste(names(c(.mean.d, .mean.h)), ".0", sep = "")
        if (.j %in% "tls") .col.names <- paste(.col.names, .j, sep = ".")
        names(.col.names) <- paste(names(c(.mean.d, .mean.h)), ".0", sep = "")
        .col.names <-
          .col.names[.col.names %in%
                       var.metr[[.j]][.var.metr[[.j]][.k, var.metr[[.j]]]]]
        if (length(.col.names) > 0) {

          if (.k %in% "fixed.area")
            .dh.0 <- fixed_area_cpp(radius_seq = .par.seq,
                                    hdist = .tree[, "h.dist"],
                                    d = .tree[, "dbh"], h = .tree[, "h"],
                                    num = .par[, "num.trees"])
          else if (.k %in% "k.tree")
            .dh.0 <- k_tree_cpp(k_seq = .par.seq,
                                radius_seq = .tree[, "radius.k"],
                                k = .tree[, "k"], d = .tree[, "dbh"],
                                h = .tree[, "h"], num = .par[, "num.trees"])
          else if (.k %in% "angle.count")
            .dh.0 <- angle_count_cpp(baf_seq = .par.seq, baf = .tree[, "BAF"],
                                     d = .tree[, "dbh"], h = .tree[, "h"],
                                     num = .par[, "num.trees"])
          rownames(.dh.0) <- .format.numb(x = .dh.0[, .plot.design[.k]],
                                          dec = .par.dec)
          colnames(.dh.0)[colnames(.dh.0) %in% names(.col.names)] <-
            .col.names[colnames(.dh.0)[colnames(.dh.0) %in% names(.col.names)]]
          .stand[[.j]][[.k]] <- cbind(.stand[[.j]][[.k]],
                                      .dh.0[rownames(.stand[[.j]][[.k]]),
                                            .col.names, drop = FALSE])

        }

        # Compute rho sequence for z coordinate metrics
        if (.k %in% "fixed.area") .rho.seq <- .par.seq
        else if (.k %in% "k.tree")
          .rho.seq <- sapply(.par.seq,
                             function(k, tree) {
                               max(tree[tree[, "k"] <= k, "radius.k"])
                             },
                             tree = .tree)
        else if (.k %in% "angle.count")
          .rho.seq <- sapply(.par.seq,
                             function(BAF, tree) {
                               max(tree[tree[, "BAF"] >= BAF, "h.dist"])
                             },
                             tree = .tree)

        # Compute percentiles of z coordinate
        .col.names <- sprintf("P%02i", .prob)
        .col.names <-
          .col.names[.col.names %in%
                       var.metr[[.j]][.var.metr[[.j]][.k, var.metr[[.j]]]]]
        if (length(.col.names) > 0) {

          .perc <- height_perc_cpp(rho_seq = .rho.seq, z = .data.tls[, "z"],
                                   rho = .data.tls[, "rho"])
          rownames(.perc) <- names(.rho.seq)
          .stand[[.j]][[.k]] <- cbind(.stand[[.j]][[.k]],
                                      .perc[rownames(.stand[[.j]][[.k]]),
                                            .col.names, drop = FALSE])

        }

        # Compute descriptive statistics of z coordinate, percentage of points
        # above mode and mean, and parameters of fitted Weibull distribution
        .col.names <- c("mean.z", "mean.q.z", "mean.g.z", "mean.h.z", "median.z", "mode.z",
                        "max.z", "min.z", "var.z", "sd.z", "CV.z", "D.z", "ID.z",
                        "kurtosis.z", "skewness.z",
                        "p.a.mean.z", "p.a.mode.z", "p.a.2m.z",
                        "p.b.mean.z", "p.b.mode.z", "p.b.2m.z", "CRR.z",
                        "L2.z", "L3.z", "L4.z", "L.CV.z", "L.skewness.z", "L.kurtosis.z",
                        "median.a.d.z", "mode.a.d.z",
                        "weibull_c.z", "weibull_b.z",

                        # Compute descriptive statistics of rho coordinate
                        "mean.rho", "mean.q.rho", "mean.g.rho", "mean.h.rho", "median.rho", "mode.rho",
                        "max.rho", "min.rho", "var.rho", "sd.rho", "CV.rho", "D.rho", "ID.rho",
                        "kurtosis.rho", "skewness.rho",
                        "p.a.mean.rho", "p.a.mode.rho",
                        "p.b.mean.rho", "p.b.mode.rho", "CRR.rho",
                        "L2.rho", "L3.rho", "L4.rho", "L.CV.rho", "L.skewness.rho", "L.kurtosis.rho",
                        "median.a.d.rho", "mode.a.d.rho",
                        "weibull_c.rho", "weibull_b.rho",

                        # Compute descriptive statistics of r coordinate
                        "mean.r", "mean.q.r", "mean.g.r", "mean.h.r", "median.r", "mode.r",
                        "max.r", "min.r", "var.r", "sd.r", "CV.r", "D.r", "ID.r",
                        "kurtosis.r", "skewness.r",
                        "p.a.mean.r", "p.a.mode.r",
                        "p.b.mean.r", "p.b.mode.r", "CRR.r",
                        "L2.r", "L3.r", "L4.r", "L.CV.r", "L.skewness.r", "L.kurtosis.r",
                        "median.a.d.r", "mode.a.d.r",
                        "weibull_c.r", "weibull_b.r")

        .col.names <-
          .col.names[.col.names %in%
                       var.metr[[.j]][.var.metr[[.j]][.k, var.metr[[.j]]]]]
        if (length(.col.names) > 0) {

          .pts.met <- .points.metrics(rho_seq = .rho.seq, data = .data.tls,
                                      metr = .col.names)
          .stand[[.j]][[.k]] <- cbind(.stand[[.j]][[.k]],
                                      .pts.met[rownames(.stand[[.j]][[.k]]),
                                               .col.names, drop = FALSE])

        }


        # Convert diameters from International System of Units (m) to cm
        .col.names <- sapply(names(.mean.d), paste,
                             c(".tls", ".0.tls", "", ".0"), sep = "")
        .col.names <-
          colnames(.stand[[.j]][[.k]])[colnames(.stand[[.j]][[.k]]) %in%
                                         .col.names]
        .stand[[.j]][[.k]][.col.names] <- .stand[[.j]][[.k]][.col.names] * 100

        # Add stratum
        if (.by.stratum.2)
          .stand[[.j]][[.k]] <- cbind(stratum = .id[.i, "stratum"],
                                      .stand[[.j]][[.k]])

      }

    }


    # Save fixed area/k-tree/angle-count plot results, and write csv files ----
    # containing them if 'save.result' is TRUE

    for (.j in plot.design) {


      # Merge fixed area/k-tree/angle-count plot results for all trees'
      # database(s)

      # Results for first trees' database(s)
      .stand.j <- .stand[[1]][[.j]]

      # Add results for the rest of trees' databases (if any)
      for (.k in names(.stand)[-1]) {

        # Restrict rows to those existing in trees' database to be added
        .stand.j <- .stand.j[rownames(.stand.j) %in% rownames(.stand[[.k]][[.j]]), ,
                             drop = FALSE]

        # Add new columns
        .stand.j <- cbind(.stand.j,
                          .stand[[.k]][[.j]][rownames(.stand.j),
                                             !colnames(.stand[[.k]][[.j]]) %in%
                                               colnames(.stand.j),
                                             drop = FALSE])
      }

      # Save reordered merged results
      rownames(.stand.j) <- NULL
      stand[[.j]] <- rbind(stand[[.j]],
                           .stand.j[, colnames(stand[[.j]]), drop = FALSE])
      stand[[.j]] <- stand[[.j]][order(stand[[.j]][, .plot.design[[.j]]],
                                       stand[[.j]][, "id"]), ,
                                 drop = FALSE]
      rownames(stand[[.j]]) <- NULL

      # Write CSV file
      if (save.result) {

        .file <- paste(switch(funct, sim = "simulations",
                              metr = "metrics.variables",
                              est = "estimation.plot.size"),
                       .j, "plot.csv", sep = ".")
        utils::write.csv(stand[[.j]], file = file.path(dir.result, .file),
                         row.names = FALSE)

      }

    }

    # Define final time, and print message
    t1 <- Sys.time()
    message(" (", format(round(difftime(t1, t0, units = "secs"), 2)), ")")

  }

  return(stand)

}


# Compute Pearson/Spearman correlations

.cor.pearson <- function(x, y){

  stats::cor.test(x = x, y = y, method = 'pearson')

}

.cor.spearman <- function(x, y){

  stats::cor.test(x = x, y = y, method = 'spearman', exact = FALSE)

}


##############################################################################
# Cálculo del índice de curvatura
##############################################################################

# This functions apply calculated ncr index for all points.
# For that purpose, voxelize point cloud by means of regular
# grid in x any z coordinates


.ncr.remove.slice.double <- function(data){

  # code <- NULL

  # Select necessary fields from original txt file of point cloud

  .data <- as.data.frame(data[,c("point", "x", "y", "z")])

  # Create x and y coordinates for grid

  .x <- seq(min(.data$x), max(.data$x)+1)
  .y <- seq(min(.data$y), max(.data$y)+1)


  if(length(.x) < 2 | length(.y) < 2){

    data$ncr <- NA
    .data <- data} else{

  # Empty data frame where coordinates neccesaries for
  # creating grid will be saved

  .grid <- data.frame(x = rep(.x, each = length(.y)),
                      y = rep(.y, times = length(.x)))

  .grid <- sf::st_as_sf(.grid, coords = c("x","y"))


  .grid <- sf::st_buffer(.grid, dist = 0.55, endCapStyle = "SQUARE")
  .grid <- sf::st_cast(.grid, "POLYGON")

  .grid.2 <- sf::st_intersects(.grid, sf::st_as_sf(.data, coords = c("x","y")))
  .grid.2 <- as.data.frame(.grid.2)
  colnames(.grid.2) <- c("id", "code")
  .data$code <- as.numeric(row.names(.data))
  .data <- merge(.data, .grid.2, by = "code", all = TRUE)
  # .data <- subset(.data, select = -code)
  .data <- .data[, 2:ncol(.data)]

  rm(.grid, .grid.2)

  .dat <- lapply(split(.data[, 1:4], .data$id), as.matrix)

  .dat <- .dat[(lapply(.dat, nrow)) > 4]
  # .dat <- .dat[(lapply(.dat, nrow)) < 40000]

  .ncr <- do.call(rbind, lapply(.dat, ncr_point_cloud_double))

  .ncr <- .ncr[.ncr$ncr > 0 & .ncr$ncr < 9999, ]


  if(is.null(.ncr)){
    data$ncr <- NA
    .data <- data
    return(.data)}

  .data <- merge(data, .ncr, by = "point", all = TRUE)

  .data <- .data[!duplicated(.data), ]

  }

  return(.data)

}


.ver.remove.slice.double <- function(data){

  # code <- NULL

  # Select necessary fields from original txt file of point cloud

  .data <- as.data.frame(data[,c("point", "x", "y", "z")])

  # Create x and y coordinates for grid

  .x <- seq(min(.data$x), max(.data$x)+1)
  .y <- seq(min(.data$y), max(.data$y)+1)


  if(length(.x) < 2 | length(.y) < 2){

    data$ncr <- NA
    .data <- data} else{

      # Empty data frame where coordinates neccesaries for
      # creating grid will be saved

      .grid <- data.frame(x = rep(.x, each = length(.y)),
                          y = rep(.y, times = length(.x)))

      .grid <- sf::st_as_sf(.grid, coords = c("x","y"))


      .grid <- sf::st_buffer(.grid, dist = 0.55, endCapStyle = "SQUARE")
      .grid <- sf::st_cast(.grid, "POLYGON")

      .grid.2 <- sf::st_intersects(.grid, sf::st_as_sf(.data, coords = c("x","y")))
      .grid.2 <- as.data.frame(.grid.2)
      colnames(.grid.2) <- c("id", "code")
      .data$code <- as.numeric(row.names(.data))
      .data <- merge(.data, .grid.2, by = "code", all = TRUE)
      # .data <- subset(.data, select = -code)
      .data <- .data[, 2:ncol(.data)]

      rm(.grid, .grid.2)

      .dat <- lapply(split(.data[, 1:4], .data$id), as.matrix)

      .dat <- .dat[(lapply(.dat, nrow)) > 4]
      # .dat <- .dat[(lapply(.dat, nrow)) < 40000]

      .ncr <- do.call(rbind, lapply(.dat, ver_point_cloud_double))

      .ncr <- .ncr[.ncr$ver > 0 & .ncr$ver < 9999, ]


      if(is.null(.ncr)){
        data$ncr <- NA
        .data <- data
        return(.data)}

      .data <- merge(data, .ncr, by = "point", all = TRUE)

      .data <- .data[!duplicated(.data), ]

    }

  return(.data)

}




.no.trees.detected.single <- function(data, d.top, plot.attributes, dir.result, save.result){

  if(is.null(data$id) & is.null(d.top)){

    # If plot identification (id) is not available

    .colnames <- c("tree", "x", "y", "phi", "phi.left", "phi.right", "h.dist", "dbh", "h", "v", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else if (is.null(data$id) & !is.null(d.top)) {

    # If plot identification (id) is not available

    .colnames <- c("tree", "x", "y", "phi", "phi.left", "phi.right", "h.dist", "dbh", "h", "h.com", "v", "v.com", "n.pts", "SS.max", "sinuosity", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else if (!is.null(data$id) & is.null(d.top)) {

    # If plot identification (id) is available

    .colnames <- c("id", "file", "tree", "x", "y", "phi", "phi.left", "phi.right", "h.dist", "dbh", "h", "v", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else {

    # If plot identification (id) is available

    .colnames <- c("id", "file", "tree", "x", "y", "phi", "phi.left", "phi.right", "h.dist", "dbh", "h", "h.com", "v", "v.com", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  }

  .tree <- data.frame(matrix(nrow = 1, ncol = length(.colnames), dimnames = list(NULL, .colnames)))

  if(!is.null(data$id)){

  .tree$id <- data$id[1]
  .tree$file <- data$file[1]

  }

  # Lastly, aggregate attributes table
  if(!is.null(plot.attributes))
    .tree <- merge(.tree, plot.attributes, by = "id", all = FALSE)


  if(isTRUE(save.result)){

    utils::write.csv(.tree,
                     file = file.path(dir.result, "tree.tls.csv"),
                     row.names = FALSE)
  }


  #####
  return(.tree)


}



.no.trees.detected.multi <- function(data, d.top, plot.attributes, dir.result, save.result){

  if(is.null(data$id) & is.null(d.top)){

    # If plot identification (id) is not available

    .colnames <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else if (is.null(data$id) & !is.null(d.top)) {

    # If plot identification (id) is not available

    .colnames <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "h.com", "v", "v.com", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else if (!is.null(data$id) & is.null(d.top)) {

    # If plot identification (id) is available

    .colnames <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else {

    # If plot identification (id) is available

    .colnames <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "h.com", "v", "v.com", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  }

  .tree <- data.frame(matrix(nrow = 1, ncol = length(.colnames), dimnames = list(NULL, .colnames)))


  if(!is.null(data$id)){

    .tree$id <- data$id[1]
    .tree$file <- data$file[1]

  }

  # Lastly, aggregate attributes table
  if(!is.null(plot.attributes))
    .tree <- merge(.tree, plot.attributes, by = "id", all = FALSE)


  if(isTRUE(save.result)){

    utils::write.csv(.tree,
                     file = file.path(dir.result, "tree.tls.csv"),
                     row.names = FALSE)
  }


  #####
  return(.tree)


}








