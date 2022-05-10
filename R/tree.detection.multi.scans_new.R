
tree.detection.multi.scans_new <- function(data, dbh.min = 4, dbh.max = 200, h.min = 1.3,
                                           ncr.threshold = 0.1, tls.precision = NULL, breaks = NULL,
                                           plot.attributes = NULL, stem.section = NULL,
                                           save.result = TRUE, dir.result = NULL){

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()

  # Converting arguments to International System of Units
  # Arguments of the forest inventory

  .dbh.min <- dbh.min / 100
  .dbh.max <- dbh.max / 100


  # Detection of stem part without shrub vegetation and crown
  stem <- data[data$prob.selec == 1, ]
  if(!is.null(stem$GLI))
    stem <- stem[stem$GLI <= 0, ]
  stem <- stem[!is.na(stem$x) & !is.na(stem$y) & !is.na(stem$z), ]


  if(is.null(stem.section)){

    den <- .getStem(stem)
    stem <- stem[stem$z > den$x & stem$z < den$x + den$diff, ]
    rm(den)

  } else {

    stem <- stem[stem$z > stem.section[1] & stem$z < stem.section[2], ]

  }


  if(is.null(tls.precision)){
  stem <- VoxR::vox(stem[, c("x", "y", "z")], res = 0.03)} else {
    stem <- VoxR::vox(stem[, c("x", "y", "z")], res = tls.precision)}

  stem <- stem[, c("x", "y", "z", "npts")]

  stem <- VoxR::project_voxels(stem)
  # plot(stem$x, stem$y, col = "grey", asp = 1, pch = 19, cex = 0.5)
  stem <- stem[stem$npts > mean(stem$npts), ]
  # points(stem$x, stem$y, pch = 19, cex = 0.5, col = "black")

  buf <- sp::SpatialPoints(cbind(stem$x,stem$y))
  buf <- raster::buffer(buf, width = 37500, dissolve = TRUE)
  buf <- buf@polygons[[1]]@Polygons
  buf <- lapply(seq_along(buf), function(i) sp::Polygons(list(buf[[i]]), ID = i))
  buf <- sp::SpatialPolygons(buf)
  # sp::plot(buf, col = "red")


  stem <- data[data$prob.selec == 1, ]
  if(!is.null(stem$GLI))
    stem <- stem[stem$GLI <= 0, ]
  stem <- stem[!is.na(stem$x) & !is.na(stem$y) & !is.na(stem$z), ]
  # stem <- stem[stem$z > den$x & stem$z < den$x + den$diff, ]


  stem$tree <- sp::over(sp::SpatialPoints(coords = cbind(stem$x,stem$y,stem$z)), buf, returnlist=TRUE)

  stem <- stem[!is.na(stem$tree), ]
  # plot(stem$x, stem$y, col = stem$tree, asp = 1)

  # Breaks argument
  if(is.null(breaks)){
    breaks <- seq(from = 0.1, to = max(stem$z), by = 0.3)
    breaks <- breaks[-length(breaks)]}


  # Defining stem axis

  stem.i <- split(stem, stem$tree)
  eje <- do.call(rbind, lapply(stem.i, .stem.axis))
  eje <- eje[eje$sec %in% as.character(breaks) & !is.na(eje$x), ]

  rm(stem.i)


  # write.csv(eje, "ejes.csv")


  # .ncr.threshold <- .ncr.threshold.double(data)

  .filteraux <- data.frame(cluster = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(),
                           n.pts = as.numeric(), n.pts.red = as.numeric(),
                           phi.left = as.numeric(), phi.right = as.numeric(),
                           circ = as.numeric(), arc.circ = as.numeric(), sec = as.numeric())



  for(cuts in breaks){


    message("Computing section: ", cuts, " m")

    .cut <- stem[which(stem$z > (cuts-0.1) & stem$z < (cuts+0.1)), , drop = FALSE]

    if(nrow(.cut) < 100)
      next

    .cut <- .ncr.remove.slice.double(.cut)

    .cut <- .cut[which(.cut$ncr < ncr.threshold | is.na(.cut$ncr)), , drop = FALSE]


    # vroom::vroom_write(.cut, path = file.path(dir.result, paste(cuts, ".txt", sep = "")), delim = ",", progress = FALSE)


    # Restrict to slice corresponding to cuts m +/- 5 cm
    .cut <- .cut[which(.cut$z > (cuts-0.05) & .cut$z < (cuts+0.05)), , drop = FALSE]

    # Dbscan parameters
    if(is.null(tls.precision)){.eps <- .dbh.min} else {.eps <- tls.precision}

    # Clustering
    .error <- try(suppressMessages(dbscan::dbscan(.cut[, c("x", "y"), drop = FALSE], eps = .eps)))
    if(class(.error)[1] == "try-error"){
      message("No computed section: ", cuts, " m")
      next} else {
    .dbscan <- dbscan::dbscan(.cut[, c("x", "y"), drop = FALSE], eps = .eps)
    .cut$cluster <- .dbscan$cluster
    .cut <- .cut[which(.cut$cluster > 0), , drop = FALSE]}

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
                          n.pts = as.numeric(), n.pts.red = as.numeric(),

                          # Circumference
                          circ = as.numeric(),

                          # Circumference arc
                          arc.circ = as.numeric(),

                          # Partial occlusion
                          occlusion = as.numeric())


    for(i in unique(.cut$cluster)){

      .pb$tick()

      # Select cluster i
      .dat <- .cut[which(.cut$cluster == i), , drop = FALSE]

      if(nrow(.dat) < 10)
        next

      # plot(.dat$x, .dat$y, asp = 1, main = i)

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
      .threeshold <- median(.density, na.rm = T)

      if(is.nan(.threeshold) | is.na(.threeshold)){next}

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

      if(nrow(.dat) < 1){next}

      # plot(.dat$x, .dat$y, asp = 1, main = i)


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

      .center.phi <- atan2(.center.y, .center.x)
      .center.phi <- ifelse(.center.phi < 0, .center.phi + (2 * pi), .center.phi)
      .center.rho <- sqrt(.center.x ^ 2 + .center.y ^ 2)
      .center.r <- sqrt(.dat$sec[1] ^ 2 + .center.rho ^ 2)
      .center.theta <- atan2(.dat$sec[1], .center.rho)

      # points(.center.x, .center.y, col = "red", pch = 19)
      # abline(h = .center.y)
      # abline(v = .center.x)


      # Distances between points and center
      .dat$dist <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE)

      # Radius value as the mean distance
      # .dat <- .dat[order(.dat$dist, decreasing = FALSE), , drop = FALSE]
      # plot(.dat$x, .dat$y, asp = 1)
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25) & .dat$dist<stats::quantile(.dat$dist, prob = 0.9)])
      if(.radio < 0 | is.na(.radio)){next}
      # xyc4<-conicfit::calculateCircle(.center.x,.center.y,.radio)
      # points(xyc4[,1],xyc4[,2],col='red',type='l')


      # c4 <- conicfit::LMcircleFit(as.matrix(.dat[.dat$dist>stats::quantile(.dat$dist, prob = 0.25), c("x", "y")]))
      # xyc4<-conicfit::calculateCircle(c4[1],c4[2],c4[3])
      # points(xyc4[,1],xyc4[,2],col='cyan',type='l')
      # .radio2 <- c4[3]
      .radio2 <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.5) & .dat$dist<stats::quantile(.dat$dist, prob = 0.9)])

      # Coefficient of variation for distances among cluster points and the estimated center
      .cv <- stats::sd(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25)]) / .radio
      if(.cv > 0.1 | length(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25)]) < 2){next}

      # At least 95 % of distances should be greater than .radio / 2
      if(stats::quantile(.dat$dist, prob = 0.05) < (.radio / 2)) {next}


      .dat.2 <- .dat[order(.dat$x, .dat$y, decreasing = F), ]

      .dat.2$a <- sqrt((.dat.2$x - .dat.2$x[1])^2+(.dat.2$y - .dat.2$y[1])^2)
      .dat.2$b <- sqrt((.dat.2$x[1] - .center.x)^2+(.dat.2$y[1] - .center.y)^2)
      .dat.2$c <- sqrt((.dat.2$x - .center.x)^2+(.dat.2$y - .center.y)^2)

      .dat.2$alpha <- acos(-(.dat.2$a^2-.dat.2$b^2-.dat.2$c^2)/(2*.dat.2$b*.dat.2$c))

      # .dat.2$x2 <- .dat.2$c * cos(.dat.2$alpha)
      # .dat.2$y2 <- .dat.2$c * sin(.dat.2$alpha)


      # Select 1st percentil, if necessary for strange points
      # It remains to be seen what happens if cluster is located in 0 +/- phi
      # .pto.left <- stats::quantile(.dat.2$alpha, prob = 0.01)
      # .x.left <- mean(.dat.2$x2[which(.dat.2$alpha <= .pto.left)])
      # .y.left <- mean(.dat.2$y2[which(.dat.2$alpha <= .pto.left)])

      # Select 99th percentil, if necessary for strange points
      # .pto.right <- stats::quantile(.dat.2$alpha, prob = 0.99)
      # .x.right <- mean(.dat.2$x2[which(.dat.2$alpha >= .pto.right)])
      # .y.right <- mean(.dat.2$y2[which(.dat.2$alpha >= .pto.right)])

      # For points in section center, select those in half the angle aperture
      # phi +/- TLS aperture .alpha
      # .phi.cent <- max(.dat.2$alpha) - ((max(.dat.2$alpha) - min(.dat.2$alpha)) / 2)
      # .x.cent <- mean(.dat.2$x2[which(round(.dat.2$alpha, 1) >= round(.phi.cent, 1) & round(.dat.2$alpha, 1) <= round(.phi.cent, 1))])
      # .y.cent <- mean(.dat.2$y2[which(round(.dat.2$alpha, 1) >= round(.phi.cent, 1) & round(.dat.2$alpha, 1) <= round(.phi.cent, 1))])
      #
      # plot(.dat.2$x2, .dat.2$y2, asp = 1, main = .cluster)
      # points(.x.cent, .y.cent, col = "green", pch = 19)
      # points(.x.left, .y.left, col = "red", pch = 19)
      # points(.x.right, .y.right, col = "blue", pch = 19)

      # if(is.nan(.rho.cent)){next}

      # Checking if tree section belong to a full circumference

      k1 <- ifelse(nrow(.dat[.dat$x > .center.x & .dat$y > .center.y, ]) > 0, 1, 0)
      k2 <- ifelse(nrow(.dat[.dat$x > .center.x & .dat$y < .center.y, ]) > 0, 1, 0)
      k3 <- ifelse(nrow(.dat[.dat$x < .center.x & .dat$y > .center.y, ]) > 0, 1, 0)
      k4 <- ifelse(nrow(.dat[.dat$x < .center.x & .dat$y < .center.y, ]) > 0, 1, 0)

      .circ <- ifelse(sum(k1, k2, k3, k4) == 4, 1, 0)
      .arc.circ <- ifelse(sum(k1, k2, k3, k4) == 3 | sum(k1, k2, k3, k4) == 2, 1, 0)

      # Convert original coordinates if cluster is located in 0 +/- phi
      # .phi.left <- ifelse(.phi.left > (2 * pi), .phi.left - (2 * pi), .phi.left)
      # .phi.right <- ifelse(.phi.right > (2 * pi), .phi.right - (2 * pi), .phi.right)

      # If the complete circumference arc is not detected (so previous criterion
      # is not satisfied), check if there is a partial occlusion.
      # If they belong to a tree section, they are found with a systematic
      # regularity so very high correlations between their correlative
      # numeration and phi must be found when they are ordered with respect to phi

      .dat.2 <- .dat.2[order(.dat.2$alpha, decreasing = F), ]
      # .cv <- stats::sd(diff(.dat.2$alpha)) / mean(diff(.dat.2$alpha))
      # if(.cv > 0.1 | length(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25)]) < 2){next}
      .dat.2$n <- c(1:nrow(.dat.2))
      .cor <- try(stats::cor.test(x = .dat.2$n, y = .dat.2$alpha, method = 'pearson'), silent = TRUE) # cor function could be used instead

      # plot(.dat.2$n, .dat.2$alpha)
      # .cor <- try(stats::cor.test(x = .dat.2$n, y = .dat.2$phi, method = 'spearman'))

      # If error, go to next iteration
      if(class(.cor) == "try-error"){next} else{

        .occlusion <- .cor[[4]]
        .occlusion.sig <- .cor["p.value"]

      }

      # plot(.dat$x, .dat$y, asp = 1, main = paste(.cluster, .circ, .arc.circ, sep = " "))


      # Zhang et al., (2019)
      .n.w.ratio <- stats::sd(.dat$z) / sqrt(stats::sd(.dat$x) ^ 2 + stats::sd(.dat$y) ^ 2)

      if(.n.w.ratio > 1 | is.nan(.n.w.ratio)){next}

      .densidad_radio <- .n.pts.red / .radio


      if(nrow(.dat) < 10){next}


      # Results
      .salida <- data.frame(cluster = .dat$cluster[1],

                            center.x = .center.x, center.y = .center.y,
                            center.phi = .center.phi, center.rho = .center.rho,
                            center.r = .center.r, center.theta = .center.theta,

                            radius = .radio, radius2 = .radio2,

                            n.pts = .n.pts, n.pts.red = .n.pts.red,

                            circ = .circ, arc.circ = .arc.circ, occlusion = .occlusion,

                            density.radio = .densidad_radio)

      if(!exists(".salida")){next} else {.filter <- rbind(.filter, .salida)}


    }

    # Arch of circumference or partial arch of circumference?
    .Q1 <- stats::quantile(.filter$density.radio, prob = 0.25)
    .Q3 <- stats::quantile(.filter$density.radio, prob = 0.75)
    .outliers <- .Q1 - 1.5 * (.Q3 - .Q1)

    # Minimum number of points
    # .filter$tree <- ifelse(.filter$n.pts.red > .outliers, 1, 0)
    # .filter <- .filter[which(.filter$tree == 1), , drop = FALSE]

    .filter$tree <- ifelse(.filter$circ == 1 & .filter$density.radio >= .outliers, 1,
                           ifelse(.filter$arc.circ == 1 & .filter$occlusion > 0.95 & .filter$density.radio >= .outliers, 1,
                                  ifelse(.filter$circ == 0 & .filter$arc.circ == 0 & .filter$occlusion > 0.975 & .filter$density.radio >= .outliers, 1, 0)))
    .filter <- .filter[which(.filter$tree == 1), , drop = FALSE]

    # Dbh maximum and minimum
    .filter$tree <- ifelse(.filter$radius > (.dbh.max / 2) | .filter$radius < (.dbh.min / 2), 0, 1)
    .filter <- subset(.filter, .filter$tree == 1)


    if(nrow(.filter) < 1){

      .filter1.0 <- data.frame(cluster = as.numeric(),
                               center.x = as.numeric(), center.y = as.numeric(),
                               center.phi = as.numeric(), center.rho = as.numeric(),
                               center.r = as.numeric(), center.theta = as.numeric(),
                               radius = as.numeric(), radius2 = as.numeric(),
                               n.pts = as.numeric(), n.pts.red = as.numeric(),
                               circ = as.numeric(), arc.cir = as.numeric(), sec = as.numeric())

    } else{

      .filter1.0 <- .filter[, c("cluster",
                                "center.x", "center.y", "center.phi", "center.rho", "center.r", "center.theta",
                                "radius", "radius2", "n.pts", "n.pts.red", "circ", "arc.circ"), drop = FALSE]
      .filter1.0$sec <- cuts

    }

    .filteraux <- rbind(.filteraux,.filter1.0)

  }# End of cuts loop

  .filter <- .filteraux
  rm(.filteraux)

  if(nrow(.filter) < 1) stop("No tree was detected")


  # Assigning sections to tree axis
  eje$sec <- as.character(eje$sec)
  .filter$sec <- as.character(.filter$sec)
  .filter <- merge(eje, .filter, by = "sec", all.y = TRUE)
  .filter$dist <- sqrt((.filter$center.x - .filter$x) ^ 2 + (.filter$center.y - .filter$y) ^ 2)
  .filter$sec <- as.numeric(.filter$sec)

  .filteraux <- data.frame(tree = as.numeric(), sec = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(), radius2 = as.numeric(),
                           n.pts = as.numeric(), n.pts.red = as.numeric(),
                           circ = as.numeric(), arc.cir = as.numeric())

  for (i in unique(.filter$tree)) {
    for (j in unique(.filter$sec)) {

      .filt <- .filter[.filter$tree == i & .filter$sec == j, ]
      .filt <- .filt[.filt$dist == min(.filt$dist), ]
      .filt <- .filt[, c("tree", "sec", "dist",
                         "center.x", "center.y", "center.phi", "center.rho", "center.r", "center.theta",
                         "radius", "radius2", "n.pts", "n.pts.red", "circ", "arc.circ")]

      .filteraux <- rbind(.filteraux, .filt)

    }

  }

  # .filteraux <- .filteraux[.filteraux$dist < 2 * .filteraux$radius, ]

  .filteraux <- .filteraux[order(.filteraux$tree, .filteraux$sec), , drop = FALSE]


  .filter <- data.frame(tree = as.numeric(), sec = as.numeric(), dist = as.numeric(),
                        center.x = as.numeric(), center.y = as.numeric(),
                        center.phi = as.numeric(), center.rho = as.numeric(),
                        center.r = as.numeric(), center.theta = as.numeric(),
                        radius = as.numeric(), radius2 = as.numeric(),
                        n.pts = as.numeric(), n.pts.red = as.numeric(),
                        circ = as.numeric(), arc.cir = as.numeric())

  for (i in unique(.filteraux$tree)) {

    .filt <- .filteraux[.filteraux$tree == i, ]

    .filt <- .filt[.filt$dist < 2 * mean(.filt$radius, na.rm = TRUE), ]

    if(nrow(.filt) < 1)
      next

    .filt$dif <- c(0, abs(diff(.filt$radius)))

    while (max(.filt$dif) > 0.1) {

      .filt <- .filt[.filt$dif <= 0.1, ]

    }

    .filt <- .filt[, -ncol(.filt)]
    .filter <- rbind(.filter, .filt)

  }

  if(nrow(.filter) < 1) stop("No tree was detected")

  # .filter <- .filteraux
  rm(.filteraux)


  # Export all section detected
  .filter <- .filter[order(.filter$tree, .filter$sec), , drop = FALSE]
  .stem <- .filter[, c("tree", "sec", "center.x",  "center.y", "radius", "radius2")]
  .stem$radius <- .stem$radius * 200
  .stem$radius2 <- .stem$radius2 * 200
  colnames(.stem) <- c("tree", "sec", "x",  "y", "dbh", "dbh2")


  # Calculate of taper coefficient as the slope coefficient of linear regression

  # .taper <- .filter[.filter$sec >= 1 & .filter$sec <= 1.6, c("tree", "sec", "radius"), drop = FALSE]
  .taper <- .filter[, c("tree", "sec", "radius", "radius2"), drop = FALSE]

  .slope.tree <- data.frame(tree = as.numeric(), slope = as.numeric(), slope2 = as.numeric())

  for(i in unique(.taper$tree)){

  .taper.i <- .taper[.taper$tree == i, ]

  if(nrow(.taper.i) < 2){

    .slope <- data.frame(tree = i, slope = NA, slope2 = NA)

  } else {

    .lm <- stats::lm(radius ~ sec, data = .taper.i)
    .lm2 <- stats::lm(radius2 ~ sec, data = .taper.i)
    .slope <- data.frame(tree = i, slope = stats::coef(.lm)[2], slope2 = stats::coef(.lm2)[2])

  }

  .slope.tree <- rbind(.slope.tree, .slope)

  }

  .filter <- merge(.filter, .slope.tree, by = "tree")



  # If a new column ("dif") is created containing the difference between dbh and
  # section from which radius is estimated, number of used sections (or cuts)
  # will not be important. An estimated radius will always exist
  .filter$dif <- 1.3 - .filter$sec
  .filter$radio.est <- ifelse(.filter$dif == 0, .filter$radius, .filter$radius + .filter$slope * .filter$dif)
  .filter$radio.est2 <- ifelse(.filter$dif == 0, .filter$radius2, .filter$radius2 + .filter$slope2 * .filter$dif)
  .filter <- .filter[!is.na(.filter$radio.est), ]

  # When there are not enough tree to the linear model, and tress were
  # detected at at different sections than 1.3 m, we assume these radius
  # as dbh. This could happen in very few situations.
  # .filter$radio.est <- ifelse(is.na(.filter$radio.est), .filter$radius, .filter$radio.est)
  # .filter$radio.est2 <- ifelse(is.na(.filter$radio.est2), .filter$radius2, .filter$radio.est2)

  .radio.est <- data.frame(radio.est = as.numeric(), radio.est2 = as.numeric())

  for (i in unique(.filter$tree)) {

    .dat <- .filter[which(.filter$tree == i), ]
    .dat$dif <- abs(.dat$dif)

    if(min(.dat$dif) > 1.3)
      next

    .dat <- .dat[order(.dat$dif), ]

    if(.dat$dif[1] == 0 & .dat$circ[1] == 1){

       .dat <- .dat[1, ]

    } else if (nrow(.dat) > 1 & .dat$dif[1] == 0) {

      .dat <- .dat[.dat$dif == 0 | .dat$dif == .dat$dif[2], ]

    } else{

      .dat <- .dat[.dat$dif == min(.dat$dif), ]

    }




    # if(max(.dat$circ) == 1){
    #
    #   .dat <- .dat[which(.dat$circ > 0), ]
    #
    # } else if(max(.dat$circ) == 0 & max(.dat$arc.circ) == 1){
    #
    #   .dat <- .dat[which(.dat$arc.circ > 0), ]
    #
    # } else {.dat <- .dat}

    .out <- data.frame(radio.est = mean(.dat$radio.est, na.rm = TRUE), radio.est2 = mean(.dat$radio.est2, na.rm = TRUE))
    .radio.est <- rbind(.radio.est, .out)

  }

  # Tree normal section coordinates
  .filteraux <- data.frame(tree = as.numeric(),

                        center.x = as.numeric(),
                        center.y = as.numeric(),
                        center.phi = as.numeric(),
                        center.rho = as.numeric(),
                        center.r = as.numeric(),
                        center.theta = as.numeric(),

                        horizontal.distance = as.numeric(), # repeated line
                        radius = as.numeric(),
                        radius2 = as.numeric(),

                        partial.occlusion = as.numeric(),

                        n.pts = as.numeric(),
                        n.pts.red = as.numeric())



  for (i in unique(.filter$tree)) {

    .dat <- .filter[which(.filter$tree == i), ]
    .dat <- .dat[order(abs(.dat$dif)), ]

    if(nrow(.dat) > 3)
      .dat <- .dat[1:3, ]

    if(min(abs(.dat$dif)) > 1.3)
      next

    .filteraux <- rbind(.filteraux, .dat)

  }

  # Dendrometric variables
  .tree <- data.frame(tree = tapply(.filteraux$tree, .filteraux$tree, mean, na.rm = TRUE),

                      x = tapply(.filteraux$center.x, .filteraux$tree, mean, na.rm = TRUE),
                      y = tapply(.filteraux$center.y, .filteraux$tree, mean, na.rm = TRUE),
                      phi = tapply(.filteraux$center.phi, .filteraux$tree, mean, na.rm = TRUE),
                      rho = tapply(.filteraux$center.rho, .filteraux$tree, mean, na.rm = TRUE),
                      r = tapply(.filteraux$center.r, .filteraux$tree, mean, na.rm = TRUE),
                      theta = tapply(.filteraux$center.theta, .filteraux$tree, mean,na.rm = TRUE),

                      horizontal.distance = tapply(.filteraux$center.rho, .filteraux$tree, mean, na.rm = TRUE), # repeated line
                      radius = .radio.est$radio.est,
                      radius2 = .radio.est$radio.est2,

                      partial.occlusion = tapply(.filteraux$arc.circ, .filteraux$tree, min, na.rm = TRUE),

                      n.pts = tapply(.filteraux$n.pts, .filteraux$tree, mean, na.rm = TRUE),
                      n.pts.red = tapply(.filteraux$n.pts.red, .filteraux$tree, mean, na.rm = TRUE))

  rm(.filteraux)

  # Numbering trees from 1 to n trees
  .tree$tree <- 1:nrow(.tree)

  # Indicate trees with partial occlusions, those for which none of the sections
  # was identified as circumference arch (ArcCirc)
  .tree$partial.occlusion <- ifelse(.tree$partial.occlusion == 0, 1, 0)

  # Compute dbh (cm)
  .tree$dbh <- .tree$radius * 200
  .tree$dbh2 <- .tree$radius2 * 200


  # Calculate points belonging to radius unit
  # Since it will be an estimation, select sections completely visible
  # (ArcCirc == 1) in section corresponding to 1.3 m (where dbh is estimated)
  .filter$filter <- ifelse(.filter$sec == 1.3 & .filter$circ == 1, 1, 0)
  .filter2 <- subset(.filter, .filter$filter == 1)

  if(nrow(.filter2) < 1)
    .filter2 <- .filter

  # Estimate number of points by cluster, with and without point cropping
  # process, corresponding to radius 1 m
  .filter2$points.radio <- .filter2$n.pts / .filter2$radius
  .filter2$points.radio.hom <- .filter2$n.pts.red / .filter2$radius

  # Average points after point cropping by m of radius
  .tree$points.m <- mean(.filter2$points.radio)
  .tree$points.m.hom <- mean(.filter2$points.radio.hom)

  # Finally, compute number of points estimated for each tree according to
  # radius
  .tree$n.pts.est <- .tree$points.m * .tree$radius
  .tree$n.pts.red.est <- .tree$points.m.hom * .tree$radius

  # Obtaining reduced point cloud
  data <- data[data$prob.selec == 1, ]
  data <- data[, c("id", "file", "x", "y", "z", "rho")]

  # If only one tree is detected, Voronoi tessellation is not working

  if(nrow(.tree) == 1){

    .P99 <- data.frame(tree = .tree$tree, P99.9 = stats::quantile(data$z, prob = 0.999))

  } else {

  # Voronoi tessellation
  # .voro <- .tree[ , c("tree", "x", "y"), drop = FALSE]
  # .voro <- ggvoronoi::voronoi_polygon(.voro, x = "x", y = "y", outline = NULL,
  #                                     data.frame = FALSE)
  #
  # sp::coordinates(data) <- ~ x + y + z
  # .voro <- sp::over(data, .voro)
  # .voro <- sp::SpatialPointsDataFrame(data, .voro)
  #
  # .voro <- as.data.frame(.voro, stringsAsFactors = FALSE)
  # .voro <- .voro[!is.na(.voro$tree), , drop = FALSE]


  # setwd("F:/TLS/Huelva/Latizal_con_gestion")
  # .tree <- read.csv("tree.tls.csv")
  .tree.2 <- .tree[ , c("tree", "x", "y"), drop = FALSE]
  .tree.2$tree <- 1:nrow(.tree.2)
  # colnames(.tree) <- c("tree", "center.x", "center.y")

  # data <- read.table("1.txt", header = TRUE, sep = ",")
  # .voro <- data[data$prob > 0.75, ]
  .voro <- data[, c("x", "y", "z")]
  .voro <- sf::st_as_sf(.voro, coords = c("x", "y"))

  .tree.2 <- sf::st_as_sf(.tree, coords = c("x", "y"))
  .voronoi <- sf::st_collection_extract(sf::st_voronoi(do.call(c, sf::st_geometry(.tree.2))))
  id <- sf::st_intersects(.voronoi, .tree.2)

  .voro$tree <- unlist(sf::st_intersects(.voro, .voronoi))
  freq <- table(.voro$tree)


  # data$pol <- pols[unlist(sf::st_intersects(data, pols))]
  # data$tree <- sf::st_intersects(data, pols)


  .tree.2 <- data.frame(tree = as.numeric())

  for (i in 1:nrow(freq)) {

    tre <- data.frame(tree = rep(id[[i]], each = freq[i]))
    .tree.2 <- rbind(.tree.2, tre)

  }


  .voro$tree <- .tree.2$tree


  rm(freq, .tree.2)

  # Compute height percentile P99.9
  .P99 <- sapply(sort(unique(.voro$tree)),
                 function(tree, voro) {
                   z <- voro$z[voro$tree == tree]
                   P99 <-
                     height_perc_cpp(rho_seq = Inf, z = z, rho = z)[, "P99.9"]
                   names(P99) <- tree
                   return(P99)
                 },
                 voro =.voro)
  .P99 <- data.frame(tree = names(.P99), h = .P99)

  }

  rm(.voro)

  # Remove possible trees above "h.min" (1.3 m by default)

  .P99 <- .P99[.P99$h >= h.min, ]
  .tree <- merge(.tree, .P99, by = "tree", all = FALSE)

  rm(.P99)


  # Assigning dbh to stem dataset

  .stem <- merge(.stem, .tree[, c("tree", "h")], by = "tree")
  .stem <- .stem[.stem$sec != 1.3, ]
  .tree$sec <- 1.3
  .stem <- rbind(.stem, .tree[, c("tree", "sec", "x", "y", "dbh", "dbh2", "h")])
  .tree <- .tree[, -ncol(.tree)]
  .stem <- .stem[order(.stem$tree, .stem$sec), , drop = FALSE]
  colnames(.stem) <- c("tree", "sec", "x", "y", "dhi", "dhi2", "h")
  .stem <- merge(.stem, .tree[, c("tree", "dbh")], by = "tree", all.x = TRUE)
  colnames(.stem) <- c("tree", "hi", "x", "y", "dhi", "dhi2", "h", "dbh")


  # Compute volume (m3)

  # .stem <- read.csv("stem.csv")

  if(length(table(.stem$hi)) > 3){

  # Stem curve

  stem.v <- .volume(.stem)
  .tree <- merge(.tree, stem.v, all = TRUE)

  } else {

  # Paraboloid volume

  .tree$v <- pi * (.tree[, "h"] ^ 2 / 2) * ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ 2)

  }

  .stem <- .stem[, c("tree", "x", "y", "dhi", "dbh", "hi", "h")]
  .stem <- .stem[order(.stem$tree, .stem$hi), , drop = FALSE]

  utils::write.csv(.stem,
                   file = file.path(dir.result, "tree.tls.stem.csv"),
                   row.names = FALSE)


  # If plot identification (id) is not available
  if(is.null(data$id)){

    .tree <- .tree[, c("tree", "x", "y", "phi", "horizontal.distance", "dbh", "dbh2", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("tree", "x", "y", "phi", "h.dist", "dbh", "dbh2", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else{

    # If plot identification (id) is available

    .tree$id <- data$id[1]
    .tree$file <- data$file[1]

    .tree <- .tree[, c("id", "file", "tree", "x", "y", "phi", "horizontal.distance", "dbh", "dbh2", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "dbh2", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  }

  # .tree$id <- as.integer(.tree$id)



  rm(.stem)

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
