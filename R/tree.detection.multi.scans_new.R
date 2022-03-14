
tree.detection.multi.scans_new <- function(data, dbh.min = 7.5, dbh.max = 200, h.min = 1.3,
                                           ncr.threshold = 0.1, tls.precision = NULL, breaks=c(1.0, 1.3, 1.6), plot.attributes = NULL,
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
  # k <- hist(stem$z, breaks = "FD")
  # hist(stem$z, ylim = c(-1000000,1000000))
  # hist(stem$z, freq = F)
  # lines(density(stem$z),lwd=2)
  # lines(density(stem$z,bw=bw.ucv(stem$z)),col=2,lwd=2)
  # lines(density(stem$z,bw=bw.SJ(stem$z)),col=4,lwd=2)
  # k <- density(stem$z)
  #
  # den <- data.frame(x = k$x, y = k$y)
  # den <- den[den$y > quantile(den$y, probs = 0.5), ]
  #
  # den$dev1 <- c(diff(den$y), 0)
  # den$dev2 <- c(diff(den$dev1), 0)
  #
  # plot(den$x, den$y, ylim = c(-0.3,0.3))
  # lines(den$x, den$dev1 * 10, col = 2)
  # lines(den$x, den$dev2 * 10, col = 3)

  # den$dev1 <- abs(den$dev1)
  # den <- den[den$dev1 < quantile(den$dev1, probs = 0.5), ]
  # points(den$x, den$dev1, col = "blue")



  den <- .getStem(stem)

  stem <- stem[stem$z > den$x & stem$z < den$x + den$diff, ]

  rm(den)

  if(is.null(tls.precision)){
  stem <- VoxR::vox(stem[, c("x", "y", "z")], res = 0.03)} else {
    stem <- VoxR::vox(stem[, c("x", "y", "z")], res = tls.precision)}

  stem <- VoxR::project_voxels(stem)
  # plot(stem$x, stem$y, col = "grey", asp = 1, pch = 19, cex = 0.5)
  stem <- stem[stem$npts > mean(stem$npts), ]
  # points(stem$x, stem$y, pch = 19, cex = 0.5, col = "black")

  buf <- sp::SpatialPoints(cbind(stem$x,stem$y))
  # buf <- sp::SpatialPointsDataFrame(coords = cbind(stem$x,stem$y), data = stem[, c("x", "y", "z")])

  buf <- raster::buffer(buf, width = 50000, dissolve = TRUE)
  # sp::plot(buf, col = "red")

  stem <- data
  if(!is.null(stem$GLI))
    stem <- stem[stem$GLI <= 0, ]
  stem <- stem[!is.na(stem$x) & !is.na(stem$y) & !is.na(stem$z), ]
  stem <- sp::SpatialPointsDataFrame(coords = cbind(stem$x,stem$y), data = stem)
  stem <- raster::intersect(stem, buf)

  stem <- stem@data

  rm(buf)


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
      .radio <- mean(.dat$dist[.dat$dist>stats::quantile(.dat$dist, prob = 0.25)])

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

                            radius = .radio,

                            n.pts = .n.pts, n.pts.red = .n.pts.red,

                            circ = .circ, arc.circ = .arc.circ, occlusion = .occlusion,

                            density.radio = .densidad_radio)

      .filter <- rbind(.filter, .salida)

    }

    # Arch of circumference or partial arch of circumference?
    .Q1 <- stats::quantile(.filter$density.radio, prob = 0.25)
    .Q3 <- stats::quantile(.filter$density.radio, prob = 0.75)
    .outliers <- .Q1 - 1.5 * (.Q3 - .Q1)

    # Minimum number of points
    # .filter$tree <- ifelse(.filter$n.pts.red > .outliers, 1, 0)
    # .filter <- .filter[which(.filter$tree == 1), , drop = FALSE]

    .filter$tree <- ifelse(.filter$circ == 1 & .filter$density.radio > .outliers, 1,
                           ifelse(.filter$arc.circ == 1 & .filter$occlusion > 0.95 & .filter$density.radio > .outliers, 1,
                                  ifelse(.filter$circ == 0 & .filter$arc.circ == 0 & .filter$occlusion > 0.975 & .filter$density.radio > .outliers, 1, 0)))
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
                               n.pts = as.numeric(), n.pts.red = as.numeric(),
                               circ = as.numeric(), arc.cir = as.numeric(), sec = as.numeric())

    } else{

      .filter1.0 <- .filter[, c("cluster",
                                "center.x", "center.y", "center.phi", "center.rho", "center.r", "center.theta",
                                "radius", "n.pts", "n.pts.red", "circ", "arc.circ"), drop = FALSE]
      .filter1.0$sec <- cuts

    }

    .filteraux<-rbind(.filteraux,.filter1.0)

  }# End of cuts loop

  .filter<-.filteraux
  rm(.filteraux)


  # Repeat this for the remaining sections!!!
  # Merge of all sections ----

  #.filter <- rbind(.filter1.0, .filter1.3, .filter1.6)

  if(nrow(.filter) < 1) stop("No tree was detected")

  .dbscan <- dbscan::dbscan(.filter[, c("center.x", "center.y"), drop = FALSE], eps = stats::quantile(.filter$radius, prob = 0.75), minPts = 1)
  .filter$cluster <- .dbscan$cluster
  .filter <- .filter[order(.filter$cluster, .filter$sec), , drop = FALSE]

  # Export all section detected
  .stem <- .filter[, c("cluster", "sec", "center.x",  "center.y", "center.phi", "center.rho", "center.r", "center.theta", "radius")]
  .stem$radius <- .stem$radius * 200
  colnames(.stem) <- c("tree", "sec", "x",  "y", "phi", "rho", "r", "theta", "dbh")

  utils::write.csv(.stem,
                   file = file.path(dir.result, "tree.tls.stem.csv"),
                   row.names = FALSE)

  rm(.stem)

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

      .dat <- .dat[which(.dat$circ | .dat$arc.circ > 0), ]

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

                      partial.occlusion = tapply(.filter$arc.circ, .filter$cluster, min, na.rm = TRUE),

                      n.pts = tapply(.filter$n.pts, .filter$cluster, mean, na.rm = TRUE),
                      n.pts.red = tapply(.filter$n.pts.red, .filter$cluster, mean, na.rm = TRUE))

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
  .filter2$points.radio <- .filter2$n.pts / .filter2$radio
  .filter2$points.radio.hom <- .filter2$n.pts.red / .filter2$radio

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
  .voro <- .tree[ , c("tree", "center.x", "center.y"), drop = FALSE]
  .voro <- ggvoronoi::voronoi_polygon(.voro, x = "center.x", y = "center.y", outline = NULL,
                                      data.frame = FALSE)

  sp::coordinates(data) <- ~ x + y + z
  .voro <- sp::over(data, .voro)
  .voro <- sp::SpatialPointsDataFrame(data, .voro)

  .voro <- as.data.frame(.voro, stringsAsFactors = FALSE)
  .voro <- .voro[!is.na(.voro$tree), , drop = FALSE]


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
  .P99 <- data.frame(tree = names(.P99), P99.9 = .P99)

  }

  # Remove possible trees above "h.min" (1.3 m by default)

  .P99 <- .P99[.P99$P99.9 >= h.min, ]
  .tree <- merge(.tree, .P99, by = "tree", all = FALSE)

  # Numbering trees again

  .tree$tree <- 1:nrow(.tree)

  # Compute volume (m3)

  # Paraboloid
  .tree$v <- pi * (.tree[, "P99.9"] ^ 2 / 2) * ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "P99.9"] - 1.3) ^ 2)



  # If plot identification (id) is not available
  if(is.null(data$id)){

    .tree <- .tree[, c("tree", "center.x", "center.y", "center.phi", "horizontal.distance", "dbh", "P99.9", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else{

    # If plot identification (id) is available

    .tree$id <- data$id[1]
    .tree$file <- data$file[1]

    .tree <- .tree[, c("id", "file", "tree", "center.x", "center.y", "center.phi", "horizontal.distance", "dbh", "P99.9", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  }

  # .tree$id <- as.integer(.tree$id)

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
