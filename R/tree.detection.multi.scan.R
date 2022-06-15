
tree.detection.multi.scan <- function(data,
                                      dbh.min = 4, dbh.max = 200, h.min = 1.3,
                                      ncr.threshold = 0.1, tls.precision = NULL,
                                      stem.section = NULL, breaks = NULL, slice = 0.1,
                                      den.type = 1, d.top = NULL,
                                      plot.attributes = NULL,
                                      save.result = TRUE, dir.result = NULL){


  #### Checking some function arguments ####

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()

  # Converting arguments to International System of Units
  # Arguments of the forest inventory

  .dbh.min <- dbh.min / 100
  .dbh.max <- dbh.max / 100


  #### Detecting possible areas with trees in the point cloud ####

  # Detection of stem part without shrub vegetation and crown

  stem <- data[data$prob.selec == 1, ]

  if(!is.null(data$GLA)){
    stem <- stem[stem$GLA <= 0, ]
    stem <- stem[!is.na(stem$x) & !is.na(stem$y) & !is.na(stem$z), ]} else {stem <- stem}


  # Defining the vertical section in which trees are detected

  if(is.null(stem.section)){

    stem.section <- .getStem(stem)
    stem.section <- c(stem.section$x, stem.section$x + stem.section$diff)
    stem <- stem[stem$z > stem.section[1] & stem$z < stem.section[2], ]

  } else {

    stem <- stem[stem$z > stem.section[1] & stem$z < stem.section[2], ]

  }


  if(is.null(tls.precision)){
  stem <- VoxR::vox(stem[, c("x", "y", "z")], res = 0.03)} else {
    stem <- VoxR::vox(stem[, c("x", "y", "z")], res = tls.precision)}

  stem <- stem[, c("x", "y", "z", "npts")]

  # dat <- data[data$prob > 0.75, ]
  # dat <- VoxR::vox(dat[, c("x", "y", "z")], res = 0.03)
  # dat <- dat[, c("x", "y", "z", "npts")]
  # dat <- VoxR::project_voxels(dat)
  # par(mar = c(4, 4, 1, 1))
  # plot(dat$x, dat$y, col = "grey", asp = 1, pch = 19, cex = 0.5,
  #      xlab = "X (m)", ylab = "Y (m)")

  stem <- VoxR::project_voxels(stem)
  # points(stem$x, stem$y, col = "green", pch = 19, cex = 0.5)
  stem <- stem[stem$npts > mean(stem$npts), ]
  # points(stem$x, stem$y, pch = 19, cex = 0.5, col = "black")
  # legend(x = "topright", legend = c("Total", "Stem section", "Above mean"),
  #        col = c("grey", "green", "black"), pch = 19, bty = "n")

  buf <- sp::SpatialPoints(cbind(stem$x,stem$y))
  buf <- suppressWarnings(raster::buffer(buf, width = 37500, dissolve = TRUE))
  buf <- buf@polygons[[1]]@Polygons
  buf <- lapply(seq_along(buf), function(i) sp::Polygons(list(buf[[i]]), ID = i))
  buf <- sp::SpatialPolygons(buf)
  # sp::plot(buf, col = "red")

  # Detection of stem part without shrub vegetation and crown

  stem <- data[data$prob.selec == 1, ]

  if(!is.null(data$GLA)){
    stem <- stem[stem$GLA <= 0, ]
    stem <- stem[!is.na(stem$x) & !is.na(stem$y) & !is.na(stem$z), ]} else {stem <- stem}



  # Assigning points to trees previously detected

  stem$tree <- sp::over(sp::SpatialPoints(coords = cbind(stem$x,stem$y,stem$z)), buf, returnlist=TRUE)
  stem <- stem[!is.na(stem$tree), ]
  # plot(stem$x, stem$y, col = stem$tree, asp = 1)

  # Breaks argument

  if(is.null(breaks)){
    breaks <- seq(from = 0.4, to = max(stem$z), by = 0.3)
    breaks <- breaks[-length(breaks)]}


  # Defining stem axis

  eje <- do.call(rbind, lapply(split(stem, stem$tree), .stem.axis, scan.approach = "multi"))
  eje <- eje[eje$sec %in% as.character(breaks) & !is.na(eje$x), ]

  rm(stem)


  # Estimating NCR threshold when RGB is available

  # if(!is.null(data$GLA)){
  #   ncr <- data[data$GLA > 0, ]
  #   ncr <- as.matrix(ncr[, c("point", "x", "y", "z")])
  #   ncr.threshold <- ncr_point_cloud_double(ncr[1:10000, ])
  #   ncr.threshold <- mean(ncr.threshold$ncr, na.rm = TRUE)}


  # Remove green parts

  if(!is.null(data$GLA)){
    data <- data[data$GLA <= 0, ]
    data <- data[!is.na(data$x) & !is.na(data$y) & !is.na(data$z), ]}


  #### Starting with clustering process ####

  .filteraux <- data.frame(cluster = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(),
                           n.pts = as.numeric(), n.pts.red = as.numeric(),
                           phi.left = as.numeric(), phi.right = as.numeric(),
                           circ = as.numeric(), arc.circ = as.numeric(), sec = as.numeric())

  slice <- slice / 2


  for(cuts in breaks){


    message("Computing section: ", cuts, " m")

    .cut <- data[which(data$z > (cuts-slice-0.05) & data$z < (cuts+slice+0.05)), , drop = FALSE]

    if(nrow(.cut) < 50){next}

    .cut <- .ncr.remove.slice.double(.cut)

    .cut <- .cut[which(.cut$ncr < ncr.threshold | is.na(.cut$ncr)), , drop = FALSE]


    # Restrict to slice corresponding to cuts m +/- 5 cm

    .cut <- .cut[which(.cut$z > (cuts-slice) & .cut$z < (cuts+slice)), , drop = FALSE]

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

    if(nrow(.cut) < 1){next}

    # Assigning section to the slice

    .cut$sec <- cuts

    # Selection of those cluster belonging to trees

    .filter <- do.call(rbind, lapply(split(.cut, .cut$cluster), .sections.multi.scan,
                                     tls.precision = tls.precision,
                                     .dbh.min = .dbh.min, .dbh.max = .dbh.max))

    .filteraux <- rbind(.filteraux, .filter)

  }# End of cuts loop


  .filter <- .filteraux
  rm(.filteraux)


  #### Assigning sections to tree axis ####

  if(nrow(.filter) < 1) {

    # Generate a warning and create empty data.frame to be returned, if no row
    # was included in .filter

    if(is.null(data$id) & is.null(d.top)){

      # If plot identification (id) is not available

      warning("No tree was detected")

      .colnames <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

    } else if (is.null(data$id) & !is.null(d.top)) {

      # If plot identification (id) is not available

      warning("No tree was detected")

      .colnames <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

    } else if (!is.null(data$id) & is.null(d.top)) {

      # If plot identification (id) is available

      warning("No tree was detected for plot ", data$id[1])

      .colnames <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

    } else {

      # If plot identification (id) is available

      warning("No tree was detected for plot ", data$id[1])

      .colnames <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

    }

    .tree <- data.frame(matrix(nrow = 0, ncol = length(.colnames), dimnames = list(NULL, .colnames)))

  } else {

  eje$sec <- as.character(eje$sec)

  eje$phi <- atan2(eje$y, eje$x)
  eje$phi <- ifelse(eje$phi < 0, eje$phi + (2 * pi), eje$phi)
  eje$rho <- sqrt(eje$x ^ 2 + eje$y ^ 2)

  .filter$cluster <- 1:nrow(.filter)
  .filter$sec <- as.character(.filter$sec)
  .filter <- merge(eje, .filter, by = "sec", all.y = TRUE)
  .filter$dist <- sqrt((.filter$center.x - .filter$x) ^ 2 + (.filter$center.y - .filter$y) ^ 2)
  .filter$dist.rho <- abs(.filter$center.rho - .filter$rho)
  .filter$dist.phi <- abs(.filter$center.phi - .filter$phi)
  .filter$sec <- as.numeric(.filter$sec)

  .filter <- .filter[!is.na(.filter$dist) | !is.na(.filter$dist.rho) | !is.na(.filter$dist.phi), ]

  .filteraux <- data.frame(tree = as.numeric(), sec = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(),
                           n.pts = as.numeric(), n.pts.red = as.numeric(),
                           circ = as.numeric(), arc.cir = as.numeric())

  for (i in unique(.filter$tree)) {
    # for (j in unique(.filter$sec)) {

    # .filt <- .filter[.filter$tree == i & .filter$sec == j, ]
    .filt.tree <- .filter[.filter$tree == i, ]

    for (j in unique(.filt.tree$sec)) {

      .filt <- .filt.tree[.filt.tree$sec == j, ]

      .filt$dist.total <- .filt$dist / max(.filt$dist) +
        .filt$dist.rho / max(.filt$dist.rho) +
        .filt$dist.phi / max(.filt$dist.phi)

      .filt$dist.total <- ifelse(is.nan(.filt$dist.total), 0, .filt$dist.total)

      .filt <- .filt[.filt$dist.total == min(.filt$dist.total), ]

      if(nrow(.filt) > 1)
        .filt <- .filt[.filt$dist.rho == min(.filt$dist.rho), ]

      if(nrow(.filt) > 1)
        .filt <- .filt[.filt$dist.phi == min(.filt$dist.phi), ]

      if(nrow(.filt) > 1)
        .filt <- .filt[1, ]

      if(.filt$dist > 2.5 | .filt$dist.rho > 2.5 | .filt$dist.phi > 0.25){next}

      .filt <- .filt[, c("tree", "sec", "dist",
                         "center.x", "center.y", "center.phi",
                         "center.rho", "center.r", "center.theta",
                         "radius", "n.pts", "n.pts.red", "circ", "arc.circ")]

      .filteraux <- rbind(.filteraux, .filt)

      # .filter <- .filter[.filter$code != .filt$code, ]

    }

  }


  # .filteraux <- .filteraux[.filteraux$dist < 2 * .filteraux$radius, ]

  .filteraux <- .filteraux[order(.filteraux$tree, .filteraux$sec), , drop = FALSE]


  .filter <- data.frame(tree = as.numeric(), sec = as.numeric(), dist = as.numeric(),
                        center.x = as.numeric(), center.y = as.numeric(),
                        center.phi = as.numeric(), center.rho = as.numeric(),
                        center.r = as.numeric(), center.theta = as.numeric(),
                        radius = as.numeric(),
                        n.pts = as.numeric(), n.pts.red = as.numeric(),
                        circ = as.numeric(), arc.cir = as.numeric())

  for (i in unique(.filteraux$tree)) {

    .filt <- .filteraux[.filteraux$tree == i, ]

    .filt <- .filt[.filt$dist < 2 * mean(.filt$radius, na.rm = TRUE), ]

    if(nrow(.filt) < 1){next}

    .filt$dif <- c(0, abs(diff(.filt$radius)))
    .filt$dif.sec <- c(0, abs(diff(.filt$sec)))
    .filt$dif <- ifelse(.filt$dif.sec > 1, .filt$dif / .filt$dif.sec, .filt$dif)

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
  .stem <- .filter[, c("tree", "sec", "center.x",  "center.y", "radius")]
  .stem$radius <- .stem$radius * 200
  colnames(.stem) <- c("tree", "sec", "x",  "y", "dbh")


  # Calculate of taper coefficient as the slope coefficient of linear regression

  # .taper <- .filter[.filter$sec >= 1 & .filter$sec <= 1.6, c("tree", "sec", "radius"), drop = FALSE]
  .taper <- .filter[, c("tree", "sec", "radius"), drop = FALSE]

  .slope.tree <- data.frame(tree = as.numeric(), slope = as.numeric(), slope2 = as.numeric())

  for(i in unique(.taper$tree)){

  .taper.i <- .taper[.taper$tree == i, ]

  if(nrow(.taper.i) < 2){

    .slope <- data.frame(tree = i, slope = NA)

  } else {

    .lm <- stats::lm(radius ~ sec, data = .taper.i)
    .slope <- data.frame(tree = i, slope = stats::coef(.lm)[2])

  }

  .slope.tree <- rbind(.slope.tree, .slope)

  }

  .filter <- merge(.filter, .slope.tree, by = "tree")



  # If a new column ("dif") is created containing the difference between dbh and
  # section from which radius is estimated, number of used sections (or cuts)
  # will not be important. An estimated radius will always exist
  .filter$dif <- 1.3 - .filter$sec
  .filter$radio.est <- ifelse(.filter$dif == 0, .filter$radius, .filter$radius + .filter$slope * .filter$dif)
  .filter <- .filter[!is.na(.filter$radio.est), ]

  # When there are not enough tree to the linear model, and tress were
  # detected at at different sections than 1.3 m, we assume these radius
  # as dbh. This could happen in very few situations.
  # .filter$radio.est <- ifelse(is.na(.filter$radio.est), .filter$radius, .filter$radio.est)
  # .filter$radio.est2 <- ifelse(is.na(.filter$radio.est2), .filter$radius2, .filter$radio.est2)

  .radio.est <- data.frame(radio.est = as.numeric())

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


    .out <- data.frame(radio.est = mean(.dat$radio.est, na.rm = TRUE))
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

                        sec.x = as.numeric(),
                        sec.y = as.numeric(),
                        sec.max = as.numeric(),

                        horizontal.distance = as.numeric(), # repeated line
                        radius = as.numeric(),

                        partial.occlusion = as.numeric(),

                        n.pts = as.numeric(),
                        n.pts.red = as.numeric())



  for (i in unique(.filter$tree)) {

    .dat <- .filter[which(.filter$tree == i), ]
    .dat <- .dat[order(abs(.dat$dif)), ]

    .sec.x <- .dat$center.x[nrow(.dat)]
    .sec.y <- .dat$center.y[nrow(.dat)]
    .sec.max <- .dat$sec[nrow(.dat)]

    .dat$sec.x <- .sec.x
    .dat$sec.y <- .sec.y
    .dat$sec.max <- .sec.max


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

                      sec.x = tapply(.filteraux$sec.x, .filteraux$tree, mean, na.rm = TRUE),
                      sec.y = tapply(.filteraux$sec.y, .filteraux$tree, mean, na.rm = TRUE),
                      sec.max = tapply(.filteraux$sec.max, .filteraux$tree, mean, na.rm = TRUE),

                      phi = tapply(.filteraux$center.phi, .filteraux$tree, mean, na.rm = TRUE),
                      rho = tapply(.filteraux$center.rho, .filteraux$tree, mean, na.rm = TRUE),
                      r = tapply(.filteraux$center.r, .filteraux$tree, mean, na.rm = TRUE),
                      theta = tapply(.filteraux$center.theta, .filteraux$tree, mean,na.rm = TRUE),

                      horizontal.distance = tapply(.filteraux$center.rho, .filteraux$tree, mean, na.rm = TRUE), # repeated line
                      radius = .radio.est$radio.est,

                      partial.occlusion = tapply(.filteraux$arc.circ, .filteraux$tree, min, na.rm = TRUE),

                      n.pts = tapply(.filteraux$n.pts, .filteraux$tree, mean, na.rm = TRUE),
                      n.pts.red = tapply(.filteraux$n.pts.red, .filteraux$tree, mean, na.rm = TRUE))

  rm(.filteraux)

  # Ordering by distance and numbering trees from 1 to n trees
  .tree <- .tree[!duplicated(.tree$x) & !duplicated(.tree$y), ]
  .tree <- .tree[.tree$radius > 0, ]
  .tree <- .tree[order(.tree$rho), ]
  .tree$tree <- 1:nrow(.tree)

  # Indicate trees with partial occlusions, those for which none of the sections
  # was identified as circumference arch (ArcCirc)
  .tree$partial.occlusion <- ifelse(.tree$partial.occlusion == 0, 1, 0)

  # Compute dbh (cm)
  .tree$dbh <- .tree$radius * 200


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


  #### Estimating tree heights ####

  # Obtaining reduced point cloud
  data <- data[data$z >= h.min & data$prob > 0.9, ]
  data <- data[, c("id", "file", "x", "y", "z", "rho")]

  # If only one tree is detected, Voronoi tessellation is not working

  if(nrow(.tree) == 1){

    .P99 <- data.frame(tree = .tree$tree, P99.9 = stats::quantile(data$z, prob = 0.999))

  } else {

  # Voronoi tessellation

  .tree.2 <- .tree[ , c("tree", "sec.x", "sec.y"), drop = FALSE]
  .tree.2 <- .tree.2[!duplicated(.tree.2$sec.x) & !duplicated(.tree.2$sec.y), ]
  colnames(.tree.2) <- c("tree", "x", "y")
  .tree.2$tree <- 1:nrow(.tree.2)

  .sec <- .tree[ , c("tree", "sec.max", "radius"), drop = FALSE]
  .sec$tree <- 1:nrow(.sec)

  .voro <- data[, c("x", "y", "z")]
  .voro <- sf::st_as_sf(.voro, coords = c("x", "y"))

  .tree.2 <- sf::st_as_sf(.tree.2, coords = c("x", "y"))
  .voronoi <- sf::st_buffer(.tree.2, dist = .sec$radius * 3)
  .voro <- suppressWarnings(sf::st_intersection(.voro, .voronoi))

  .voronoi <- sf::st_collection_extract(sf::st_voronoi(do.call(c, sf::st_geometry(.tree.2))))
  # .voronoi <- sf::st_collection_extract(sf::st_buffer(dist = .sec$radius * 3, do.call(c, sf::st_geometry(.tree.2))))

  .tree.3 <- sf::st_intersects(.voronoi, .tree.2)
  .tree.3 <- data.frame(id = 1:length(.tree.3),
                        tree = unlist(.tree.3, recursive = TRUE, use.names = TRUE))
  .tree.3 <- merge(.tree.3, .sec, by = "tree")
  .voro$tree <- unlist(sf::st_intersects(.voro, .voronoi))


  # Compute height percentile P99.9
  .P99 <- sapply(sort(unique(.tree.3$id)),
                 function(id, voro, tree.3) {
                   sec.max <- .tree.3[.tree.3$id == id, "sec.max"]
                   z <- voro$z[voro$tree == id & voro$z > sec.max]
                   P99 <-
                     height_perc_cpp(rho_seq = Inf, z = z, rho = z)[, "P99.9"]
                   names(P99) <- tree.3[tree.3$id == id, "tree"]
                   return(P99)
                 },
                 voro = .voro, tree.3 = .tree.3)
  .P99 <- data.frame(tree = names(.P99), h = .P99)

  }

  rm(.tree.2, .tree.3, .voro)

  # Remove possible trees above "h.min" (1.3 m by default)

  .P99 <- .P99[.P99$h >= h.min, ]
  .tree <- merge(.tree, .P99, by = "tree", all = FALSE)

  rm(.P99)


  #### Estimating stem volume ####

  # Assigning dbh to stem dataset

  .stem <- merge(.stem, .tree[, c("tree", "h")], by = "tree")
  .stem <- .stem[.stem$sec != 1.3, ]
  .tree$sec <- 1.3
  .stem <- rbind(.stem, .tree[, c("tree", "sec", "x", "y", "dbh", "h")])
  .tree <- .tree[, -ncol(.tree)]
  .stem <- .stem[order(.stem$tree, .stem$sec), , drop = FALSE]
  colnames(.stem) <- c("tree", "sec", "x", "y", "dhi", "h")
  .stem <- merge(.stem, .tree[, c("tree", "dbh")], by = "tree", all.x = TRUE)
  colnames(.stem) <- c("tree", "hi", "x", "y", "dhi", "h", "dbh")

  # Cheking the trees h
  .stem <- merge(.stem[, c("tree", "hi", "x", "y", "dhi", "dbh")],
                 data.frame(tree = unique(.stem$tree),
                            h = do.call(rbind, lapply(split(.stem[, c("hi", "h")], .stem$tree), max))),
                 all.x = TRUE, by = "tree")

  # Compute volume (m3)

  if(length(table(.stem$hi)) > 3 & is.null(d.top)){

  # Stem curve

  stem.v <- .volume(.stem)
  .tree <- merge(.tree, stem.v, all = TRUE)

  } else if (length(table(.stem$hi)) > 3 & !is.null(d.top)) {

    stem.v <- .volume(.stem, d.top)
    .tree <- merge(.tree, stem.v, all = TRUE)


  } else if (length(table(.stem$hi)) <= 3 & !is.null(d.top)) {

  den.type <- 1
  n <- den.type
  # Paraboloid volume

  .tree$v <- pi * (.tree[, "h"] ^ (n + 1) / (n + 1)) * ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ n)
  h.lim <- (((d.top / 200) ^ 2) / ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ n)) ^ (1 / n)
  .tree$v.com <- pi * ((.tree[, "h"] ^ (n + 1) - h.lim ^ (n + 1)) / (n + 1)) * ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ n)
  .tree$v.com <- ifelse(.tree$v.com < 0, 0, .tree$v.com)

  } else {

    den.type <- 1
    n <- den.type

  .tree$v <- pi * (.tree[, "h"] ^ (n + 1) / (n + 1)) * ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ n)

  }

  .stem <- .stem[, c("tree", "x", "y", "dhi", "dbh", "hi", "h")]
  .stem <- .stem[order(.stem$tree, .stem$hi), , drop = FALSE]

  utils::write.csv(.stem,
                   file = file.path(dir.result, "tree.tls.stem.csv"),
                   row.names = FALSE)


  # If plot identification (id) is not available
  if(is.null(data$id) & is.null(.tree$v.com)){

    .tree <- .tree[, c("tree", "x", "y", "phi", "horizontal.distance", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else if (is.null(data$id) & !is.null(.tree$v.com)) {

    .tree <- .tree[, c("tree", "x", "y", "phi", "horizontal.distance", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")


  } else if (!is.null(data$id) & is.null(.tree$v.com)) {

    # If plot identification (id) is available

    .tree$id <- data$id[1]
    .tree$file <- data$file[1]

    .tree <- .tree[, c("id", "file", "tree", "x", "y", "phi", "horizontal.distance", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else {

    # If plot identification (id) is available

    .tree$id <- data$id[1]
    .tree$file <- data$file[1]

    .tree <- .tree[, c("id", "file", "tree", "x", "y", "phi", "horizontal.distance", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  }

  rm(.stem)

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
