
tree.detection.multi.scan <- function(data, single.tree = NULL,
                                      dbh.min = 4, dbh.max = 200, h.min = 1.3,
                                      ncr.threshold = 0.1, tls.precision = NULL,
                                      stem.section = NULL, stem.range = NULL, breaks = NULL,
                                      slice = 0.1, understory = NULL, bark.roughness = 1,
                                      den.type = 1, d.top = NULL,
                                      plot.attributes = NULL,
                                      save.result = TRUE, dir.result = NULL){


  set.seed(123)

  #### Checking some function arguments ####

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()

  # Converting arguments to International System of Units
  # Arguments of the forest inventory

  .dbh.min <- dbh.min / 100
  .dbh.max <- dbh.max / 100

  # Obtaining Cartesian coordinates (x,y) from center

  kk <- data[data$phi < (pi/2) & data$prob.selec == 1, c("x", "y", "phi", "rho")]

  x.center <- mean(kk$x - sin(kk$phi) * kk$rho)
  y.center <- mean(kk$y - cos(kk$phi) * kk$rho)

  rm(kk)


  #### Detecting possible areas with trees in the point cloud ####

  if(!is.null(data$GLA)){
    woody <- data[data$GLA <= 0, ]
    woody <- woody[!is.na(woody$x) & !is.na(woody$y) & !is.na(woody$z), ]} else {woody <- data}


  # Statistical filtering of a point cloud
  # Implements the Statistical Outliers Removal (SOR)

  message("Statistical filtering of the whole point cloud")

  woody <- woody[, c("x", "y", "z")]
  woody <- VoxR::filter_noise(data = data.table::setDT(woody), store_noise = TRUE, message = FALSE)

  noise <- woody[woody$Noise == 2, ]
  woody <- woody[woody$Noise == 1, ]

  woody <- merge(data, woody[, c("x", "y", "z")], by = c("x", "y", "z"), all = FALSE)
  noise <- merge(data, noise, by = c("x", "y", "z"), all = FALSE)


  # Detection of stem part without shrub vegetation and crown

  message("Detecting tree stem axes")

  stem <- woody[woody$prob.selec == 1, ]

  if(!is.null(data$intensity)){

    stem <- stem[stem$intensity > mean(stem$intensity, na.rm = T), ]

  }


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

  # Creation of raster with projected voxels

  stem <- VoxR::project_voxels(stem)

  # Filtering pixels - double branch peeling

  stem.2 <- NULL

  stem <- stem[stem$npts > mean(stem$npts) & stem$ratio > mean(stem$ratio) & stem$nvox > mean(stem$nvox), ]

  if(!is.null(understory)){

    if(is.null(single.tree))
      stem.2 <- stem[stem$npts > mean(stem$npts), ]

    stem <- stem[stem$npts > mean(stem$npts) & stem$ratio > mean(stem$ratio) & stem$nvox > mean(stem$nvox), ]

    }


  # Creation polygon to extract those projected areas in the original point cloud
  # where trees are probably located

  buf <- sp::SpatialPoints(cbind(stem$x,stem$y))
  raster::crs(buf) <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +unit=m"
  buf <- suppressWarnings(raster::buffer(buf, width = max(.dbh.min, 0.5), dissolve = TRUE))
  buf <- buf@polygons[[1]]@Polygons
  buf <- lapply(seq_along(buf), function(i) sp::Polygons(list(buf[[i]]), ID = i))
  buf <- sp::SpatialPolygons(buf)
  # raster::plot(buf)


  if(!is.null(understory) & is.null(single.tree)){

  buf.2 <- sp::SpatialPoints(cbind(stem.2$x,stem.2$y))
  raster::crs(buf.2) <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +unit=m"
  buf.2 <- suppressWarnings(raster::buffer(buf.2, width = 0.1, dissolve = TRUE))

  buf.1 <- sf::st_as_sf(buf)
  sf::st_crs(buf.1) <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +unit=m"
  buf.2 <- sf::st_as_sf(buf.2)
  sf::st_crs(buf.2) <- "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +unit=m"


  buf.2 <- sf::st_difference(sf::st_union(buf.2), sf::st_buffer(sf::st_union(buf.1), 0.1))

  rm(buf.1)


  # Detection of stem part without shrub vegetation and crown

  stem.2 <- woody[woody$prob.selec == 1, ]

  buf.2 <- sf::as_Spatial(buf.2)
  buf.2 <- buf.2@polygons[[1]]@Polygons
  retained <- data.frame(tree = 1:length(buf.2),
                         T = unlist(lapply(seq_along(buf.2), function(i) buf.2[[i]]@area>(pi/4)*(.dbh.min+0.1)^2)))
  retained <- retained[retained$T == TRUE, ]

  buf.2 <- lapply(seq_along(buf.2), function(i) sp::Polygons(list(buf.2[[i]]), ID = i))
  buf.2 <- sp::SpatialPolygons(buf.2)
  stem.2$tree <- sp::over(sp::SpatialPoints(coords = cbind(stem.2$x,stem.2$y,stem.2$z)), buf.2, returnlist=TRUE)
  stem.2 <- stem.2[!is.na(stem.2$tree), ]
  stem.2 <- stem.2[stem.2$tree %in% retained$tree, ]

  rm(buf.2, retained)

  }


  # Assigning points to trees previously detected

  stem <- woody[woody$prob.selec == 1, ]

  stem$tree <- sp::over(sp::SpatialPoints(coords = cbind(stem$x,stem$y,stem$z)), buf, returnlist=TRUE)
  stem <- stem[!is.na(stem$tree), ]


  # Filtering stems axis

  stem.i <- do.call(rbind, lapply(split(stem, stem$tree), .n.w.ratio))
  .Q1 <- stats::quantile(stem.i$n.w.ratio, prob = 0.25, na.rm = TRUE)
  .Q3 <- stats::quantile(stem.i$n.w.ratio, prob = 0.75, na.rm = TRUE)
  stem.i <- stem.i[stem.i$n.w.ratio >= .Q1 - 1.5 * (.Q3 - .Q1) & stem.i$n.w.ratio <= .Q3 + 1.5 * (.Q3 - .Q1), ]
  stem <- stem[stem$tree %in% stem.i$tree, ]

  rm(stem.i)


  # If there is only one tree in the point cloud

  if(!is.null(single.tree)){

    filter <- data.frame(table(stem$tree))
    filter <- filter[order(filter$Freq, decreasing = TRUE), ]
    stem <- stem[stem$tree == filter$Var1[1], ]
    # stem$tree <- 1

  }


  # Breaks argument

  if(is.null(breaks)){
    breaks <- seq(from = 0.4, to = max(stem$z), by = 0.3)
    breaks <- breaks[-length(breaks)]}


  # Defining stem axis

  eje <- do.call(rbind, lapply(split(stem, stem$tree), .stem.axis, scan.approach = "multi"))
  .Q1 <- stats::quantile(eje$n.w.ratio, prob = 0.25, na.rm = TRUE)
  .Q3 <- stats::quantile(eje$n.w.ratio, prob = 0.75, na.rm = TRUE)
  eje <- eje[eje$n.w.ratio >= .Q1 - 1.5 * (.Q3 - .Q1), ]
  eje <- eje[eje$tree %in% eje$tree, ]
  eje <- eje[eje$sec %in% as.character(breaks) & !is.na(eje$x), ]

  rm(stem)

  # Removing those axis ver very closed (0.25 m) to each other

  eje.2 <- eje[eje$sec == 1.3, ]
  dbscan <- dbscan::dbscan(eje.2[, c("x", "y"), drop = FALSE], eps = 0.25, minPts = 1)
  eje.2$cluster <- dbscan$cluster
  # eje.2 <- eje.2[!duplicated(eje.2$cluster), ]
  n.w.ratio <- tapply(eje.2$n.w.ratio, eje.2$cluster, max)
  eje <- eje[eje$n.w.ratio %in% n.w.ratio, ]

  rm(eje.2)


  # Estimating NCR threshold when RGB is available

  # if(!is.null(data$GLA)){
  #   ncr <- data[data$GLA > 0, ]
  #   ncr <- as.matrix(ncr[, c("point", "x", "y", "z")])
  #   ncr.threshold <- ncr_point_cloud_double(ncr[1:10000, ])
  #   ncr.threshold <- mean(ncr.threshold$ncr, na.rm = TRUE)}


  # Assigning points to trees previously detected

  woody$tree <- sp::over(sp::SpatialPoints(coords = cbind(woody$x,woody$y,woody$z)), buf, returnlist=TRUE)
  woody <- woody[woody$tree %in% eje$tree, ]
  woody <- woody[!is.na(woody$tree), ]

  rm(buf)


  # If there is only one tree in the point cloud

  # if(!is.null(single.tree)){
  #
  #   filter <- data.frame(table(woody$tree))
  #   filter <- filter[order(filter$Freq, decreasing = TRUE), ]
  #   woody <- woody[woody$tree == filter$Var1[1], ]
  #   woody$tree <- 1
  #
  # }


  #### Starting with clustering process ####

  .filteraux <- data.frame(cluster = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(),
                           n.pts = as.numeric(), n.pts.red = as.numeric(),
                           phi.left = as.numeric(), phi.right = as.numeric(),
                           circ = as.numeric(), arc.circ = as.numeric(), sec = as.numeric())

  .filteraux.2 <- data.frame(cluster = as.numeric(),
                             center.x = as.numeric(), center.y = as.numeric(),
                             center.phi = as.numeric(), center.rho = as.numeric(),
                             center.r = as.numeric(), center.theta = as.numeric(),
                             radius = as.numeric(),
                             n.pts = as.numeric(), n.pts.red = as.numeric(),
                             phi.left = as.numeric(), phi.right = as.numeric(),
                             circ = as.numeric(), arc.circ = as.numeric(), sec = as.numeric())

  # Defining slice

  slice <- slice / 2


  for(cuts in breaks){


    message("Computing section: ", cuts, " m")

    .cut <- woody[which(woody$z > (cuts-slice-0.05) & woody$z < (cuts+slice+0.05)), , drop = FALSE]

    if(nrow(.cut) < 50){next}

    .cut <- .ncr.remove.slice.double(.cut)

    .cut <- .cut[which(.cut$ncr < ncr.threshold | is.na(.cut$ncr)), , drop = FALSE]


    # Restrict to slice corresponding to cuts m +/- 5 cm

    .cut <- .cut[which(.cut$z > (cuts-slice) & .cut$z < (cuts+slice)), , drop = FALSE]

    # Dbscan parameters

    if(is.null(tls.precision)){.eps <- 0.03} else {.eps <- tls.precision}

    # Clustering

    .error <- try(suppressMessages(dbscan::dbscan(.cut[, c("x", "y"), drop = FALSE], eps = .eps, minPts = 15)))
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
                                     .dbh.min = .dbh.min, .dbh.max = .dbh.max,
                                     slice = slice * 2, bark.roughness = bark.roughness,
                                     x.center = x.center, y.center = y.center))

    .filteraux <- rbind(.filteraux, .filter)


    # Second...


    if(!is.null(stem.2)){
    .cut <- stem.2[which(stem.2$z > (cuts-slice-0.05) & stem.2$z < (cuts+slice+0.05)), , drop = FALSE]

    if(nrow(.cut) < 50){next}

    .cut <- .ncr.remove.slice.double(.cut)

    .cut <- .cut[which(.cut$ncr < ncr.threshold | is.na(.cut$ncr)), , drop = FALSE]


    # Restrict to slice corresponding to cuts m +/- 5 cm

    .cut <- .cut[which(.cut$z > (cuts-slice) & .cut$z < (cuts+slice)), , drop = FALSE]

    # Dbscan parameters

    if(is.null(tls.precision)){.eps <- 0.03} else {.eps <- tls.precision}

    # Clustering

    .error <- try(suppressMessages(dbscan::dbscan(.cut[, c("x", "y"), drop = FALSE], eps = .eps, minPts = 15)))
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
                                     .dbh.min = .dbh.min, .dbh.max = .dbh.max,
                                     slice = slice * 2, bark.roughness = bark.roughness,
                                     x.center = x.center, y.center = y.center))

    .filteraux.2 <- rbind(.filteraux.2, .filter)}


  }# End of cuts loop


  #### Assigning sections to tree axis ####

  if(nrow(.filteraux) < 1 & nrow(.filteraux.2) < 1) {

    # Generate a warning and create empty data.frame to be returned, if no row
    # was included in .filter

    warning("No tree was detected")

    .tree <- .no.trees.detected.multi(data, d.top, plot.attributes, dir.result, save.result)
    return(.tree)

  } else {

    if(nrow(.filteraux.2) > 0 & sum(c(.filteraux.2$circ, .filteraux.2$arc.circ)) > 0){

    # Assigning sections to tree axis

    eje.2 <- .stem.assignment.single.scan(.filteraux.2, stem.section, breaks)

    eje.2$tree <- eje.2$tree +  max(eje$tree)

    eje <- rbind(eje[, c("tree", "sec", "x", "y")], eje.2[, c("tree", "sec", "x", "y")])
    rm(eje.2)

    .filter <- rbind(.filteraux, .filteraux.2)
    rm(.filteraux, .filteraux.2)

    }  else {

    .filter <- .filteraux
    rm(.filteraux)

    }


  .filter <- .stem.assignment.multi.scan(.filter, eje, stem.section, x.center, y.center, single.tree)


  if(nrow(.filter) < 1){

    warning("No tree was detected")

    .tree <- .no.trees.detected.multi(data, d.top, plot.attributes, dir.result, save.result)
    return(.tree)

  }


  # Export all section detected

  .filter <- .filter[order(.filter$tree, .filter$sec), , drop = FALSE]
  .stem <- .filter[, c("tree", "sec", "center.x",  "center.y", "radius")]
  .stem$radius <- .stem$radius * 200
  colnames(.stem) <- c("tree", "sec", "x",  "y", "dbh")



  ####################################################################################
  ### Calculate of taper coefficient as the slope coefficient of linear regression ###
  ####################################################################################

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

  top.lim <- max(1.3, abs(max(stem.section) - 1.3))

  .radio.est <- data.frame(radio.est = as.numeric())

  for (i in unique(.filter$tree)) {

    .dat <- .filter[which(.filter$tree == i), ]
    .dat$dif <- abs(.dat$dif)

    if(min(.dat$dif) > top.lim)
      next

    .dat <- .dat[order(.dat$dif), ]

    if(.dat$dif[1] == 0 & .dat$circ[1] == 1){

       .dat <- .dat[1, ]

    } else if (nrow(.dat) > 1 & .dat$dif[1] == 0) {

      .dat <- .dat[.dat$dif == 0 | .dat$dif == .dat$dif[2], ]

    } else {

      .dat <- .dat[.dat$dif == min(.dat$dif), ]

    }


    .out <- data.frame(radio.est = mean(.dat$radio.est, na.rm = TRUE))
    .radio.est <- rbind(.radio.est, .out)

  }



  # Tree normal section coordinates

  .filteraux <- data.frame(tree = as.numeric(),

                           filter = as.numeric(),

                        center.x = as.numeric(),
                        center.y = as.numeric(),
                        center.phi = as.numeric(),
                        center.rho = as.numeric(),
                        center.r = as.numeric(),
                        center.theta = as.numeric(),

                        sec.x = as.numeric(),
                        sec.y = as.numeric(),
                        sec.max = as.numeric(),

                        horizontal.distance = as.numeric(),
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

    .dat$filter <- nrow(.dat)


    if(nrow(.dat) > 3)
      .dat <- .dat[1:3, ]

    if(min(abs(.dat$dif)) > top.lim)
      next

    .filteraux <- rbind(.filteraux, .dat)

  }


  if(nrow(.filteraux) < 1){

    warning("No tree was detected")

    .tree <- .no.trees.detected.multi(data, d.top, plot.attributes, dir.result, save.result)
    return(.tree)

  }

  # Genrating dendrometric variables

  .tree <- data.frame(tree = tapply(.filteraux$tree, .filteraux$tree, mean, na.rm = TRUE),
                      filter = tapply(.filteraux$filter, .filteraux$tree, mean, na.rm = TRUE),

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



  # Cheking if there are negative radious

  .tree <- .tree[.tree$radius > 0, ]

  if(nrow(.tree) < 1){

    warning("No tree was detected")

    .tree <- .no.trees.detected.multi(data, d.top, plot.attributes, dir.result, save.result)
    return(.tree)

  }


  # Cheking minimum and maximum radius defined in the arguments

  .tree <- .tree[.tree$radius >= .dbh.min / 2 & .tree$radius <= .dbh.max / 2, ]

  if(nrow(.tree) < 1){

    warning("No tree was detected")

    .tree <- .no.trees.detected.multi(data, d.top, plot.attributes, dir.result, save.result)
    return(.tree)

  }


  # Selecting only those trees with more than one section detected when more than two breaks have been specified


  if(length(breaks) > 3)
    .tree <- .tree[.tree$filter >  1, ]

  if(nrow(.tree) < 1){

    warning("No tree was detected")

    .tree <- .no.trees.detected.multi(data, d.top, plot.attributes, dir.result, save.result)
    return(.tree)

  }


  # Ordering by distance and numbering trees from 1 to n trees

  .tree <- .tree[!duplicated(.tree$x) & !duplicated(.tree$y), ]
  .tree <- .tree[!duplicated(.tree$sec.x) & !duplicated(.tree$sec.y), ]

  if(nrow(.tree) < 1){

    warning("No tree was detected")

    .tree <- .no.trees.detected.multi(data, d.top, plot.attributes, dir.result, save.result)
    return(.tree)

  }


  .tree <- .tree[order(.tree$rho), ]
  .tree$tree <- 1:nrow(.tree)


  # Detecting possible trees overlaped

  # if(nrow(.tree) < 2 & !is.null(single.tree)){.tree.2 <- .tree}
  if(nrow(.tree) < 2){.tree.2 <- .tree}


  if(nrow(.tree) > 1 & is.null(single.tree)){

  .tree.2 <- data.frame(tree = as.numeric(), filter = as.numeric(),

                        x = as.numeric(), y = as.numeric(),

                        sec.x = as.numeric(), sec.y = as.numeric(), sec.max = as.numeric(),

                        phi = as.numeric(), rho = as.numeric(), r = as.numeric(), theta = as.numeric(),

                        horizontal.distance = as.numeric(), radius = as.numeric(),

                        partial.occlusion = as.numeric(),

                        n.pts = as.numeric(), n.pts.red = as.numeric())


  for (i in unique(.tree$tree)) {


    if(nrow(.tree[.tree$tree == i, ]) < 1)
      next

    .filt <- .tree[.tree$tree == i, ]
    .filteraux <- .tree[.tree$tree != i, ]

    if(nrow(.filteraux) < 1)
      next

    .filteraux$dist <- sqrt((.filteraux$x - .filt$x) ^ 2 + (.filteraux$y - .filt$y) ^ 2) - .filteraux$radius - .filt$radius
    .filteraux$rho.dist <- abs(.filteraux$rho - .filt$rho) - .filteraux$radius - .filt$radius
    .filteraux$phi.dist <- abs(.filteraux$phi - .filt$phi)

    if(min(.filteraux$dist) < 0 |
       .filteraux$rho.dist[.filteraux$dist == min(.filteraux$dist)] < mean(.tree$radius) &
       .filteraux$phi.dist[.filteraux$dist == min(.filteraux$dist)] < 0.1){

      .filteraux <- .filteraux[.filteraux$dist < 0 | .filteraux$rho.dist < mean(.tree$radius) & .filteraux$phi.dist < 0.1, ]

      .filteraux <- rbind(.filt, .filteraux[ , 1:(ncol(.filteraux)-3)])

      .filt <- .filteraux[.filteraux$filter == max(.filteraux$filter) & .filteraux$sec.max ==  min(.filteraux$sec.max), ]


      if(nrow(.filt) > 1)
        .filt <- .filt[.filt$partial.occlusion > 0, ]


      if(nrow(.filt) > 1)
        .filt <- .filt[1, ]

      .tree.remove <- .filteraux[.filteraux$tree != .filt$tree, ]$tree


      if(length(.tree.remove) < 1){.tree <- .tree} else {
        suppressWarnings(.tree <- .tree[.tree$tree != .tree.remove, ])
      }


      .tree.2 <- rbind(.tree.2, .filt)


    } else {

      .tree.2 <- rbind(.tree.2, .filt)
      .tree <- .tree[.tree$tree != .filt$tree, ]

    }

  }

  }

  .tree <- .tree.2
  rm(.tree.2)


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

  data <- data[data$z >= h.min, ]
  s <- sample(nrow(data), round(nrow(data)*0.1))
  data <- data[s, ]
  data <- data[, c("id", "file", "x", "y", "z", "rho")]

  # If only one tree is detected, Voronoi tessellation is not working

  if(nrow(.tree) == 1){

    .P99 <- data.frame(tree = .tree$tree, h = stats::quantile(data$z, prob = 0.9999999999))

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
                   if(length(z) < 1) {
                     P99 <- 0
                     names(P99) <- tree.3[tree.3$id == id, "tree"]
                   } else {
                     P99 <-
                       height_perc_cpp(rho_seq = Inf, z = z, rho = z)[, "P99.9"]
                     names(P99) <- tree.3[tree.3$id == id, "tree"]}
                   return(P99)
                 },
                 voro = .voro, tree.3 = .tree.3)
  .P99 <- data.frame(tree = names(.P99), h = .P99)

  rm(.tree.2, .tree.3, .voro)

  }


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

  stem.v <- .volume(.stem, id = data$id[1])
  .tree <- merge(.tree, stem.v, all = TRUE)

  } else if (length(table(.stem$hi)) > 3 & !is.null(d.top)) {

  stem.v <- .volume(.stem, d.top, id = data$id[1])
  .tree <- merge(.tree, stem.v, all = TRUE)

  } else if (length(table(.stem$hi)) <= 3 & !is.null(d.top)) {

  den.type <- 1
  n <- den.type

  # Estimating volume according to dendrometric type

  .tree$v <- pi * (.tree[, "h"] ^ (n + 1) / (n + 1)) * ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ n)
  h.lim <- (((d.top / 200) ^ 2) / ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ n)) ^ (1 / n)
  .tree$v.com <- pi * ((.tree[, "h"] ^ (n + 1) - h.lim ^ (n + 1)) / (n + 1)) * ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ n)
  .tree$v.com <- ifelse(.tree$v.com < 0, 0, .tree$v.com)
  .tree$h.com <- h.lim

  } else {

  den.type <- 1
  n <- den.type

  .tree$v <- pi * (.tree[, "h"] ^ (n + 1) / (n + 1)) * ((.tree[, "dbh"] / 200) ^ 2 / (.tree[, "h"] - 1.3) ^ n)

  }

  .stem <- .stem[, c("tree", "x", "y", "dhi", "dbh", "hi", "h")]
  .stem <- .stem[order(.stem$tree, .stem$hi), , drop = FALSE]

  if(!is.null(data$id)){
    .stem$id <- data$id[1]
    .stem <- .stem[, c("id", "tree", "x", "y", "dhi", "dbh", "hi", "h")]}



  # Straighness analysis

  straightness <- do.call(rbind, lapply(split(.stem, .stem$tree), .straightness, stem.range = stem.range))

  .tree <- merge(.tree, straightness, by = "tree")

  utils::write.csv(.stem,
                   file = file.path(dir.result, "stem.curve.csv"),
                   row.names = FALSE)

  rm(.stem)

  # If plot identification (id) is not available

  if(is.null(data$id) & is.null(.tree$v.com)){

    .tree <- .tree[, c("tree", "x", "y", "phi", "horizontal.distance", "dbh", "h", "v", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else if (is.null(data$id) & !is.null(.tree$v.com)) {

    .tree <- .tree[, c("tree", "x", "y", "phi", "horizontal.distance", "dbh", "h", "h.com", "v", "v.com", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("tree", "x", "y", "phi", "h.dist", "dbh", "h", "h.com", "v", "v.com", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")


  } else if (!is.null(data$id) & is.null(.tree$v.com)) {

    # If plot identification (id) is available

    .tree$id <- data$id[1]
    .tree$file <- data$file[1]

    .tree <- .tree[, c("id", "file", "tree", "x", "y", "phi", "horizontal.distance", "dbh", "h", "v", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "v", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

  } else {

    # If plot identification (id) is available

    .tree$id <- data$id[1]
    .tree$file <- data$file[1]

    .tree <- .tree[, c("id", "file", "tree", "x", "y", "phi", "horizontal.distance", "dbh", "h", "h.com", "v", "v.com", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
    colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "h.dist", "dbh", "h", "h.com", "v", "v.com", "SS.max", "sinuosity", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

    }

  }

  .tree <- .tree[order(.tree$h.dist), ]
  .tree$tree <- 1:nrow(.tree)

  # Removing values of 0 in n.pts
  .tree$n.pts <- ifelse(.tree$n.pts < 1, 0.01, .tree$n.pts)
  .tree$n.pts.red <- ifelse(.tree$n.pts.red < 1, 0.01, .tree$n.pts.red)
  .tree$n.pts.est <- ifelse(.tree$n.pts.est < 1, 0.01, .tree$n.pts.est)
  .tree$n.pts.red.est <- ifelse(.tree$n.pts.red.est < 1, 0.01, .tree$n.pts.red.est)


  # Lastly, aggregate attributes table
  if(!is.null(plot.attributes))
    .tree <- merge(.tree, plot.attributes, by = "id", all = FALSE)


  if(isTRUE(save.result)){

    utils::write.csv(.tree,
                     file = file.path(dir.result, "tree.tls.csv"),
                     row.names = FALSE)
  }


  if(isTRUE(save.result)){

    .data.red <- noise[which(noise$prob.selec == 1), , drop = FALSE]

    vroom::vroom_write(.data.red, path = file.path(dir.result, paste("noise_", .data.red$file[1], sep = "")), delim = ",", progress = FALSE)

  }


  #####

  return(.tree)

}
