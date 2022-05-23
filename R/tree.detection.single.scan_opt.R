
tree.detection.single.scan_opt <- function(data, dbh.min = 4, dbh.max = 200, h.min = 1.3,
                                       ncr.threshold = 0.1, tls.resolution = NULL,
                                       stem.section = NULL, breaks = NULL,
                                       plot.attributes = NULL,
                                       d.top = NULL,
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

  if(is.null(tls.resolution))
    stop("The argument tls.resolution must be defined")

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

  # Detection of stem part without shrub vegetation and crown

  stem <- data[data$prob.selec == 1, ]

  if(!is.null(data$GLA)){
    stem <- data[data$GLA <= 0, ]
    stem <- stem[!is.na(stem$x) & !is.na(stem$y) & !is.na(stem$z), ]} else {

      stem <- stem

    }


  # Defining the vertical section in which trees are detected

  if(is.null(stem.section)){

    den <- .getStem(stem)
    stem.section <- c(den$x, den$x + den$diff)
    stem <- stem[stem$z > stem.section[1] & stem$z < stem.section[2], ]
    rm(den)

  } else {

    stem <- stem[stem$z > stem.section[1] & stem$z < stem.section[2], ]

  }

  # # Selecting those regions coontaining trees
  #
  # res <- (sin(.alpha.h / 2) * (max(data$rho) / cos(mean(data$slope, na.rm = TRUE))) * 2)
  #
  # stem <- VoxR::vox(stem[, c("x", "y", "z")], res = res * 2)
  # stem <- stem[, c("x", "y", "z", "npts")]
  # plot(stem$x, stem$y, asp = 1, col = "grey")
  # stem <- VoxR::project_voxels(stem)
  # stem <- stem[stem$npts > mean(stem$npts), ]
  # points(stem$x, stem$y, asp = 1, col = "black")
  #
  # buf <- sp::SpatialPoints(cbind(stem$x,stem$y))
  # buf <- raster::buffer(buf, width = 37500, dissolve = TRUE)
  # buf <- buf@polygons[[1]]@Polygons
  # buf <- lapply(seq_along(buf), function(i) sp::Polygons(list(buf[[i]]), ID = i))
  # buf <- sp::SpatialPolygons(buf)
  #
  #
  # # Detection of stem part without shrub vegetation and crown
  #
  stem <- data[data$prob.selec == 1, ]

  if(!is.null(data$GLA)){
    stem <- data[data$GLA <= 0, ]
    stem <- stem[!is.na(stem$x) & !is.na(stem$y) & !is.na(stem$z), ]} else {

      stem <- stem

    }


  # # Assigning points to trees previously detected
  #
  # stem$tree <- sp::over(sp::SpatialPoints(coords = cbind(stem$x,stem$y,stem$z)), buf, returnlist=TRUE)
  # stem <- stem[!is.na(stem$tree), ]


  # Breaks argument

  if(is.null(breaks)){
    breaks <- seq(from = 0.1, to = max(stem$z), by = 0.3)
    breaks <- breaks[-length(breaks)]}


  # # Defining stem axis
  #
  # stem.i <- split(stem, stem$tree)
  # eje <- do.call(rbind, lapply(stem.i, .stem.axis))
  # eje <- eje[eje$sec %in% as.character(breaks) & !is.na(eje$x), ]
  #
  # rm(stem.i)

  # Estimating NCR threshold when RGB are available

  if(!is.null(stem$GLA)){
    ncr <- data[data$GLA > 0, ]
    ncr <- as.matrix(ncr[, c("point", "x", "y", "z")])

    ncr.threshold <- ncr_point_cloud_double(ncr[1:10000, ])
    ncr.threshold <- mean(ncr.threshold$ncr, na.rm = TRUE)}


  # Remove green parts

  if(!is.null(data$GLA)){
    data <- data[data$GLA <= 0, ]
    data <- data[!is.na(data$x) & !is.na(data$y) & !is.na(data$z), ]}


  # Starting with clustering process

  .filteraux <- data.frame(cluster = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(),
                           n.pts = as.numeric(), n.pts.red = as.numeric(),
                           phi.left = as.numeric(), phi.right = as.numeric(),
                           arc.circ = as.numeric(), sec = as.numeric())


  for(cuts in breaks){

    message("Computing section: ", cuts, " m")

    .cut <- data[which(data$z > (cuts-0.1) & data$z < (cuts+0.1)), , drop = FALSE]

    .cut <- .ncr.remove.slice.double(.cut)

    .cut <- .cut[which(.cut$ncr < ncr.threshold | is.na(.cut$ncr)), , drop = FALSE]

    # Restrict to slice corresponding to cuts m +/- 5 cm
    .cut <- .cut[which(.cut$z > (cuts-0.05) & .cut$z < (cuts+0.05)), , drop = FALSE]

    # DBSCAN parameters
    .eps <- (tan(.alpha.h / 2) * (max(.cut$r) / cos(mean(.cut$slope, na.rm = TRUE))) * 2)

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

    .filter <- do.call(rbind, lapply(split(.cut, .cut$cluster), .sections.single.scan,
                                     .alpha.v = .alpha.v, .alpha.h = .alpha.h,
                                     .dbh.min = .dbh.min, .dbh.max = .dbh.max))

    .filteraux<-rbind(.filteraux, .filter)

  }# End of cuts loop


  .filter<-.filteraux
  rm(.filteraux)

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
                     "h.dist", "dbh", "h", "v", "n.pts",
                     "n.pts.red", "n.pts.est", "n.pts.red.est",
                     "partial.occlusion")
    } else {

      # If plot identification (id) is available

      warning("No tree was detected for plot ", data$id[1])

      .colnames <- c("id", "file", "tree", "x", "y", "phi", "phi.left",
                     "phi.right", "h.dist", "dbh", "h", "v", "n.pts",
                     "n.pts.red", "n.pts.est", "n.pts.red.est",
                     "partial.occlusion")
    }

    .tree <- data.frame(matrix(nrow = 0, ncol = length(.colnames),
                               dimnames = list(NULL, .colnames)))
  }

  else {


    # Assigning sections to tree axis
    kk <- data.frame(sec = table(.filter$sec))
    kk <- kk[kk[,2] == max(kk[,2]), ]
    eje <- .filter[.filter$sec < as.numeric(as.character(kk[, 1])) + 0.6 &
                     .filter$sec > as.numeric(as.character(kk[, 1])) - 0.6,
                   c("center.x", "center.y", "center.rho", "center.phi", "radius", "sec")]
    .dbscan <- dbscan::dbscan(eje[, c("center.x", "center.y"), drop = FALSE], eps = mean(eje$radius), minPts = 1)
    eje$tree <- .dbscan$cluster
    plot(eje$center.x, eje$center.y, asp = 1, col = eje$tree)
    eje <- eje[, c("tree", "sec", "center.x", "center.y", "center.rho", "center.phi")]
    colnames(eje) <- c("tree", "sec", "x", "y", "rho", "phi")
    eje <- eje[order(eje$tree), ]

    eje.2 <- data.frame(tree = as.numeric(), sec = as.numeric(),
                        x = as.numeric(), y = as.numeric(),
                        rho = as.numeric(), phi = as.numeric())

    for (i in unique(eje$tree)) {

      kk <- eje[eje$tree == i, ]
      kk <- merge(kk, data.frame(sec = breaks), all.y = TRUE)

      kk$tree <- i
      kk$x <- ifelse(is.na(kk$x), mean(kk$x, na.rm = TRUE), kk$x)
      kk$y <- ifelse(is.na(kk$y), mean(kk$y, na.rm = TRUE), kk$y)
      kk$rho <- ifelse(is.na(kk$rho), mean(kk$rho, na.rm = TRUE), kk$rho)
      kk$phi <- ifelse(is.na(kk$phi), mean(kk$phi, na.rm = TRUE), kk$phi)

      eje.2 <- rbind(eje.2, kk)

    }

    eje <- eje.2
    rm(eje.2)



    eje$sec <- as.character(eje$sec)
    .filter$sec <- as.character(.filter$sec)
    .filter <- merge(eje, .filter, by = "sec", all.y = TRUE)
    .filter$dist <- sqrt((.filter$center.x - .filter$x) ^ 2 + (.filter$center.y - .filter$y) ^ 2)
    .filter$dist.rho <- abs(.filter$center.rho - .filter$rho)
    .filter$dist.phi <- abs(.filter$center.phi - .filter$phi)
    .filter$sec <- as.numeric(.filter$sec)
    # .filter$code <- 1:nrow(.filter)

    .filteraux <- data.frame(tree = as.numeric(), sec = as.numeric(),
                             center.x = as.numeric(), center.y = as.numeric(),
                             center.phi = as.numeric(), phi.left = as.numeric(), phi.right = as.numeric(),
                             center.rho = as.numeric(), center.r = as.numeric(), center.theta = as.numeric(),
                             radius = as.numeric(),
                             n.pts = as.numeric(), n.pts.red = as.numeric(),
                             arc.cir = as.numeric())


    for (i in unique(.filter$tree)) {
      # for (j in unique(.filter$sec)) {

        # .filt <- .filter[.filter$tree == i & .filter$sec == j, ]
        .filt.tree <- .filter[.filter$tree == i, ]

        for (j in unique(.filt.tree$sec)) {

        .filt <- .filt.tree[.filt.tree$sec == j, ]

        .filt$dist.total <- .filt$dist / max(.filt$dist) +
          .filt$dist.rho / max(.filt$dist.rho) +
          .filt$dist.phi / max(.filt$dist.phi)
        .filt <- .filt[.filt$dist.total == min(.filt$dist.total), ]

        if(nrow(.filt) > 1)
          .filt <- .filt[.filt$dist.rho == min(.filt$dist.rho), ]

        if(nrow(.filt) > 1)
          .filt <- .filt[.filt$dist.phi == min(.filt$dist.phi), ]

        if(nrow(.filt) > 1)
          .filt <- .filt[1, ]

        if(.filt$dist > 2.5){next}
        if(.filt$dist.rho > 2.5){next}
        if(.filt$dist.phi > 0.25){next}



        .filt <- .filt[, c("tree", "sec", "dist",
                           "center.x", "center.y", "center.phi", "phi.left", "phi.right",
                           "center.rho", "center.r", "center.theta",
                           "radius", "n.pts", "n.pts.red", "arc.circ")]

        .filteraux <- rbind(.filteraux, .filt)

        # .filter <- .filter[.filter$code != .filt$code, ]

      }

    }

    # .filteraux <- .filteraux[.filteraux$dist < 2 * .filteraux$radius, ]

    .filteraux <- .filteraux[order(.filteraux$tree, .filteraux$sec), , drop = FALSE]


    .filter <- data.frame(tree = as.numeric(), sec = as.numeric(), dist = as.numeric(),
                          center.x = as.numeric(), center.y = as.numeric(),
                          center.phi = as.numeric(), phi.left = as.numeric(), phi.right = as.numeric(),
                          center.rho = as.numeric(),
                          center.r = as.numeric(), center.theta = as.numeric(),
                          radius = as.numeric(),
                          n.pts = as.numeric(), n.pts.red = as.numeric(),
                          arc.cir = as.numeric())

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
    .stem <- .filter[, c("tree", "sec", "center.x",  "center.y", "radius")]
    .stem$radius <- .stem$radius * 200
    colnames(.stem) <- c("tree", "sec", "x",  "y", "dbh")


    # Calculate of taper coefficient as the slope coefficient of linear regression

    # .taper <- .filter[.filter$sec >= 1 & .filter$sec <= 1.6, c("tree", "sec", "radius"), drop = FALSE]
    .taper <- .filter[, c("tree", "sec", "radius"), drop = FALSE]

    .slope.tree <- data.frame(tree = as.numeric(), slope = as.numeric())

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

      if(.dat$dif[1] == 0 & .dat$arc.circ[1] == 1){

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

      .out <- data.frame(radio.est = mean(.dat$radio.est, na.rm = TRUE))
      .radio.est <- rbind(.radio.est, .out)

    }

    # Tree normal section coordinates
    .filteraux <- data.frame(tree = as.numeric(),

                             center.x = as.numeric(),
                             center.y = as.numeric(),
                             center.phi = as.numeric(),
                             phi.left = as.numeric(),
                             phi.right = as.numeric(),
                             center.rho = as.numeric(),
                             center.r = as.numeric(),
                             center.theta = as.numeric(),

                             horizontal.distance = as.numeric(), # repeated line
                             radius = as.numeric(),

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
                        phi.left = tapply(.filteraux$phi.left, .filteraux$tree, mean, na.rm = TRUE),
                        phi.right = tapply(.filteraux$phi.right, .filteraux$tree, mean, na.rm = TRUE),
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
    .filter$filter <- ifelse(.filter$sec == 1.3 & .filter$arc.circ == 1, 1, 0)
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

      .tree.2 <- .tree[ , c("tree", "x", "y"), drop = FALSE]
      .tree.2$tree <- 1:nrow(.tree.2)

      .voro <- data[ , c("x", "y", "z")]
      .voro <- sf::st_as_sf(.voro, coords = c("x", "y"))

      .tree.2 <- sf::st_as_sf(.tree, coords = c("x", "y"))
      .voronoi <- sf::st_collection_extract(sf::st_voronoi(do.call(c, sf::st_geometry(.tree.2))))
      id <- sf::st_intersects(.voronoi, .tree.2)

      # .voro <- sf::st_intersection(.voro, .voronoi)
      .voro$tree <- unlist(sf::st_intersects(.voro, .voronoi))

      freq <- table(.voro$tree)


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
    if(is.null(data$id) & is.null(.tree$v.com)){

      .tree <- .tree[, c("tree", "x", "y", "phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
      colnames(.tree) <- c("tree", "x", "y", "phi", "phi.left", "phi.right", "h.dist", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

    } else if (is.null(data$id) & !is.null(.tree$v.com)) {

      .tree <- .tree[, c("tree", "x", "y", "phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
      colnames(.tree) <- c("tree", "x", "y", "phi", "phi.left", "phi.right", "h.dist", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")


    } else if (!is.null(data$id) & is.null(.tree$v.com)) {

      # If plot identification (id) is available

      .tree$id <- data$id[1]
      .tree$file <- data$file[1]

      .tree <- .tree[, c("id", "file", "tree", "x", "y", "phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
      colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "phi.left", "phi.right", "h.dist", "dbh", "h", "v", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

    } else {

      # If plot identification (id) is available

      .tree$id <- data$id[1]
      .tree$file <- data$file[1]

      .tree <- .tree[, c("id", "file", "tree", "x", "y", "phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion"), drop = FALSE]
      colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "phi.left", "phi.right", "h.dist", "dbh", "h", "v", "v.com", "n.pts", "n.pts.red", "n.pts.est", "n.pts.red.est", "partial.occlusion")

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
