
# This function assigns detected sections to their stem axis

.stem.assignment.single.scan <- function(data, stem.section, breaks){


    kk <- data.frame(sec = table(data$sec[data$arc.circ == 1]))


    kk[, 1] <- as.numeric(as.character(kk[, 1]))
    kk <- kk[kk[,2] == max(kk[,2]), ]


    if(nrow(kk) > 1){
      kk$dif <- abs(kk[, 1] - 1.3)
      kk <- kk[kk$dif == min(kk$dif), ]}


    if(nrow(kk) > 1){
      kk <- kk[kk[, 1] == min(kk[, 1]), ]}


    if(kk[, 1] > 1){
      eje <- data[data$sec < kk[, 1] + 1 & data$sec > kk[, 1] - 1,
                    c("center.x", "center.y", "center.rho", "center.phi", "radius", "sec")]} else {

      eje <- data[data$sec < as.numeric(as.character(kk[, 1])) + 1,
                    c("center.x", "center.y", "center.rho", "center.phi", "radius", "sec")]}


    .dbscan <- dbscan::dbscan(eje[, c("center.x", "center.y"), drop = FALSE],
                              eps = mean(eje$radius[eje$sec == as.character(kk[, 1])], na.rm = TRUE), minPts = 1)

    eje$tree <- .dbscan$cluster
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

      eje.2 <- rbind(eje.2, kk)}

    return(eje.2)

}


.stem.assignment.single.scan.2 <- function(filter, eje, stem.section, x.center, y.center, single.tree){

  eje$sec <- as.character(eje$sec)

  eje$phi <- atan2(eje$y - y.center, eje$x - x.center)
  eje$phi <- ifelse(eje$phi < 0, eje$phi + (2 * pi), eje$phi)
  eje$rho <- sqrt((eje$x - x.center) ^ 2 + (eje$y - y.center) ^ 2)

  filter$cluster <- 1:nrow(filter)
  filter$sec <- as.character(filter$sec)
  filter <- merge(eje, filter, by = "sec", all.y = TRUE)
  filter$dist <- sqrt((filter$center.x - filter$x) ^ 2 + (filter$center.y - filter$y) ^ 2)
  filter$dist.rho <- abs(filter$center.rho - filter$rho)
  filter$dist.phi <- abs(filter$center.phi - filter$phi)
  filter$sec <- as.numeric(filter$sec)

  filter <- filter[!is.na(filter$dist) | !is.na(filter$dist.rho) | !is.na(filter$dist.phi), ]

  filteraux <- data.frame(tree = as.numeric(), sec = as.numeric(),
                          center.x = as.numeric(), center.y = as.numeric(),
                          center.phi = as.numeric(), center.rho = as.numeric(),
                          center.r = as.numeric(), center.theta = as.numeric(),
                          radius = as.numeric(),
                          n.pts = as.numeric(), n.pts.red = as.numeric(),
                          circ = as.numeric(), arc.cir = as.numeric())

  if(!is.null(single.tree)){

    for (i in unique(filter$tree)) {

      .filt.tree <- filter[filter$tree == i, ]

      for (j in unique(.filt.tree$sec)) {

        .filt <- .filt.tree[.filt.tree$sec == j, ]

        .filt$dist.total <- .filt$dist / max(.filt$dist)

        .filt$dist.total <- ifelse(is.nan(.filt$dist.total), 0, .filt$dist.total)

        .filt <- .filt[.filt$dist.total == min(.filt$dist.total), ]

        if(nrow(.filt) > 1)
          .filt <- .filt[1, ]

        if(.filt$dist > 2.5){next}

        .filt <- .filt[, c("tree", "cluster", "sec", "dist",
                           "center.x", "center.y", "center.phi",
                           "center.rho", "center.r", "center.theta",
                           "radius", "n.pts", "n.pts.red",
                           "phi.left",  "phi.right", "arc.circ")]

        filteraux <- rbind(filteraux, .filt)

        filter <- filter[filter$cluster != .filt$cluster, ]

      }

    }

  } else {

    for (i in unique(filter$tree)) {

      .filt.tree <- filter[filter$tree == i, ]

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


        .filt <- .filt[, c("tree", "cluster", "sec", "dist",
                           "center.x", "center.y", "center.phi",
                           "center.rho", "center.r", "center.theta",
                           "radius", "n.pts", "n.pts.red",
                           "phi.left",  "phi.right", "arc.circ")]

        filteraux <- rbind(filteraux, .filt)

        filter <- filter[filter$cluster != .filt$cluster, ]

      }

    }

  }

  # Removing repeated cluster

  filteraux <- filteraux[!duplicated(filteraux$cluster), ]


  filteraux <- filteraux[order(filteraux$tree, filteraux$sec), , drop = FALSE]


  filter <- data.frame(tree = as.numeric(), sec = as.numeric(), dist = as.numeric(),
                       center.x = as.numeric(), center.y = as.numeric(),
                       center.phi = as.numeric(), center.rho = as.numeric(),
                       center.r = as.numeric(), center.theta = as.numeric(),
                       radius = as.numeric(),
                       n.pts = as.numeric(), n.pts.red = as.numeric(),
                       circ = as.numeric(), arc.cir = as.numeric())

  for (i in unique(filteraux$tree)) {

    .filt <- filteraux[filteraux$tree == i, ]

    .filt <- .filt[.filt$dist < mean(.filt$dist, na.rm = TRUE) + mean(.filt$radius, na.rm = TRUE), ]


    if(nrow(.filt) < 1){next}


    if(nrow(.filt) == 1){

      .filt$dif <- NA
      filter <- rbind(filter, .filt)
      next

    }


    .filt$dif <- c(diff(.filt$radius), diff(.filt$radius)[length(diff(.filt$radius))])
    .filt$dif.sec <- c(abs(diff(.filt$sec)), abs(diff(.filt$sec))[length(abs(diff(.filt$sec)))])
    .filt$dif <- ifelse(.filt$dif.sec > 1, .filt$dif * .filt$dif.sec, .filt$dif)

    threshold <- mean(.filt$radius) / 2

    .filt <- .filt[.filt$dif < threshold & .filt$dif > - threshold, ]

    if(nrow(.filt) < 1){next}

    if(nrow(.filt) == 1){

      .filt <- .filt[, -ncol(.filt)]
      filter <- rbind(filter, .filt)
      next

    }


    while (max(.filt$dif) >= threshold / 2 & min(.filt$dif) <= - threshold) {

      .filt <- .filt[.filt$dif < threshold / 2 & .filt$dif > - threshold, ]


      .filt$dif <- c(diff(.filt$radius), diff(.filt$radius)[length(diff(.filt$radius))])
      .filt$dif.sec <- c(abs(diff(.filt$sec)), abs(diff(.filt$sec))[length(abs(diff(.filt$sec)))])
      .filt$dif <- ifelse(.filt$dif.sec > 1, .filt$dif * .filt$dif.sec, .filt$dif)


    }


    .filt <- .filt[, -ncol(.filt)]
    filter <- rbind(filter, .filt)



  }

  return(filter)

}


# This function assigns detected sections to their stem axis

.stem.assignment.multi.scan <- function(.filter, eje, stem.section, x.center, y.center, single.tree){


  # eje$sec <- as.character(eje$sec)

  eje$phi <- atan2(eje$y - y.center, eje$x - x.center)
  eje$phi <- ifelse(eje$phi < 0, eje$phi + (2 * pi), eje$phi)
  eje$rho <- sqrt((eje$x - x.center) ^ 2 + (eje$y - y.center) ^ 2)


  .filter$cluster <- 1:nrow(.filter)
  # .filter$sec <- as.character(.filter$sec)
  .filter <- merge(eje, .filter, by = "sec", all.y = TRUE)
  .filter$dist <- sqrt((.filter$center.x - .filter$x) ^ 2 + (.filter$center.y - .filter$y) ^ 2)
  .filter$dist.rho <- abs(.filter$center.rho - .filter$rho)
  .filter$dist.phi <- abs(.filter$center.phi - .filter$phi)
  # .filter$sec <- as.numeric(.filter$sec)


  .filter <- .filter[!is.na(.filter$dist) | !is.na(.filter$dist.rho) | !is.na(.filter$dist.phi), ]

  filteraux <- data.frame(tree = as.numeric(), sec = as.numeric(),
                          center.x = as.numeric(), center.y = as.numeric(),
                          center.phi = as.numeric(), center.rho = as.numeric(),
                          center.r = as.numeric(), center.theta = as.numeric(),
                          radius = as.numeric(),
                          n.pts = as.numeric(), n.pts.red = as.numeric(),
                          circ = as.numeric(), arc.cir = as.numeric())

  if(!is.null(single.tree)){

    for (i in unique(.filter$tree)) {

      .filt.tree <- .filter[.filter$tree == i, ]

      for (j in unique(.filt.tree$sec)) {

        .filt <- .filt.tree[.filt.tree$sec == j, ]

        .filt$dist.total <- .filt$dist / max(.filt$dist)

        .filt$dist.total <- ifelse(is.nan(.filt$dist.total), 0, .filt$dist.total)

        .filt <- .filt[.filt$dist.total == min(.filt$dist.total), ]

        if(nrow(.filt) > 1)
          .filt <- .filt[1, ]

        if(.filt$dist > 2.5){next}

        .filt <- .filt[, c("tree", "cluster", "sec", "dist",
                           "center.x", "center.y", "center.phi",
                           "center.rho", "center.r", "center.theta",
                           "radius", "n.pts", "n.pts.red", "circ", "arc.circ")]

        filteraux <- rbind(filteraux, .filt)

        .filter <- .filter[.filter$cluster != .filt$cluster, ]

      }

    }

  } else {


    for (i in unique(.filter$tree)) {

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


        .filt <- .filt[, c("tree", "cluster", "sec", "dist",
                           "center.x", "center.y", "center.phi",
                           "center.rho", "center.r", "center.theta",
                           "radius", "n.pts", "n.pts.red", "circ", "arc.circ")]

        filteraux <- rbind(filteraux, .filt)

        .filter <- .filter[.filter$cluster != .filt$cluster, ]

      }

    }

  }

  # Removing repeated cluster

  filteraux <- filteraux[!duplicated(filteraux$cluster), ]


  filteraux <- filteraux[order(filteraux$tree, filteraux$sec), , drop = FALSE]


  .filter <- data.frame(tree = as.numeric(), sec = as.numeric(), dist = as.numeric(),
                        center.x = as.numeric(), center.y = as.numeric(),
                        center.phi = as.numeric(), center.rho = as.numeric(),
                        center.r = as.numeric(), center.theta = as.numeric(),
                        radius = as.numeric(),
                        n.pts = as.numeric(), n.pts.red = as.numeric(),
                        circ = as.numeric(), arc.cir = as.numeric())

  for (i in unique(filteraux$tree)) {

    .filt <- filteraux[filteraux$tree == i, ]


    .filt <- .filt[.filt$dist < mean(.filt$dist, na.rm = TRUE) + mean(.filt$radius, na.rm = TRUE), ]



    if(nrow(.filt) < 1){next}


    if(nrow(.filt) == 1){

      .filt$dif <- NA
      .filter <- rbind(.filter, .filt)
      next

      }


    .filt$dif <- c(diff(.filt$radius), diff(.filt$radius)[length(diff(.filt$radius))])
    .filt$dif.sec <- c(abs(diff(.filt$sec)), abs(diff(.filt$sec))[length(abs(diff(.filt$sec)))])
    .filt$dif <- ifelse(.filt$dif.sec > 1, .filt$dif * .filt$dif.sec, .filt$dif)

    threshold <- mean(.filt$radius) / 2

    .filt <- .filt[.filt$dif < threshold & .filt$dif > - threshold, ]

    if(nrow(.filt) < 1){next}

    if(nrow(.filt) == 1){

        .filt <- .filt[, -ncol(.filt)]
        .filter <- rbind(.filter, .filt)
        next

      }


    while (max(.filt$dif) >= threshold / 2 & min(.filt$dif) <= -threshold) {

      .filt <- .filt[.filt$dif < threshold / 2 & .filt$dif > -threshold, ]


      .filt$dif <- c(diff(.filt$radius), diff(.filt$radius)[length(diff(.filt$radius))])
      .filt$dif.sec <- c(abs(diff(.filt$sec)), abs(diff(.filt$sec))[length(abs(diff(.filt$sec)))])
      .filt$dif <- ifelse(.filt$dif.sec > 1, .filt$dif * .filt$dif.sec, .filt$dif)


    }


    .filt <- .filt[, -ncol(.filt)]
    .filter <- rbind(.filter, .filt)


  }

  return(.filter)

}



# Function to fit a circle through three points
.fit_circle <- function(points) {

  x <- points[, 1]
  y <- points[, 2]

  D <- 2 * (x[1] * (y[2] - y[3]) + x[2] * (y[3] - y[1]) + x[3] * (y[1] - y[2]))

  Ux <- ((x[1]^2 + y[1]^2) * (y[2] - y[3]) + (x[2]^2 + y[2]^2) * (y[3] - y[1]) + (x[3]^2 + y[3]^2) * (y[1] - y[2])) / D
  Uy <- ((x[1]^2 + y[1]^2) * (x[3] - x[2]) + (x[2]^2 + y[2]^2) * (x[1] - x[3]) + (x[3]^2 + y[3]^2) * (x[2] - x[1])) / D

  center <- c(Ux, Uy)
  radius <- sqrt((x[1] - Ux)^2 + (y[1] - Uy)^2)

  return(list(center = center, radius = radius))
}



.RANSAC <- function(data){

  dat <- data

  dist <- 0.05

  inliers <- dat[sample(nrow(dat), 3), ]

  # fit <- .fit_circle(inliers)
  fit <- fit_circle_cpp(inliers)

  radio <- sqrt((dat[, "x"] - fit$center[1]) ^ 2 + (dat[, "y"] - fit$center[2]) ^ 2)
  dat <- cbind(dat, radio)
  inliers <- abs(fit$radius-dat[, "radio"])
  dat <- cbind(dat, inliers)
  dat <- dat[dat[, "inliers"] < dist, ]

  if(length(dat) == 0 | nrow(dat) < 2){

    out <- data.frame(x = as.numeric(), y = as.numeric(),
                      radio = as.numeric(), n = as.numeric(), mae = as.numeric(), cv = as.numeric())

    return(out)}


  dat <- dat[, c("x", "y")]

  inliers <- dat[sample(nrow(dat), 3), ]

  # fit <- .fit_circle(inliers)
  fit <- fit_circle_cpp(inliers)

  radio <- sqrt((dat[, "x"] - fit$center[1]) ^ 2 + (dat[, "y"] - fit$center[2]) ^ 2)
  dat <- cbind(dat, radio)
  inliers <- abs(fit$radius-dat[, "radio"])
  dat <- cbind(dat, inliers)
  dat <- dat[dat[, "inliers"] < dist, ]

  if(length(dat) == 0 | nrow(dat) < 2){

    out <- data.frame(x = as.numeric(), y = as.numeric(),
                      radio = as.numeric(), n = as.numeric(), mae = as.numeric(), cv = as.numeric())

    return(out)}

  mae <- abs(sum(coef(fit)[1]-dat[, 3])) / nrow(dat)
  cv <- stats::sd(raster::pointDistance(dat[, c("x", "y")], c(fit$center[1], fit$center[2]), lonlat = FALSE)) / fit$radius

  out <- data.frame(x = fit$center[1], y = fit$center[2], radio = fit$radius, n = nrow(dat), mae = mae, cv = cv)

  return(out)

}


