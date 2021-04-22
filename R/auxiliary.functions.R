
# Compute several weighted mean functions for a numeric vector

.wmean.calculation <- function(data, w, mean.names) {

  sapply(mean.names,
         function(x, data, w) {
           get(paste("weighted_mean", x, sep = "_"))(data, w)
         },
         data = data, w = w)

}


# Compute number of decimals places

.decimals <- function(x){

  x <- format(x, scient = FALSE)
  ifelse(!base::grepl(".", x, fix = TRUE), 0,
         nchar(base::strsplit(x, ".", fix = TRUE)[[1]][2]))

}


# Custom rounding of numbers

.customCeiling <- function(x, Decimals = 1) {

  x2 <- x * 10 ^ Decimals
  ceiling(x2) / 10 ^ Decimals

}

.customFloor <- function(x, Decimals = 1) {

  x2 <- x * 10 ^ Decimals
  floor(x2) / 10 ^ Decimals

}


# Custom format for numbers

.format.numb <- function(x, dec) {

  format(x, trim = TRUE, nsmall = dec)

}


# Select last value of a vector

.select.last.value <- function(x) {

  x <- x[!is.na(x)]
  x <- ifelse(length(x) == 0, NA, x[length(x)])

}


# Order trees' database by horizontal distance, and compute variables/metrics:
# density, basal area, volume, mean diameters, mean heights, and/or BAF

.metrics.calculation <- function(tree, fixed.area = TRUE, k.tree = TRUE,
                                 angle.count = TRUE, mean.names = NULL,
                                 case = NULL) {

  # Order trees by horizontal distance
  tree <- tree[order(tree[, "horizontal.distance"], decreasing = F), ,
               drop = FALSE]

  # Variables/metrics for fixed area plot and/or k-tree plot
  if (fixed.area | k.tree) {

    # Density (trees/plot)
    tree <- cbind(tree, n = 1:nrow(tree))

    # Basal area (cross-sectional area of trees at 1.3 m from the ground),
    # accumulated from ntree = 1 to ntree = Ntree (m2/ha)
    .g <- (pi / 4) * tree[, "dbh"] ^ 2
    tree <- cbind(tree, G.acum = cumsum(.g))

    if (!is.null(case)) {

      # Height values
      .column.height <- switch(case, tls =  tree[, "P99"],
                               field = tree[, "total.height"])

      # Volume, accumulated from ntree = 1 to ntree = Ntree (m3/ha)

      # Coeficient of 0.45
      # .v <- (pi / 4) * tree[, "dbh"] ^ 2 * .column.height * 0.45

      # Paraboloid
      .v <- pi * (.column.height ^ 2 / 2) * ((tree[, "dbh"] / 2) ^ 2 /
                                               (.column.height - 1.3) ^ 2)

      tree <- cbind(tree, V.acum = cumsum(.v))

      # Mean diameters and heights
      if (!is.null(mean.names)) {

        # Auxiliary function to compute mean diameters and heights
        .wmean <-
          function(n, data, mean.names) {
            t(sapply(n,
                     function(n, data, mean.names) {
                       .wmean.calculation(mean.names = mean.names,
                                          data = data[1:n], w = rep(1, n))
                     },
                     data = data, mean.names = mean.names))
          }

        # Mean diameters (m)
        .mean.names <- mean.names[substr(names(mean.names),1, 1) == "d"]
        .d <- .wmean(n = tree[, "n"], data = tree[, "dbh"],
                     mean.names = .mean.names)
        tree <- cbind(tree, .d)

        # Mean heights (m)
        .mean.names <- mean.names[substr(names(mean.names),1, 1) == "h"]
        .h <- .wmean(n = tree[, "n"], data = .column.height,
                     mean.names = .mean.names)
        tree <- cbind(tree, .h)

      }

    }

  }

  # Variables and metrics for angle-count plot
  if (angle.count) {

    # BAF threshold for each tree
    .BAF <- 2500 / (tree[, "horizontal.distance"] / tree[, "dbh"]) ^ 2
    tree <- cbind(tree, BAF = .BAF)

  }

  return(tree)

}


# Compute a radius sequence, select trees according to argument 'radius.max',
# compute accumulated number of points, and create a data.frame containing the
# trees' data for each radius value - FIXED AREA PLOTS

.radius.fixed.area.calculation <- function(radius.min, radius.increment,
                                           radius.max, tree, case = NULL,
                                           num.dec = 1) {

  # Create a radius sequence according to 'radius.min', 'radius.increment', and
  # 'radius.max' arguments
  .radius.seq <- seq(from = radius.max, to = radius.min,
                     by = - radius.increment)
  .radius.seq <- sort(unique(round(.radius.seq, num.dec)))

  # Select trees according to the value radius.max
  tree <- tree[tree[, "horizontal.distance"] <= radius.max, , drop = FALSE]

  # Compute minimum radius which include each tree
  .radius <- .radius.seq[sapply(tree[, "horizontal.distance"],
                                function(radius, radius.seq){
                                  which(radius<=radius.seq)[1]
                                },
                                radius.seq = .radius.seq)]
  tree <- cbind (radius = .radius, tree)


  # Calculate accumulated number of points
  if (!is.null(case) && case == "tls") {

    tree[, "num.points"] <- cumsum(tree[, "num.points"])
    tree[, "num.points.est"] <- cumsum(tree[, "num.points.est"])
    tree[, "num.points.hom"] <- cumsum(tree[, "num.points.hom"])
    tree[, "num.points.hom.est"] <- cumsum(tree[, "num.points.hom.est"])

  }


  # Create a data.frame containing the trees' data corresponding to each value
  # in the previous radius sequence
  tree <-
    lapply(1:length(.radius.seq),
           function(n, x, y, data) {
             ind <- y == x[n]
             if (sum(ind) > 0) result <- data[ind, , drop = FALSE]
             else {
               result <- data.frame(matrix(NA, nrow = 1, ncol = ncol(data),
                                           dimnames = list(NULL,
                                                           colnames(data))),
                                    stringsAsFactors = FALSE)
               result[, "radius"] <- as.numeric(x[n])
             }
             return(result)
           },
           x = .format.numb(x = .radius.seq, dec = num.dec),
           y = .format.numb(x = tree[, "radius"], dec =num.dec), data = tree)
  tree <- do.call(rbind, tree)

  # Fill missing values in the data.frame using previous values

  # Tree
  .tree <- tree %>% tidyr::fill("tree", .direction = "down")
  tree$tree <- .tree$tree

  # Density (trees/ha)
  .n <- tree %>% tidyr::fill("n", .direction = "down")
  tree$n <- .n$n

  # Basal area (m2/ha)
  .G.acum <- tree %>% tidyr::fill("G.acum", .direction = "down")
  tree$G.acum <- .G.acum$G.acum

  if (is.null(case)) {

    # Stratum
    .stratum<- tree %>% tidyr::fill("stratum", .direction = "down")
    tree$stratum <- .stratum$stratum

  } else {

    # Volume (m3/ha)
    .V.acum <- tree %>% tidyr::fill("V.acum", .direction = "down")
    tree$V.acum <- .V.acum$V.acum

    # Dbh (m)
    .dbh <- tree %>% tidyr::fill("dbh", .direction = "down")
    tree$dbh <- .dbh$dbh

    # Mean diameters (cm)
    .d <- tree %>% tidyr::fill("d", .direction = "down")
    tree$d <- .d$d

    .dg <- tree %>% tidyr::fill("dg", .direction = "down")
    tree$dg <- .dg$dg

    .dgeom <- tree %>% tidyr::fill("dgeom", .direction = "down")
    tree$dgeom <- .dgeom$dgeom

    .dharm <- tree %>% tidyr::fill("dharm", .direction = "down")
    tree$dharm <- .dharm$dharm

    # Mean heights (m)
    .h <- tree %>% tidyr::fill("h", .direction = "down")
    tree$h <- .h$h

    .hg <- tree %>% tidyr::fill("hg", .direction = "down")
    tree$hg <- .hg$hg

    .hgeom <- tree %>% tidyr::fill("hgeom", .direction = "down")
    tree$hgeom <- .hgeom$hgeom

    .hharm <- tree %>% tidyr::fill("hharm", .direction = "down")
    tree$hharm <- .hharm$hharm

    if (case == "tls") {

      # Number of points
      .num.points <- tree %>% tidyr::fill("num.points", .direction = "down")
      tree$num.points <- .num.points$num.points

      .num.points.hom <- tree %>% tidyr::fill("num.points.hom",
                                              .direction = "down")
      tree$num.points.hom <- .num.points.hom$num.points.hom

      .num.points.est <- tree %>% tidyr::fill("num.points.est",
                                              .direction = "down")
      tree$num.points.est <- .num.points.est$num.points.est

      .num.points.hom.est <- tree %>% tidyr::fill("num.points.hom.est",
                                                  .direction = "down")
      tree$num.points.hom.est <- .num.points.hom.est$num.points.hom.est

    }

  }

  return(tree)

}


# Compute radius for each k-tree plot, select trees according to argument
# 'k.tree', and compute accumulated number of points - K-TREE PLOTS

.radius.k.tree.calculation <- function(k.tree.max, tree, case = NULL) {

  # Compute radius for each k-tree
  .dist2 <- c(tree[-1, "horizontal.distance"],
              utils::tail(tree[, "horizontal.distance"], n = 1))
  tree <- cbind(radius = (tree[, "horizontal.distance"] + .dist2) / 2, tree)

  # Select trees according to the value k.tree.max
  tree <- tree[1:k.tree.max, , drop = FALSE]

  # Calculate accumulated number of points
  if (!is.null(case) && case == "tls") {

    tree[, "num.points"] <- cumsum(tree[, "num.points"])
    tree[, "num.points.est"] <- cumsum(tree[, "num.points.est"])
    tree[, "num.points.hom"] <- cumsum(tree[, "num.points.hom"])
    tree[, "num.points.hom.est"] <- cumsum(tree[, "num.points.hom.est"])

  }

  return(tree)

}


# Compute expansion factors, and estimate variables per ha - FIXED AREA PLOTS
# AND K-TREE PLOTS

.fixed.area.k.tree.calculation <- function(data, distance.sampling = NULL,
                                           case = NULL){


  # Compute expansion factor, expansion factor with distance-sampling based
  # correction (if 'distance.sampling' is not NULL), and expansion factor with
  # correction under spetial consideration of shadowing


  # Compute expansion factor for each radius (to convert units per plot into
  # units per ha)
  data <- cbind(data, EF = 10000 / (pi * data[, "radius"] ^ 2))

  if (!is.null(case) && case == "tls") {

    # Compute expansion factor for each radius using Distance-Sampling based
    # correction (Astrup et al., 2014), if 'distance.sampling' is not NULL
    if (!is.null(distance.sampling)) {

      .col.names <- c("P.hn", "P.hr", "P.hn.cov", "P.hr.cov")
      distance.sampling <-
        distance.sampling[match(data[, "tree"], distance.sampling[, "tree"]),
                          .col.names, drop = FALSE]
      colnames(distance.sampling) <- gsub("P", "EF", .col.names, fixed = TRUE)
      distance.sampling <- 10000 /
        (distance.sampling * pi * data[, "radius"] ^ 2)
      data <- cbind(data, distance.sampling)

    }

    # Compute expansion factor for each radius with correction under spetial
    # consideration of shadowing (Seidel and Ammer, 2014)

    data <- cbind(data, EF.sh = NA)

    for (.j in unique(data[, "radius"])) {

      # Compute shadows area
      .col.names <- c("horizontal.distance", "dbh", "partial.occlusion",
                      "wide")
      .shadow <- as.matrix(cbind(data[data[, "radius"] <= .j, .col.names,
                                      drop = FALSE],
                                 shadow = NA))

      .shadow[, "shadow"] <- ifelse(
        .shadow[, "partial.occlusion"] == 0,

        # Non-occluded trees
        ((((pi * .j ^ 2) - (pi * .shadow[, "horizontal.distance"] ^ 2)) /
            (2 * pi)) * atan2(.shadow[, "dbh"],
                              .shadow[, "horizontal.distance"])) -
          ((pi * (.shadow[, "dbh"] / 2) ^ 2) / 2),

        # Occluded trees
        ((((pi * .j ^ 2) - (pi * .shadow[, "horizontal.distance"] ^ 2)) /
            (2 * pi)) * (.shadow[, "wide"])) -
          (((pi * (.shadow[, "dbh"] / 2) ^ 2) *
              .shadow[, "wide"]) / atan2(.shadow[, "dbh"],
                                         .shadow[, "horizontal.distance"])))

      # Correction: shadow = 0 if tree is not completely inside the plot
      .shadow[.shadow[, "horizontal.distance"] + .shadow[, "dbh"] / 2 > .j,
              "shadow"] <- 0

      # Compute expansion factor
      data[data[, "radius"] == .j, "EF.sh"] <-
        1 + (sum(.shadow[, "shadow"], na.rm = TRUE) / (pi * .j ^ 2))

    }

    data[, "EF.sh"] <- data[, "EF"] * data[, "EF.sh"]

  }


  # Estimate variables per ha: density (trees/ha), basal area (m2/ha) and
  # volume (m3/ha)


  # Select columns containing expansion factors
  .col.names <- colnames(data)[substr(colnames(data), 1, 2) == "EF"]

  # Density (trees/ha)
  .N <- data[, "n"] * data[, .col.names, drop = FALSE]
  colnames(.N) <- gsub("EF", "N", colnames(.N), fixed = TRUE)

  # Basal area (m2/ha)
  .G <- data[, "G.acum"] * data[, .col.names, drop = FALSE]
  colnames(.G) <- gsub("EF", "G", colnames(.G), fixed = TRUE)

  if (!is.null(case)) {

    # Volume (m3/ha)
    .V <- data[, "V.acum"] * data[, .col.names, drop = FALSE]
    colnames(.V) <- gsub("EF", "V", colnames(.V), fixed = TRUE)

    # Mean diameters (cm), and mean heights (m)
    .col.names <- c("d", "dg", "dgeom", "dharm", "h", "hg", "hgeom", "hharm")
    .dh <- data[, .col.names, drop = FALSE]

    if (case == "tls"){

      colnames(.N)[colnames(.N) == "N"] <- paste("N", case, sep = ".")
      colnames(.G)[colnames(.G) == "G"] <- paste("G", case, sep = ".")
      colnames(.V)[colnames(.V) == "V"] <- paste("V", case, sep = ".")
      colnames(.dh) <- paste(colnames(.dh), case, sep = ".")

    }

  }


  # Final matrix with stratum, radius, k, N (trees/ha), number of points,
  # G (m2/ha), V (m3/ha), mean diameters (cm), and/or mean heights (m)
  .Plot <- matrix(nrow = nrow(data), ncol = 0)
  if (is.null(case)) .Plot <- cbind(.Plot, data[, "stratum", drop = FALSE])
  .Plot <- cbind(.Plot, data[, "radius", drop = FALSE], k = data[, "n"], .N)
  if (!is.null(case) && case == "tls") {

    .Plot <- cbind(.Plot, data[, c("num.points", "num.points.hom",
                                   "num.points.est", "num.points.hom.est"),
                               drop = FALSE])

  }
  .Plot <- cbind(.Plot, .G)
  if (!is.null(case)) .Plot <- cbind(.Plot, .V, .dh)
  rownames(.Plot) <- NULL

  return(.Plot)

}


# Compute correction of occlusion, estimate variables per ha, and compute
# mean diameters and heights - ANGLE-COUNT
# PLOTS

.angle.count.calculation <- function(BAF, data, mean.names = NULL,
                                     case = NULL){

  # Select data according to BAF
  data <- data[data[, "BAF"] >= BAF, , drop = FALSE]


  # Estimate variables per ha: density (trees / ha), basal area (m2 / ha) and
  # volume (m3 / ha)

  # Density (trees / ha)
  .N <- cbind(N = BAF / (pi * (data[, "dbh"] / 2) ^ 2))

  # Basal area (m2 / ha)
  .G <- cbind(G = rep(BAF, nrow(data)))

  if (!is.null(case)) {

    # Volume (m3 / ha)
    .column.height <- switch(case, tls = data[, "P99"],
                             field = data[, "total.height"])
    .V <- cbind(V = pi * (data[, "dbh"] / 2) ^ 2 * .column.height * 0.45 *
                  .N[, "N"])

    if (case == "tls") {

      # Correction of occlusion - Strahler et al. (2008)
      .De <- mean(data[, "dbh"]) *
        (1 + (stats::sd(data[, "dbh"], na.rm = TRUE) /
                mean(data[, "dbh"], na.rm = TRUE)) ^ 2) ^ 0.5
      .t <- ((.N[, "N"] / 10000) * .De * data[, "dbh"]) /
        (2 * sqrt(BAF / 10000))
      .Ft <- (2 / .t ^ 2) * (1 - exp(-.t) * (1 + .t))

      # Density (trees / ha), basal area (m2 / ha) and volume (m3 / ha) with
      # correction of occlusion
      .N <- cbind(.N, N.pam = .N[, "N"] / .Ft)
      colnames(.N)[colnames(.N) == "N"] <- paste("N", case, sep = ".")
      .G <- cbind(.G, G.pam = .G[, "G"] / .Ft)
      colnames(.G)[colnames(.G) == "G"] <- paste("G", case, sep = ".")
      .V <- cbind(.V, V.pam = .V[, "V"] / .Ft)
      colnames(.V)[colnames(.V) == "V"] <- paste("V", case, sep = ".")

    }

    # Mean diameters and heights
    if (!is.null(mean.names)) {

      # Mean diameters (m)
      .mean.names <- mean.names[substr(names(mean.names),1, 1) == "d"]
      .d <- .wmean.calculation(data = data[, "dbh"],
                               w = rep(1, nrow(data)),
                               mean.names = .mean.names)

      # Mean heights (m)
      .mean.names <- mean.names[substr(names(mean.names),1, 1) == "h"]
      .h <- .wmean.calculation(data = .column.height, w = rep(1, nrow(data)),
                               mean.names = .mean.names)

      if (case == "tls") {

        names(.d) <- paste(names(.d), case, sep = ".")
        names(.h) <- paste(names(.h), case, sep = ".")

      }

    }

  }


  # Final vector with BAF, N (trees/ha), number of points, G (m2/ha),
  # V (m3/ha), mean diameters (cm), and/or mean heights (m)
  .Plot <- c()
  if (is.null(case)) .Plot <- c(.Plot, data[1, "stratum"])
  .Plot <- c(.Plot, BAF = BAF, apply(.N, 2, sum))
  if (!is.null(case) && case == "tls") {

    .Plot <- c(.Plot, apply(data[, c("num.points", "num.points.hom",
                                     "num.points.est", "num.points.hom.est"),
                                 drop = FALSE], 2, sum))
  }
  .Plot <- c(.Plot, apply(.G, 2, sum))
  if (!is.null(case)) .Plot <- c(.Plot, apply(.V, 2, sum), .d, .h)
  rownames(.Plot) <- NULL

  return(.Plot)

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


# .pol create polygons as regular squares overlapped 0.05

.pol <- function(data, dist){

  .xpol <- c(data$x.1-dist, data$x.2+dist, data$x.3+dist, data$x.4-dist)
  .ypol <- c(data$y.1-dist, data$y.2-dist, data$y.3+dist, data$y.4+dist)

  return(sp::Polygons(list(sp::Polygon(cbind(.xpol,.ypol))), ID = data$id))

}


# This functions apply calculated ncr index for all points.
# For that purpose, voxelize point cloud by means of regular
# grid in x any z coordinates

.ncr.remove.slice.double <- function(data){

  # Select necessary fields from original txt file of point cloud

  .data <- data[,c("point", "x", "y", "z")]


  # Create x and y coordinates for grid

  .x = seq(min(.data$x), max(.data$x))
  .y = seq(min(.data$y), max(.data$y))

  # Empty data frame where coordinates neccesaries for
  # creating grid will be saved

  .grid <- data.frame(id = as.character(),
                      x.1 = as.numeric(), x.2 = as.numeric(),
                      x.3 = as.numeric(), x.4 = as.numeric(),
                      y.1 = as.numeric(), y.2 = as.numeric(),
                      y.3 = as.numeric(), y.4 = as.numeric())


  # Fill the empty data frame .grid created just before

  for (i in 1:(length(.x)-1)) {
    for (j in 1:(length(.y)-1)) {

      .out <- data.frame(id = paste("x", i, "y", j, sep="."),
                         x.1 = .x[i], x.2 = .x[i+1], x.3 = .x[i+1], x.4 = .x[i],
                         y.1 = .y[j], y.2 = .y[j], y.3 = .y[j+1], y.4 = .y[j+1])

      .grid <- rbind(.grid, .out)

    }
  }

  .row.names <- .grid$id

  # Split .grid by id into a list for using lapply function and
  # applying .pol function over all elements of the list

  .grid <- split(.grid, .grid$id)

  # Apply .pol function to every element of the list and create
  # and object SpatialPolygons

  .grid <- sp::SpatialPolygons(lapply(.grid, .pol, dist = 0.05))


  # Generate and SpatialPoints object to extract those points
  # which are avor the polygons of .grid

  .pts <- .data[, c("x", "y")]
  # dimnames(.pts)[[1]] <- c(.data$point)
  .pts <- sp::SpatialPoints(.pts)

  .attributes <- .data[, c("point", "x", "y", "z")]

  .pts = sp::SpatialPointsDataFrame(.pts, .attributes)

  .attributes <- data.frame(row.names = .row.names)

  .grid = sp::SpatialPolygonsDataFrame(.grid, .attributes)

  .pts <- sp::over(.grid, .pts, returnList = TRUE)

  .dat <- lapply(.pts, as.matrix, ncol = 4)
  .dat <- .dat[names(which(lapply(.dat, length) > 4))]
  # .dat <- .dat[names(which(lapply(.dat, length) < 40000))]


  .ncr <- do.call(rbind, lapply(.dat, ncr_point_cloud_double))

  .ncr <- .ncr[which(.ncr$ncr > 0 & .ncr$ncr < 9999), ]

  .data <- merge(data, .ncr, by = "point", all = TRUE)

  .data <- .data[!duplicated(.data), ]

  return(.data)

}


