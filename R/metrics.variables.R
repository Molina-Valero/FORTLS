
metrics.variables <- function(tree.list.tls, distance.sampling = NULL,
                              plot.parameters, dir.data = NULL,
                              save.result = TRUE, dir.result = NULL) {

  # Check if plot.parameter argument is included, otherwise, function will stop
  if(is.null(plot.parameters))
    stop('No plot design specified in plot.parameters argument')


  # Define working directory as directory by default for loading/saving files
  if (is.null(dir.data)) dir.data <- getwd()
  if (is.null(dir.result)) dir.result <- getwd()


  # Convert dbh (cm) to International System of Units (m)
  tree.list.tls$dbh <- tree.list.tls$dbh / 100


  # Define values for certain plot parameters, and create empty data.frames
  # where results will be saved


  # Define an auxiliary list containing metrics names

  # Define expansion factors/correction of occlusion names
  .ef.names <- character()
  if (!is.null(distance.sampling))
    .ef.names <- c(.ef.names, "hn", "hr", "hn.cov", "hr.cov")
  .ef.names <- list(fixed.area = c(.ef.names, "sh"),
                    k.tree = c(.ef.names, "sh"), angle.count = "pam")

  # Define mean names
  .mean.names <- c(d = "arit", dg = "sqrt", dgeom = "geom", dharm = "harm",
                   h = "arit", hg = "sqrt", hgeom = "geom", hharm = "harm")

  # Define probabilities for percentiles calculation
  .prob = c(1, 5, 10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90, 95, 99)/100

  # Define list containing metrics names for fixed area, k-tree and angle-count
  # plots
  .metrics.names <- sapply(c("fixed.area", "k.tree", "angle.count"),
                           function(x, y, z) {
                             c(
                               # Density (trees/ha)
                               c("N.tls", paste("N", y[[x]], sep = ".")),

                               # Number of points
                               "num.points", "num.points.est", "num.points.hom",
                               "num.points.hom.est",

                               # Basal area (m2/ha)
                               c("G.tls", paste("G", y[[x]], sep = ".")),

                               # Volume (m3/ha)
                               c("V.tls", paste("V", y[[x]], sep = ".")),

                               # Mean diameters (cm), and mean heights (m)
                               paste(names(z), "tls", sep = "."),

                               # Mean dominant diameters (cm), and mean dominant
                               # heights (m)
                               paste(names(z), "0", "tls", sep = "."),

                               # Height percentiles (m)
                               sprintf("P%02i", 100 * .prob))
                           },
                           y = .ef.names, z = .mean.names, simplify = FALSE)

  # Create an empty data.frame where results will be saved for fixed area plots
  fixed.area.plot <- NULL
  if (!is.null(plot.parameters$radius)) {

    # Create an empty data.frame
    .col.names <- c("stratum", "id", "radius", .metrics.names$fixed.area)
    fixed.area.plot <- data.frame(matrix(numeric(), ncol = length(.col.names),
                                         dimnames = list(NULL, .col.names)),
                                  stringsAsFactors = FALSE)

  }

  # Create an empty data.frame where results will be saved for k-tree plots
  k.tree.plot <- NULL
  if (!is.null(plot.parameters$k.tree)) {

    # Create an empty data.frame
    .col.names <- c("stratum", "id", "k", .metrics.names$k.tree)
    k.tree.plot <- data.frame(matrix(numeric(), ncol = length(.col.names),
                                     dimnames = list(NULL, .col.names)),
                              stringsAsFactors = FALSE)

  }

  # Create an empty data.frame where results will be saved for angle-count plots
  angle.count.plot <- NULL
  if (!is.null(plot.parameters$BAF)) {

    # Create an empty data.frame
    .col.names <- c("stratum", "id", "BAF", .metrics.names$angle.count)
    angle.count.plot <- data.frame(matrix(numeric(), ncol = length(.col.names),
                                          dimnames = list(NULL, .col.names)),
                                   stringsAsFactors = FALSE)

  }

  # Define number of trees per plot for dominant metrics
  if (is.null(plot.parameters$num.trees)) {

    .num <- rep(100, nrow(plot.parameters))

  } else {

    .num <- plot.parameters$num.trees

  }
  names(.num) <- rownames(plot.parameters)


  # Define TXT files list

  # Obtain initial TXT files list
  .files <- unique(tree.list.tls[, "file"])

  # Check if TXT files exist
  .files.exists <- file.exists(file.path(dir.data, .files))
  if (all(!.files.exists))
    warning("None of the TXT files in 'tree.list.tls' argument is available ",
            "in 'dir.data', so no computation will be done")
  else if (any(!.files.exists))
    warning(sum(!.files.exists), " TXT file(s) in 'tree.list.tls' argument ",
            "is(are) missing in 'dir.data'. This(these) plot(s) was(were) not ",
            "taken into account during the execution")
  .files <- .files[.files.exists]


  # Loop for each TLS plot
  for (.i in .files) {

    # Select TLS plot stratum
    .stratum <- ifelse(!is.null(tree.list.tls$stratum),
                       as.character(unique(
                         tree.list.tls[tree.list.tls$file == .i, "stratum"])),
                       rownames(plot.parameters)[1])


    # Select TLS plot id
    .id <- unique(tree.list.tls[tree.list.tls$file == .i, "id"])

    # Define initial time, and print message with the id of plot being
    # processing
    t0 <- Sys.time()
    message("Computing metrics for plot: ", .id)


    # Select 'distance.sampling' info for plot .id ----


    .distSampling <- distance.sampling
    if (!is.null(.distSampling)) {

      .distSampling <-
        distance.sampling$tree[distance.sampling$tree$id == .id, , drop = FALSE]

    }


    # Create points' database and trees' database for plot .id ----


    # Read the points' database for the TLS plot from the file .i
    .data.tls <- suppressMessages(
      vroom::vroom(file.path(dir.data, .i),
                   col_select = c("x", "y", "z", "rho"), progress = FALSE)
    )
    .data.tls <- as.data.frame(.data.tls, stringsAsFactors = FALSE)

    # Select data corresponding to the TLS plot from the trees' database
    .tree.tls <- tree.list.tls[tree.list.tls$file == .i, , drop = FALSE]

    # Voronoi tessellation
    .voro <- .tree.tls[ , c("tree", "x", "y"), drop = FALSE]
    .voro <- ggvoronoi::voronoi_polygon(.voro, x = "x", y = "y", outline = NULL,
                                        data.frame = FALSE)

    sp::coordinates(.data.tls) <- ~ x + y + z
    .voro <- sp::over(.data.tls, .voro)
    .voro <- sp::SpatialPointsDataFrame(.data.tls, .voro)

    .voro <- as.data.frame(.voro, stringsAsFactors = FALSE)
    .voro <- .voro[!is.na(.voro$tree), , drop = FALSE]

    # Order .data.tls by rho, select columns required for height percentiles
    # calculations below, and convert to matrix
    .data.tls <- as.data.frame(.data.tls, stringsAsFactors = FALSE)
    .data.tls <- .data.tls[order(.data.tls$rho, decreasing = F), c("z", "rho"),
                           drop = FALSE]
    .data.tls <- as.matrix(.data.tls)
    rownames(.data.tls) <- NULL

    # Compute height percentile P99
    .P99 <- sapply(sort(unique(.voro$tree)),
                   function(tree, voro) {
                     z <- voro$z[voro$tree == tree]
                     P99 <-
                       height_perc_cpp(rho_seq = Inf, z = z, rho = z)[, "P99"]
                     names(P99) <- tree
                     return(P99)
                   },
                   voro =.voro)
    .P99 <- data.frame(tree = names(.P99), P99 = .P99)
    .tree.tls <- merge(.tree.tls, .P99, by = "tree", all = FALSE)

    # Compute angular aperture
    .wide <- .tree.tls$phi.right - .tree.tls$phi.left
    .wide <- ifelse(.wide < 0, (2 * pi) + .wide, .wide)
    .tree.tls <- cbind(.tree.tls, wide = .wide)

    # Select only columns required for calculations below
    .col.names <- c("tree", "horizontal.distance", "dbh", "num.points",
                    "num.points.hom", "num.points.est", "num.points.hom.est",
                    "partial.occlusion", "wide", "P99")
    .tree.tls <- .tree.tls[ , .col.names, drop = FALSE]
    rownames(.tree.tls) <- NULL

    # Order by horizontal distance, and compute variables/metrics: density,
    # basal area, volume, mean diameters, mean heights, and BAF
    .tree.tls <-
      .metrics.calculation(tree = .tree.tls,
                           fixed.area = !is.null(plot.parameters$radius),
                           k.tree = !is.null(plot.parameters$k.tree),
                           angle.count = !is.null(plot.parameters$BAF),
                           mean.names = .mean.names, case = "tls")


    # Compute fixed area plot metrics for plot .id ----


    if (!is.null(plot.parameters$radius)) {

      # Define number of decimals places to be considered
      .num.dec <- .decimals(plot.parameters[.stratum, "radius"])

      # Define minimum and maximum value for radius sequence according to
      # 'plot.parameters' argument, and horizontal distances in trees' database
      .radius.min <- min(.customCeiling(.tree.tls[, "horizontal.distance"],
                                        Decimals = .num.dec))
      .radius.max <- plot.parameters[.stratum, "radius"]
      if (.radius.min > .radius.max) {

        .radius.max <- .radius.min
        warning("For plot ", .id, ", radius was increased to ", .radius.min,
                " to ensure that at least one tree is included")

      }
      .radius.min <- .radius.max

      # Compute a radius sequence, select trees according to maximum radius,
      # compute accumulated number of points, and create a data.frame containing
      # the trees' data for each radius value
      .fixedAreaPlot <-
        .radius.fixed.area.calculation(radius.min = .radius.min,
                                       radius.increment = 10 ^ (-.num.dec),
                                       radius.max = .radius.max,
                                       tree = .tree.tls, case = "tls",
                                       num.dec = .num.dec)

      # Compute expansion factors, and estimate variables per ha
      .fixedAreaPlot <-
        .fixed.area.k.tree.calculation(data = .fixedAreaPlot,
                                       distance.sampling = .distSampling,
                                       case = "tls")

      # Select last row for each radius value
      .fixedAreaPlot <- matrix(apply(.fixedAreaPlot, 2, tapply,
                                     .fixedAreaPlot[, "radius"],
                                     .select.last.value),
                               ncol = ncol(.fixedAreaPlot),
                               dimnames = list(NULL, colnames(.fixedAreaPlot)))
      rownames(.fixedAreaPlot) <- .format.numb(x = .fixedAreaPlot[, "radius"],
                                               dec = .num.dec)

      # Filter plot by maximum radius
      .row.names <- .format.numb(x = .radius.max, dec = .num.dec)
      .fixedAreaPlot <- .fixedAreaPlot[.row.names, , drop = FALSE]

      # Compute mean dominant diameters and heights for each radius value
      .dh.0 <- fixed_area_cpp(radius_seq = .fixedAreaPlot[, "radius"],
                              hdist = .tree.tls[, "horizontal.distance"],
                              d = .tree.tls[, "dbh"], h = .tree.tls[, "P99"],
                              num = .num[.stratum])
      rownames(.dh.0) <- .format.numb(x = .dh.0[, "radius"], dec = .num.dec)
      .dh.0 <- as.matrix(.dh.0[, colnames(.dh.0) != "radius", drop = FALSE])
      colnames(.dh.0) <- paste(colnames(.dh.0), "tls", sep = ".")
      .fixedAreaPlot <- cbind(.fixedAreaPlot,
                              .dh.0[rownames(.fixedAreaPlot), , drop = FALSE])

      # Compute height percentiles for each radius value
      .perc <- height_perc_cpp(rho_seq = .fixedAreaPlot[, "radius"],
                               z = .data.tls[, "z"], rho = .data.tls[, "rho"])
      rownames(.perc) <- .format.numb(x = .perc[, "rho_seq"], dec = .num.dec)
      .perc <- as.matrix(.perc[, colnames(.perc) != "rho_seq", drop = FALSE])
      .fixedAreaPlot <- cbind(.fixedAreaPlot,
                              .perc[rownames(.fixedAreaPlot), , drop = FALSE])

      # Convert diameters from International System of Units (m) to cm
      .col.names <- names(.mean.names)[substr(names(.mean.names), 1, 1) == "d"]
      .col.names <- c(.col.names, paste(.col.names, "0", sep = "."))
      .col.names <- paste(.col.names, "tls", sep = ".")
      .fixedAreaPlot[, .col.names] <- .fixedAreaPlot[, .col.names] * 100

      # Save fixed area plot results
      .row.names <- paste(.id, rownames(.fixedAreaPlot), sep = "-")
      .fixedAreaPlot <- data.frame(stratum = .stratum, id = .id, .fixedAreaPlot,
                                   row.names = .row.names,
                                   stringsAsFactors = FALSE)
      .fixedAreaPlot <- .fixedAreaPlot[, colnames(fixed.area.plot),
                                       drop = FALSE]
      fixed.area.plot <- rbind(fixed.area.plot, .fixedAreaPlot)

    }


    # Compute k-tree plot metrics for plot .id ----


    if (!is.null(plot.parameters$k.tree)) {

      # Define number of decimals places to be considered
      .num.dec <- 0

      # Define maximum number of trees according to 'plot.parameters' argument,
      # and number of trees' database
      .k.tree.max <- nrow(.tree.tls)
      if (.k.tree.max < plot.parameters[.stratum, "k.tree"]) {

        warning("For plot ", .id, ", k-tree was reduced to ", .k.tree.max,
                " since it is the number of trees in the plot")

      } else {

        .k.tree.max <- round(plot.parameters[.stratum, "k.tree"])
        if (.k.tree.max != plot.parameters[.stratum, "k.tree"])
          warning("k-tree was rounded to ", .k.tree.max)

      }

      # Compute radius for each k-tree plot, select trees according to
      # maximum number of trees, and compute accumulated number of points
      .kTreePlot <- .radius.k.tree.calculation(k.tree.max = .k.tree.max,
                                               tree = .tree.tls, case = "tls")

      # Compute expansion factors, and estimate variables per ha
      .kTreePlot <-
        .fixed.area.k.tree.calculation(data = .kTreePlot,
                                       distance.sampling = .distSampling,
                                       case = "tls")
      rownames(.kTreePlot) <- .format.numb(x = .kTreePlot[, "k"],
                                           dec = .num.dec)

      # Filter plot by maximum k-tree
      .row.names <- .format.numb(x = .k.tree.max, dec = .num.dec)
      .kTreePlot <- .kTreePlot[.row.names, , drop = FALSE]

      # Create radius sequence for computing dominant diameters and heights,
      # and height percentiles for each k value
      .rad.seq <- .kTreePlot[, "radius"]
      names(.rad.seq) <- .format.numb(x = .kTreePlot[, "k"], dec = .num.dec)
      .rad.seq <- sapply(1:length(.rad.seq),
                         function(n, x) {
                           x.max <- max(x[1:n])
                           names(x.max) <- names(x)[n]
                           return(x.max)
                         },
                         x = .rad.seq)
      # Compute mean dominant diameters and heights for each k value
      .dh.0 <- k_tree_cpp(k_seq = .kTreePlot[, "k"], radius_seq = .rad.seq,
                          k = .tree.tls[, "n"], d = .tree.tls[, "dbh"],
                          h = .tree.tls[, "P99"], num = .num[.stratum])
      rownames(.dh.0) <- .format.numb(x = .dh.0[, "k"], dec = .num.dec)
      .dh.0 <- as.matrix(.dh.0[, colnames(.dh.0) != "k", drop = FALSE])
      colnames(.dh.0) <- paste(colnames(.dh.0), "tls", sep = ".")
      .kTreePlot <- cbind(.kTreePlot,
                          .dh.0[rownames(.kTreePlot), , drop = FALSE])

      # Compute height percentiles for each k value
      .perc <- height_perc_cpp(rho_seq = .rad.seq, z = .data.tls[, "z"],
                               rho = .data.tls[, "rho"])
      rownames(.perc) <- names(.rad.seq)
      .perc <- as.matrix(.perc[, colnames(.perc) != "rho_seq", drop = FALSE])
      .kTreePlot <- cbind(.kTreePlot,
                          .perc[rownames(.kTreePlot), , drop = FALSE])

      # Convert diameters from International System of Units (m) to cm
      .col.names <- names(.mean.names)[substr(names(.mean.names), 1, 1) == "d"]
      .col.names <- c(.col.names, paste(.col.names, "0", sep = "."))
      .col.names <- paste(.col.names, "tls", sep = ".")
      .kTreePlot[, .col.names] <- .kTreePlot[, .col.names] * 100

      # Save k-tree plot results
      .row.names <- paste(.id, rownames(.kTreePlot), sep = "-")
      .kTreePlot <- data.frame(stratum = .stratum, id = .id, .kTreePlot,
                               row.names = .row.names, stringsAsFactors = FALSE)
      .kTreePlot <- .kTreePlot[, colnames(k.tree.plot), drop = FALSE]
      k.tree.plot <- rbind(k.tree.plot, .kTreePlot)

    }


    # Compute angle-count plot metrics for plot .id ----


    if (!is.null(plot.parameters$BAF)) {

      # Define number of decimals places to be considered
      .num.dec <- .decimals(plot.parameters[.stratum, "BAF"])

      # Define minimum and maximum value for BAF sequence according to
      # 'plot.parameters' argument, and BAF values in trees' database
      .BAF.rang <- range(.customFloor(.tree.tls[, "BAF"], Decimals = .num.dec))
      .BAF.min <- .BAF.rang[1]
      .BAF.max <- plot.parameters[.stratum, "BAF"]
      if (.BAF.rang[2] < .BAF.max) {

        .BAF.max <- .BAF.rang[2]
        warning("For plot ", .id, ", BAF was reduced to ", .BAF.max,
                " to ensure that at least one tree is included")

      }
      .BAF.min <- .BAF.max

      # Compute a BAF sequence for angle-count plots
      .BAF.seq <- seq(from = .BAF.max, to = .BAF.min, by = - 10 ^ (-.num.dec))
      .BAF.seq <- sort(unique(round(.BAF.seq, .num.dec)))
      names(.BAF.seq) <- .format.numb(x = .BAF.seq, dec = .num.dec)

      # Compute correction of occlusion, estimate variables per ha, and
      # compute mean diameters and heights
      .angleCountPlot <- lapply(.BAF.seq, .angle.count.calculation,
                                data = .tree.tls, mean.names = .mean.names,
                                case = "tls")
      .angleCountPlot <- do.call(rbind, .angleCountPlot)
      rownames(.angleCountPlot) <- .format.numb(x = .angleCountPlot[, "BAF"],
                                                dec = .num.dec)

      # Compute mean dominant diameters and heights for each BAF value
      .BAF.order <- order(.tree.tls[, "BAF"], decreasing = TRUE)
      .dh.0 <- angle_count_cpp(baf_seq = .BAF.seq,
                               baf = .tree.tls[.BAF.order, "BAF"],
                               d = .tree.tls[.BAF.order, "dbh"],
                               h = .tree.tls[.BAF.order, "P99"],
                               num = .num[.stratum])
      rownames(.dh.0) <- .format.numb(x = .dh.0[, "BAF"], dec = .num.dec)
      .dh.0 <- as.matrix(.dh.0[, colnames(.dh.0) != "BAF", drop = FALSE])
      colnames(.dh.0) <- paste(colnames(.dh.0), "tls", sep = ".")
      .angleCountPlot <- cbind(.angleCountPlot,
                               .dh.0[rownames(.angleCountPlot), , drop = FALSE])

      # Compute height percentiles for each BAF value
      .rho.seq <- sapply(.BAF.seq,
                         function(BAF, tree){
                           max(tree[tree[, "BAF"] >= BAF,
                                    "horizontal.distance"])
                         },
                         tree = .tree.tls)
      .perc <- height_perc_cpp(rho_seq = .rho.seq, z = .data.tls[, "z"],
                               rho = .data.tls[, "rho"])
      rownames(.perc) <- names(.BAF.seq)
      .perc <- as.matrix(.perc[, colnames(.perc) != "rho_seq", drop = FALSE])
      .angleCountPlot <- cbind(.angleCountPlot,
                               .perc[rownames(.angleCountPlot), , drop = FALSE])

      # Convert diameters from International System of Units (m) to cm
      .col.names <- names(.mean.names)[substr(names(.mean.names), 1, 1) == "d"]
      .col.names <- c(.col.names, paste(.col.names, "0", sep = "."))
      .col.names <- paste(.col.names, "tls", sep = ".")
      .angleCountPlot[, .col.names] <- .angleCountPlot[, .col.names] * 100

      # Save angle-count plot results
      .row.names <- paste(.id, rownames(.angleCountPlot), sep = "-")
      .angleCountPlot <- data.frame(stratum = .stratum, id = .id,
                                    .angleCountPlot, row.names = .row.names,
                                    stringsAsFactors = FALSE)
      .angleCountPlot <- .angleCountPlot[, colnames(angle.count.plot),
                                         drop = FALSE]
      angle.count.plot <- rbind(angle.count.plot, .angleCountPlot)

    }


    # Define final time, and print message
    t1 <- Sys.time()
    message(" (", format(round(difftime(t1, t0, units = "secs"), 2)), ")")

  }


  # Save fixed area, k-tree and/or angle-count plot results, and write csv files
  # containing them if 'save.result' is TRUE


  if (!is.null(plot.parameters$radius)) {

    fixed.area.plot <- fixed.area.plot[order(fixed.area.plot[, "id"],
                                             fixed.area.plot[, "radius"]), ,
                                       drop = FALSE]
    rownames(fixed.area.plot) <- NULL
    if (isTRUE(save.result))
      utils::write.csv(fixed.area.plot,
                       file =
                         file.path(dir.result,
                                   "metrics.variables.fixed.area.plot.csv"),
                       row.names = FALSE)

  }

  if (!is.null(plot.parameters$k.tree)) {

    k.tree.plot <- k.tree.plot[order(k.tree.plot[, "id"], k.tree.plot[, "k"]), ,
                               drop = FALSE]
    rownames(k.tree.plot) <- NULL
    if (isTRUE(save.result))
      utils::write.csv(k.tree.plot,
                       file =
                         file.path(dir.result,
                                   "metrics.variables.k.tree.plot.csv"),
                       row.names = FALSE)

  }

  if (!is.null(plot.parameters$BAF)) {

    angle.count.plot <- angle.count.plot[order(angle.count.plot[, "id"],
                                               angle.count.plot[, "BAF"]), ,
                                         drop = FALSE]
    rownames(angle.count.plot) <- NULL
    if (isTRUE(save.result))
      utils::write.csv(angle.count.plot,
                       file =
                         file.path(dir.result,
                                   "metrics.variables.angle.count.plot.csv"),
                       row.names = FALSE)

  }


  return(list(fixed.area.plot = fixed.area.plot, k.tree.plot = k.tree.plot,
              angle.count.plot = angle.count.plot))

}
