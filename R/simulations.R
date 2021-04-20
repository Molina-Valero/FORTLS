
simulations <- function(tree.list.tls, distance.sampling = NULL,
                        tree.list.field,
                        plot.parameters = list(radius.max = 25, k.tree.max = 50,
                                               BAF.max = 4),
                        dir.data = NULL, save.result = TRUE,
                        dir.result = NULL) {


  # Check arguments


  # 'tree.list.tls' must be a data.frame with at least one row, certain
  # mandatory columns and most of them numeric
  .col.names <-  c("id", "file", "tree", "x", "y", "phi.left", "phi.right",
                   "horizontal.distance", "dbh", "num.points", "num.points.hom",
                   "num.points.est", "num.points.hom.est", "partial.occlusion")
  if (!is.data.frame(tree.list.tls))
    stop("'tree.list.tls' must be a data.frame")
  if (nrow(tree.list.tls) < 1)
    stop("'tree.list.tls' must have at least one row")
  for (.i in .col.names){

    if (!.i %in% colnames(tree.list.tls))
      stop("'tree.list.tls' must have a column named '", .i,"'")
    if (!.i %in% c("id", "file") && !is.numeric(tree.list.tls[, .i]))
      stop("Column '", .i,"' of 'tree.list.tls' must be numeric")

  }

  # 'distance.sampling' must be NULL (by default) or a list with at least an
  # element named 'tree', which must be a data.frame with at least one row,
  # certain mandatory columns and most of them numeric
  if (!is.null(distance.sampling)) {

    .col.names <-  c("id", "tree", "P.hn", "P.hn.cov", "P.hr", "P.hr.cov")
    if (!is.list(distance.sampling)) stop("'distance.sampling' must be a list")
    if (!"tree" %in% names(distance.sampling))
      stop("'distance.sampling' must have an element named 'tree'")
    if (!is.data.frame(distance.sampling$tree))
      stop("'distance.sampling$tree' must be a data.frame")
    if (nrow(distance.sampling$tree) < 1)
      stop("'distance.sampling$tree' must have at least one row")
    for (.i in .col.names){

      if (!.i %in% colnames(distance.sampling$tree))
        stop("'distance.sampling$tree' must have a column named '", .i,"'")
      if (.i != "id" && !is.numeric(distance.sampling$tree[, .i]))
        stop("Column '", .i, "' of 'distance.sampling$tree' must be numeric")

    }
    if (any(!unique(tree.list.tls[, "id"]) %in%
            unique(distance.sampling$tree[, "id"])))
      warning("There is any plot identification in 'tree.list.tls' argument ",
              "which is missing in 'distance.sampling' argument")
    for (.i in unique(tree.list.tls[, "id"])) {

      if (any(!unique(tree.list.tls[tree.list.tls[, "id"] == .i, "tree"]) %in%
              unique(distance.sampling$tree[distance.sampling$tree[, "id"] ==
                                            .i, "tree"])))
        warning("There is any tree numbering for plot ", .id, " in ",
                "'tree.list.tls' argument which is missing in ",
                "'distance.sampling' argument")

    }

  }

  # 'tree.list.field' must be a data.frame with at least one row, certain
  # mandatory columns and most of them numeric
  .col.names <-  c("id", "tree", "dbh", "horizontal.distance", "total.height",
                   "dead")
  if (!is.data.frame(tree.list.field))
    stop("'tree.list.field' must be a data.frame")
  if (nrow(tree.list.field) < 1)
    stop("'tree.list.field' must have at least one row")
  for (.i in .col.names){

    if (!.i %in% colnames(tree.list.field))
      stop("'tree.list.field' must have a column named '", .i,"'")
    if (.i != "id" && !is.numeric(tree.list.field[, .i]))
      stop("Column '", .i,"' of 'tree.list.field' must be numeric")

  }
  if (any(!unique(tree.list.tls[, "id"]) %in% unique(tree.list.field[, "id"])))
    stop("There is any plot identification in 'tree.list.tls' argument ",
         "which is missing in 'tree.list.field' argument")

  # 'plot.parameters' must be list(radius.max = 25, k.tree.max = 50,
  # BAF.max = 4) (by default) or a list with all or any of the numeric and
  # strictly positive preset parameters, including at least one of the
  # following: radius.max, k.tree.max or BAF.max
  .col.names <-  c("radius.max" = TRUE, "radius.increment" = FALSE,
                   "k.tree.max" = TRUE, "BAF.max" = TRUE,
                   "BAF.increment" = FALSE, "num.trees" = FALSE)
  if (!is.list(plot.parameters)) stop("'plot.parameters' must be a list")
  if (is.null(names(plot.parameters)) ||
      all(!names(.col.names)[.col.names] %in% names(plot.parameters)))
    stop("'plot.parameters' must have at least one of the following elements:",
         "'radius.max', 'k.tree.max' or 'BAF.max'")
  if (any(!names(plot.parameters) %in% names(.col.names))) {

    plot.parameters <-
      plot.parameters[names(plot.parameters) %in% names(.col.names)]
    warning("There is any element in 'plot.parameters' which do not ",
            "correspond with preset ones. It was not taken into account ",
            "during the execution")

  }
  for (.i in names(plot.parameters)){

    if (!is.numeric(plot.parameters[[.i]]))
      stop("Element '", .i, "' of 'plot.parameters' must be numeric")
    if (length(plot.parameters[[.i]]) > 1) {

      plot.parameters[[.i]] <- plot.parameters[[.i]][1]
      warning("Only first value in element '", .i, "' of 'plot.parameters' ",
              "was taken into account during the execution")

    }
    if (plot.parameters[[.i]] <= 0)
      stop("Element '", .i, "' of 'plot.parameters' must be strictly positive")

  }

  # 'dir.data' must be NULL (by default) or a character string containing the
  # absolute path to a existing directory
  if (!is.null(dir.data)) {

    if (!is.character(dir.data)) stop("'dir.data' must be a character string")
    if (length(dir.data) > 1) {

      dir.data <- dir.data[1]
      warning("Only first value in 'dir.data' was taken into account during ",
              "the execution")

    }
    if (!dir.exists(dir.data)) stop("'dir.data' directory must exist")

  } else {

    # Define working directory as directory by default for loading files
    dir.data <- getwd()

  }

  # 'save.result' must be TRUE (by default) or a logical
  if (!is.logical(save.result)) stop("'save.result' must be a logical")
  if (length(save.result) > 1) {

    save.result <- save.result[1]
    warning("Only first value in 'save.result' was taken into account during ",
            "the execution")

  }

  # If 'save.result' is TRUE, 'dir.result' must be NULL (by default) or a
  # character string containing the absolute path to a existing directory
  if (save.result) {

    if (!is.null(dir.result)) {

      if (!is.character(dir.result))
        stop("'dir.result' must be a character string")
      if (length(dir.result) > 1) {

        dir.result <- dir.result[1]
        warning("Only first value in 'dir.result' was taken into account ",
                "during the execution")

      }
      if (!dir.exists(dir.result)) stop("'dir.result' directory must exist")

    } else {

      # Define working directory as directory by default for saving files
      dir.result <- getwd()

    }

  }


  # Convert dbh (cm) to International System of Units (m)
  tree.list.tls$dbh <- tree.list.tls$dbh / 100
  tree.list.field$dbh <- tree.list.field$dbh / 100


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
  # plots - TLS and field data
  .metrics.names <- list(
    tls = sapply(c("fixed.area", "k.tree", "angle.count"),
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
                 y = .ef.names, z = .mean.names, simplify = FALSE),
    field = c(
      # Density (trees/ha), basal area (m2/ha) and volume (m3/ha)
      "N", "G", "V",

      # Mean diameters (cm), and mean heights (m)
      paste(names(.mean.names), sep = "."),

      # Mean dominant diameters (cm), and mean dominant heights (m)
      paste(names(.mean.names), "0", sep = ".")))

  # Define radius increment, and create a list containing empty data.frames
  # where results will be saved for fixed area plots
  fixed.area.plot <- NULL
  if (!is.null(plot.parameters$radius.max)) {

    # Define radius increment
    .radius.increment <- ifelse(is.null(plot.parameters$radius.increment), 0.1,
                                plot.parameters$radius.increment)
    # Create list containing empty data.frame - TLS and field data
    .col.names <- list(tls = c("id", "radius", .metrics.names$tls$fixed.area),
                       field = c("id", "radius", .metrics.names$field))
    fixed.area.plot <- lapply(.col.names,
                              function(x){
                                data.frame(matrix(numeric(), ncol = length(x),
                                                  dimnames = list(NULL, x)),
                                           stringsAsFactors = FALSE)
                              })

  }

  # Create a list containing empty data.frames where results will be saved for
  # k-tree plots
  k.tree.plot <- NULL
  if (!is.null(plot.parameters$k.tree.max)) {

    # Create list containing empty data.frames - TLS and field data
    .col.names <- list(tls = c("id", "k", .metrics.names$tls$k.tree),
                       field = c("id", "k", .metrics.names$field))
    k.tree.plot <- lapply(.col.names,
                          function(x){
                            data.frame(matrix(numeric(), ncol = length(x),
                                              dimnames = list(NULL, x)),
                                       stringsAsFactors = FALSE)
                          })

  }

  # Define BAF increment, and create a list containing empty data.frames where
  # results will be saved for angle-count plots
  angle.count.plot <- NULL
  if (!is.null(plot.parameters$BAF.max)) {

    # Define BAF increment
    .BAF.increment <- ifelse(is.null(plot.parameters$BAF.increment), 0.1,
                             plot.parameters$BAF.increment)
    # Create list containing empty data.frames - TLS and field data
    .col.names <- list(tls = c("id", "BAF", .metrics.names$tls$angle.count),
                       field = c("id", "BAF", .metrics.names$field))
    angle.count.plot <- lapply(.col.names,
                               function(x){
                                 data.frame(matrix(numeric(), ncol = length(x),
                                                   dimnames = list(NULL, x)),
                                            stringsAsFactors = FALSE)
                               })

  }

  # Define number of trees per plot for dominant metrics
  .num <- ifelse(is.null(plot.parameters$num.trees), 100,
                 plot.parameters$num.trees)


  # Create trees' database - Field data


  # Select field plots with available TLS data from the trees' database, and
  # remove dead trees, and trees with missing height or missing dbh
  tree.list.field <-
    tree.list.field[tree.list.field$id %in% unique(tree.list.tls$id) &
                      is.na(tree.list.field$dead) &
                      !is.na(tree.list.field$total.height) &
                      !is.na(tree.list.field$dbh), , drop = FALSE]

  # Select only columns required for calculations below
  tree.list.field <- tree.list.field[, c("id", "tree", "horizontal.distance",
                                         "dbh", "total.height"), drop = FALSE]

  # Define TXT files list
  .files <- unique(tree.list.tls$file)
  .files.exists <- .files %in% list.files(pattern = "txt", path = dir.data)
  if (all(!.files.exists)) {

    warning("None of the TXT files in 'tree.list.tls' is available in ",
            "'dir.data', so no computation will be done")
  }
  else if (any(!.files.exists)) {

    warning(sum(!.files.exists), " TXT files in 'tree.list.tls' are ",
            "missing in 'dir.data', so no computation will be done for them")

  }
  .files <- .files[.files.exists]

  # Loop for each TLS plot
  for (.i in .files) {


    # Select TLS plot id
    .id <- unique(tree.list.tls[tree.list.tls$file == .i, "id"])

    # Define initial time, and print message with the id of plot being
    # processing
    t0 <- Sys.time()
    message("Computing simulations for plot: ", .id)

    # Define list for saving trees' databases - TLS and field data
    .tree <- list(tls = NULL, field = NULL)


    # Select 'distance.sampling' info for plot .id - TLS data ----


    .distSampling <- distance.sampling
    if (!is.null(.distSampling))
      .distSampling <-
      as.matrix(distance.sampling$tree[distance.sampling$tree$id == .id, ,
                                       drop = FALSE])


    # Create points' database and trees' database for plot .id - TLS data ----


    # Read the points' database for the TLS plot from the file .i
    .data.tls <- suppressMessages(
      vroom::vroom(file.path(dir.data, .i),
                   col_select = c("x", "y", "z", "rho"), progress = FALSE)
    )
    .data.tls <- as.data.frame(.data.tls, stringsAsFactors = FALSE)

    # Select data corresponding to the TLS plot from the trees' database
    .tree$tls <- tree.list.tls[tree.list.tls$file == .i, , drop = FALSE]

    # Voronoi tessellation
    .voro <- .tree$tls[ , c("tree", "x", "y"), drop = FALSE]
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
    .tree$tls <- merge(.tree$tls, .P99, by = "tree", all = FALSE)

    # Compute angular aperture
    .wide <- .tree$tls$phi.right - .tree$tls$phi.left
    .wide <- ifelse(.wide < 0, (2 * pi) + .wide, .wide)
    .tree$tls <- cbind(.tree$tls, wide = .wide)

    # Select only columns required for calculations below, and convert to matrix
    .col.names <- c("tree", "horizontal.distance", "dbh", "num.points",
                    "num.points.hom", "num.points.est", "num.points.hom.est",
                    "partial.occlusion", "wide", "P99")
    .tree$tls <- as.matrix(.tree$tls[ , .col.names, drop = FALSE])
    rownames(.tree$tls) <- NULL

    # Order by horizontal distance, and compute variables/metrics: density,
    # basal area, volume, mean diameters, mean heights, and BAF
    .tree$tls <-
      .metrics.calculation(tree = .tree$tls,
                           fixed.area = !is.null(plot.parameters$radius.max),
                           k.tree = !is.null(plot.parameters$k.tree.max),
                           angle.count = !is.null(plot.parameters$BAF.max),
                           mean.names = .mean.names, case = "tls")


    # Create trees' database for plot .id - Field data ----


    # Select only columns required for calculations below, and convert to matrix
    .tree$field <- tree.list.field[tree.list.field$id == .id, , drop = FALSE]
    .col.names <- c("tree", "horizontal.distance", "dbh", "total.height")
    .tree$field <- as.matrix(.tree$field[ , .col.names, drop = FALSE])
    rownames(.tree$field) <- NULL

    # Order by horizontal distance, and compute variables/metrics: density,
    # basal area, volume, mean diameters, mean heights, and BAF
    .tree$field <-
      .metrics.calculation(tree = .tree$field,
                           fixed.area = !is.null(plot.parameters$radius.max),
                           k.tree = !is.null(plot.parameters$k.tree.max),
                           angle.count = !is.null(plot.parameters$BAF.max),
                           mean.names = .mean.names, case = "field")


    # Compute fixed area plot simulations for plot .id ----


    if (!is.null(plot.parameters$radius.max)) {

      # Define number of decimals places to be considered
      .num.dec <- .decimals(.radius.increment)

      # Define minimum and maximum value for radius sequence according to
      # 'plot.parameters' argument, and horizontal distances in trees' databases
      .radius.min <- sapply(.tree,
                            function(tree, num.dec){
                              min(.customCeiling(tree[, "horizontal.distance"],
                                                 Decimals = num.dec))
                            },
                            num.dec = .num.dec)
      .radius.min <- max(.radius.min)
      .radius.max <- plot.parameters$radius.max
      .rho.max <- .customCeiling(max(.data.tls[, "rho"]), Decimals = .num.dec)
      if (.radius.max > .rho.max) {

        .radius.max <- .rho.max
        warning("For plot ", .id, ", 'plot.parameters$radius.max' was reduced ",
                "to ", .radius.max, " since it is the maximum distance in ",
                "point cloud data")

      }
      if (.radius.min > .radius.max) {

        .radius.max <- .radius.min
        warning("For plot ", .id, ", 'plot.parameters$radius.max' was ",
                "increased to ", .radius.min, " to ensure that at least one ",
                "tree is included in all the simulated plots")

      }

      # Loop for TLS and field cases
      for (.j in names(fixed.area.plot)) {

        # Compute a radius sequence, select trees according to maximum radius,
        # compute accumulated number of points, and create a matrix containing
        # the trees' data for each radius value
        .fixedAreaPlot <-
          .radius.fixed.area.calculation(radius.min = .radius.min,
                                         radius.increment = .radius.increment,
                                         radius.max = .radius.max,
                                         tree = .tree[[.j]], case = .j,
                                         num.dec = .num.dec)

        # Compute expansion factors, and estimate variables per ha
        .fixedAreaPlot <- .fixed.area.k.tree.calculation(data = .fixedAreaPlot,
                                                         distance.sampling =
                                                           .distSampling,
                                                         case = .j)

        # Select last row for each radius value
        .fixedAreaPlot <- matrix(apply(.fixedAreaPlot, 2, tapply,
                                       .fixedAreaPlot[, "radius"],
                                       .select.last.value),
                                 ncol = ncol(.fixedAreaPlot),
                                 dimnames = list(NULL,
                                                 colnames(.fixedAreaPlot)))
        rownames(.fixedAreaPlot) <- .format.numb(x = .fixedAreaPlot[, "radius"],
                                                 dec = .num.dec)

        # Compute mean dominant diameters and heights for each radius value
        .dh.0 <- fixed_area_cpp(radius_seq = .fixedAreaPlot[, "radius"],
                                hdist = .tree[[.j]][, "horizontal.distance"],
                                d = .tree[[.j]][, "dbh"],
                                h = .tree[[.j]][, switch(.j, tls = "P99",
                                                         field =
                                                           "total.height")],
                                num = .num)
        rownames(.dh.0) <- .format.numb(x = .dh.0[, "radius"], dec = .num.dec)
        .dh.0 <- as.matrix(.dh.0[, colnames(.dh.0) != "radius", drop = FALSE])
        if (.j == "tls")
          colnames(.dh.0) <- paste(colnames(.dh.0), .j, sep = ".")
        .fixedAreaPlot <- cbind(.fixedAreaPlot,
                                .dh.0[rownames(.fixedAreaPlot), , drop = FALSE])

        # Compute height percentiles for each radius value (only for 'tls' case)
        if (.j == "tls") {

          .perc <- height_perc_cpp(rho_seq = .fixedAreaPlot[, "radius"],
                                   z = .data.tls[, "z"],
                                   rho = .data.tls[, "rho"])
          rownames(.perc) <- .format.numb(x = .perc[, "rho_seq"],
                                          dec = .num.dec)
          .perc <- as.matrix(.perc[, colnames(.perc) != "rho_seq",
                                   drop = FALSE])
          .fixedAreaPlot <- cbind(.fixedAreaPlot,
                                  .perc[rownames(.fixedAreaPlot), ,
                                        drop = FALSE])
        }

        # Convert diameters from International System of Units (m) to cm
        .col.names <- names(.mean.names)[substr(names(.mean.names), 1, 1) ==
                                           "d"]
        .col.names <- c(.col.names, paste(.col.names, "0", sep = "."))
        if (.j == "tls") .col.names <- paste(.col.names, .j, sep = ".")
        .fixedAreaPlot[, .col.names] <- .fixedAreaPlot[, .col.names] * 100

        # Save fixed area plot results
        .row.names <- paste(.id, rownames(.fixedAreaPlot), sep = "-")
        .fixedAreaPlot <- data.frame(id = .id, .fixedAreaPlot,
                                     row.names = .row.names,
                                     stringsAsFactors = FALSE)
        .fixedAreaPlot <- .fixedAreaPlot[, colnames(fixed.area.plot[[.j]]),
                                         drop = FALSE]
        fixed.area.plot[[.j]] <- rbind(fixed.area.plot[[.j]], .fixedAreaPlot)

      }

    }


    # Compute k-tree plot simulations for plot .id ----


    if (!is.null(plot.parameters$k.tree.max)) {

      # Define number of decimals places to be considered
      .num.dec <- 0

      # Define maximum number of trees according to 'plot.parameters' argument,
      # and number of trees' databases
      .k.tree.max <- min(sapply(.tree, nrow))
      if (.k.tree.max < plot.parameters$k.tree.max) {

        warning("For plot ", .id, ", 'plot.parameters$k.tree.max' was reduced ",
                "to ", .k.tree.max, " since it is the number of trees in ",
                "TLS/field plot")

      } else {

        .k.tree.max <- round(plot.parameters$k.tree.max)
        if (.k.tree.max != plot.parameters$k.tree.max)
          warning("'plot.parameters$k.tree.max' was rounded to ", .k.tree.max)

      }

      # Loop for TLS and field cases
      for (.j in names(k.tree.plot)) {

        # Compute radius for each k-tree plot, select trees according to
        # maximum number of trees, and compute accumulated number of points
        .kTreePlot <- .radius.k.tree.calculation(k.tree.max = .k.tree.max,
                                                 tree = .tree[[.j]], case = .j)

        # Compute expansion factors, and estimate variables per ha
        .kTreePlot <- .fixed.area.k.tree.calculation(data = .kTreePlot,
                                                     distance.sampling =
                                                       .distSampling,
                                                     case = .j)
        rownames(.kTreePlot) <- .format.numb(x = .kTreePlot[, "k"],
                                             dec = .num.dec)

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
                            k = .tree[[.j]][, "n"], d = .tree[[.j]][, "dbh"],
                            h = .tree[[.j]][, switch(.j, tls = "P99",
                                                     field = "total.height")],
                            num = .num)
        rownames(.dh.0) <- .format.numb(x = .dh.0[, "k"], dec = .num.dec)
        .dh.0 <- as.matrix(.dh.0[, colnames(.dh.0) != "k", drop = FALSE])
        if (.j == "tls")
          colnames(.dh.0) <- paste(colnames(.dh.0), .j, sep = ".")
        .kTreePlot <- cbind(.kTreePlot,
                            .dh.0[rownames(.kTreePlot), , drop = FALSE])

        # Compute height percentiles for each k value (only for 'tls' case)
        if (.j == "tls") {

          .perc <- height_perc_cpp(rho_seq = .rad.seq, z = .data.tls[, "z"],
                                   rho = .data.tls[, "rho"])
          rownames(.perc) <- names(.rad.seq)
          .perc <- as.matrix(.perc[, colnames(.perc) != "rho_seq",
                                   drop = FALSE])
          .kTreePlot <- cbind(.kTreePlot,
                              .perc[rownames(.kTreePlot), , drop = FALSE])

        }

        # Convert diameters from International System of Units (m) to cm
        .col.names <- names(.mean.names)[substr(names(.mean.names), 1, 1) ==
                                           "d"]
        .col.names <- c(.col.names, paste(.col.names, "0", sep = "."))
        if (.j == "tls") .col.names <- paste(.col.names, .j, sep = ".")
        .kTreePlot[, .col.names] <- .kTreePlot[, .col.names] * 100

        # Save k-tree plot results
        .row.names <- paste(.id, rownames(.kTreePlot), sep = "-")
        .kTreePlot <- data.frame(id = .id, .kTreePlot, row.names = .row.names,
                                 stringsAsFactors = FALSE)
        .kTreePlot <- .kTreePlot[, colnames(k.tree.plot[[.j]]), drop = FALSE]
        k.tree.plot[[.j]] <- rbind(k.tree.plot[[.j]], .kTreePlot)

      }

    }


    # Compute angle-count plot simulations for plot .id ----


    if (!is.null(plot.parameters$BAF.max)) {

      # Define number of decimals places to be considered
      .num.dec <- .decimals(.BAF.increment)

      # Define minimum and maximum value for BAF sequence according to
      # 'plot.parameters' argument, and BAF values in trees' databases
      .BAF.rang <- sapply(.tree,
                          function(tree, num.dec){
                            range(.customFloor(tree[, "BAF"],
                                               Decimals = num.dec))
                          },
                          num.dec = .num.dec)
      .BAF.min <- min(.BAF.rang[1, ])
      .BAF.max <- plot.parameters$BAF.max
      if (min(.BAF.rang[2, ]) < .BAF.max) {

        .BAF.max <- min(.BAF.rang[2, ])
        warning("For plot ", .id, ", 'plot.parameters$BAF.max' was reduced to ",
                .BAF.max, " to ensure that at least one tree is included in ",
                "all the simulated plots")

      }

      if (min(.BAF.rang[1, ]) > .BAF.max) {

        .BAF.max <- max(.BAF.rang[1, ])
        warning("For plot ", .id, ", 'plot.parameters$BAF.max' was increased to ",
                .BAF.max, " to ensure that at least one tree is included in ",
                "all the simulated plots")

      }

      # Loop for TLS and field cases
      for (.j in names(angle.count.plot)) {

        # Compute a BAF sequence for angle-count plots
        .BAF.seq <- seq(from = .BAF.max, to = .BAF.min, by = - .BAF.increment)
        .BAF.seq <- sort(unique(round(.BAF.seq, .num.dec)))
        names(.BAF.seq) <- .format.numb(x = .BAF.seq, dec = .num.dec)

        # Compute correction of occlusion, estimate variables per ha, and
        # compute mean diameters and heights
        .angleCountPlot <- lapply(.BAF.seq, .angle.count.calculation,
                                  data = .tree[[.j]], mean.names = .mean.names,
                                  case = .j)
        .angleCountPlot <- do.call(rbind, .angleCountPlot)
        rownames(.angleCountPlot) <- .format.numb(x = .angleCountPlot[, "BAF"],
                                                  dec = .num.dec)

        # Compute mean dominant diameters and heights for each BAF value
        .BAF.order <- order(.tree[[.j]][, "BAF"], decreasing = TRUE)
        .dh.0 <- angle_count_cpp(baf_seq = .BAF.seq,
                                 baf = .tree[[.j]][.BAF.order, "BAF"],
                                 d = .tree[[.j]][.BAF.order, "dbh"],
                                 h = .tree[[.j]][.BAF.order,
                                                 switch(.j, tls = "P99",
                                                        field =
                                                          "total.height")],
                                 num = .num)
        rownames(.dh.0) <- .format.numb(x = .dh.0[, "BAF"], dec = .num.dec)
        .dh.0 <- as.matrix(.dh.0[, colnames(.dh.0) != "BAF", drop = FALSE])
        if (.j == "tls")
          colnames(.dh.0) <- paste(colnames(.dh.0), .j, sep = ".")
        .angleCountPlot <- cbind(.angleCountPlot,
                                 .dh.0[rownames(.angleCountPlot), ,
                                       drop = FALSE])

        # Compute height percentiles for each BAF value (only for 'tls' case)
        if (.j == "tls") {

          .rho.seq <- sapply(.BAF.seq,
                             function(BAF, tree){
                               max(tree[tree[, "BAF"] >= BAF,
                                        "horizontal.distance"])
                             },
                             tree = .tree[[.j]])
          .perc <- height_perc_cpp(rho_seq = .rho.seq, z = .data.tls[, "z"],
                                   rho = .data.tls[, "rho"])
          rownames(.perc) <- names(.BAF.seq)
          .perc <- as.matrix(.perc[, colnames(.perc) != "rho_seq",
                                   drop = FALSE])
          .angleCountPlot <- cbind(.angleCountPlot,
                                   .perc[rownames(.angleCountPlot), ,
                                         drop = FALSE])

        }

        # Convert diameters from International System of Units (m) to cm
        .col.names <- names(.mean.names)[substr(names(.mean.names), 1, 1) ==
                                           "d"]
        .col.names <- c(.col.names, paste(.col.names, "0", sep = "."))
        if (.j == "tls") .col.names <- paste(.col.names, .j, sep = ".")
        .angleCountPlot[, .col.names] <- .angleCountPlot[, .col.names] * 100

        # Save angle-count plot results
        .row.names <- paste(.id, rownames(.angleCountPlot), sep = "-")
        .angleCountPlot <- data.frame(id = .id, .angleCountPlot,
                                      row.names = .row.names,
                                      stringsAsFactors = FALSE)
        .angleCountPlot <- .angleCountPlot[, colnames(angle.count.plot[[.j]]),
                                           drop = FALSE]
        angle.count.plot[[.j]] <- rbind(angle.count.plot[[.j]], .angleCountPlot)

      }

    }


    # Define final time, and print message
    t1 <- Sys.time()
    message(" (", format(round(difftime(t1, t0, units = "secs"), 2)), ")")

  }


  # Save fixed area, k-tree and/or angle-count plot results, and write csv files
  # containing them if 'save.result' is TRUE


  if (!is.null(plot.parameters$radius.max)) {

    .row.names <-
      rownames(fixed.area.plot$field)[rownames(fixed.area.plot$field) %in%
                                        rownames(fixed.area.plot$tls)]
    .col.names <-
      colnames(fixed.area.plot$tls)[!colnames(fixed.area.plot$tls) %in%
                                      colnames(fixed.area.plot$field)]
    fixed.area.plot <- cbind(fixed.area.plot$field[.row.names, , drop = FALSE],
                             fixed.area.plot$tls[.row.names, .col.names,
                                                 drop = FALSE])
    fixed.area.plot <- fixed.area.plot[order(fixed.area.plot[, "radius"],
                                             fixed.area.plot[, "id"]), ,
                                       drop = FALSE]
    rownames(fixed.area.plot) <- NULL
    if (save.result)
      utils::write.csv(fixed.area.plot,
                       file = file.path(dir.result,
                                        "simulations.fixed.area.plot.csv"),
                       row.names = FALSE)

  }

  if (!is.null(plot.parameters$k.tree.max)) {

    .row.names <- rownames(k.tree.plot$field)[rownames(k.tree.plot$field) %in%
                                                rownames(k.tree.plot$tls)]
    .col.names <- colnames(k.tree.plot$tls)[!colnames(k.tree.plot$tls) %in%
                                              colnames(k.tree.plot$field)]
    k.tree.plot <- cbind(k.tree.plot$field[.row.names, , drop = FALSE],
                         k.tree.plot$tls[.row.names, .col.names, drop = FALSE])
    k.tree.plot <- k.tree.plot[order(k.tree.plot[, "k"], k.tree.plot[, "id"]), ,
                               drop = FALSE]
    rownames(k.tree.plot) <- NULL
    if (save.result)
      utils::write.csv(k.tree.plot,
                       file = file.path(dir.result,
                                        "simulations.k.tree.plot.csv"),
                       row.names = FALSE)

  }

  if (!is.null(plot.parameters$BAF.max)) {

    .row.names <-
      rownames(angle.count.plot$field)[rownames(angle.count.plot$field) %in%
                                         rownames(angle.count.plot$tls)]
    .col.names <-
      colnames(angle.count.plot$tls)[!colnames(angle.count.plot$tls) %in%
                                       colnames(angle.count.plot$field)]
    angle.count.plot <- cbind(angle.count.plot$field[.row.names, ,
                                                     drop = FALSE],
                              angle.count.plot$tls[.row.names, .col.names,
                                                   drop = FALSE])
    angle.count.plot <- angle.count.plot[order(angle.count.plot[, "BAF"],
                                               angle.count.plot[, "id"]), ,
                                         drop = FALSE]
    rownames(angle.count.plot) <- NULL
    if (save.result)
      utils::write.csv(angle.count.plot,
                       file = file.path(dir.result,
                                        "simulations.angle.count.plot.csv"),
                       row.names = FALSE)

  }


  return(list(fixed.area.plot = fixed.area.plot, k.tree.plot = k.tree.plot,
              angle.count.plot = angle.count.plot))

}
