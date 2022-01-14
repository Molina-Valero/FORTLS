
estimation.plot.size <- function(tree.list.tls,
                                 plot.parameters = list(radius.max = 25,
                                                        k.tree.max = 50,
                                                        BAF.max = 4),
                                 average = FALSE, all.plot.designs = FALSE) {


  # Define id as 1 (by default) when tree.list$id is missed

  if (is.null(tree.list.tls$id))
    tree.list.tls$id <- 1

  # Define id as 1 (by default) when tree.list$id is missed

  if (is.null(tree.list.tls$stratum) | all.plot.designs)
    tree.list.tls$stratum <- 1

  # Convert dbh (cm) to International System of Units (m)
  tree.list.tls$dbh <- tree.list.tls$dbh / 100

  # Define character vector containing metrics names for fixed area, k-tree
  # and angle-count plots: density (trees/ha), and basal area (m2/ha)
  .metrics.names <- c("N", "G")


  # Define radius increment, and create an empty data.frame where results will
  # be saved for fixed area plots
  .fixed.area.plot <- NULL
  if (!is.null(plot.parameters$radius.max)) {

    # Define radius increment
    .radius.increment <- 0.1

    # Create an empty data.frame
    .col.names <- c("stratum", "id", "radius", .metrics.names)
    .fixed.area.plot <- data.frame(matrix(numeric(), ncol = length(.col.names),
                                          dimnames = list(NULL, .col.names)),
                                   stringsAsFactors = FALSE)

    .fixed.area.plot.axis <- list(main = "Fixed area plots",
                                  xlab = "Radius (m)")

  }


  # Create an empty data.frame where results will be saved for k-tree plots
  .k.tree.plot <- NULL
  if (!is.null(plot.parameters$k.tree.max)) {

    # Create an empty data.frame
    .col.names <- c("stratum", "id", "k", .metrics.names)
    .k.tree.plot <- data.frame(matrix(numeric(), ncol = length(.col.names),
                                      dimnames = list(NULL, .col.names)),
                               stringsAsFactors = FALSE)

    .k.tree.plot.axis <- list(main = "k-tree plots",
                              xlab = "k-tree (trees)")

  }


  # Define BAF increment, and create an empty data.frame where results will be
  # saved for angle-count plots
  .angle.count.plot <- NULL
  if (!is.null(plot.parameters$BAF.max)) {

    # Define BAF increment
    .BAF.increment <- 0.1

    # Create an empty data.frame
    .col.names <- c("stratum", "id", "BAF", .metrics.names)
    .angle.count.plot <- data.frame(matrix(numeric(), ncol = length(.col.names),
                                           dimnames = list(NULL, .col.names)),
                                    stringsAsFactors = FALSE)

    .angle.count.plot.axis <- list(main = "Angle-count plots",
                                   xlab = expression("BAF" ~ (m^{2}/ha)))

  }


  # Loop for each TLS plot
  for (.i in unique(tree.list.tls$id)) {

    # Select data corresponding to the TLS plot from the trees' database
    .tree.tls <- tree.list.tls[tree.list.tls$id == .i, , drop = FALSE]

    # Select only columns required for calculations below
    .col.names <- c("stratum", "tree", "h.dist", "dbh")
    .tree.tls <- .tree.tls[ , .col.names, drop = FALSE]
    rownames(.tree.tls) <- NULL

    # Order by horizontal distance, and compute variables/metrics: density,
    # basal area, and/or BAF
    .tree.tls <-
      .metrics.calculation(tree = .tree.tls,
                           fixed.area = !is.null(plot.parameters$radius.max),
                           k.tree = !is.null(plot.parameters$k.tree.max),
                           angle.count = !is.null(plot.parameters$BAF.max))


    # Compute fixed area plot simulations for plot .i ----


    if (!is.null(plot.parameters$radius.max)) {

      # Define number of decimals places to be considered
      .num.dec <- .decimals(.radius.increment)

      # Define minimum and maximum value for radius sequence according to
      # 'plot.parameters' argument, and horizontal distances in trees' database
      .radius.min <- min(.customCeiling(.tree.tls[, "h.dist"],
                                        Decimals = .num.dec))
      .radius.max <- min(plot.parameters$radius.max,
                         max(.tree.tls[, "h.dist"]))
      if (.radius.min > .radius.max) .radius.max <- .radius.min

      # Compute a radius sequence, select trees according to maximum radius,
      # and create a data.frame containing the trees' data for each radius value
      .fixedAreaPlot <-
        .radius.fixed.area.calculation(radius.min = .radius.min,
                                       radius.increment = .radius.increment,
                                       radius.max = .radius.max,
                                       tree = .tree.tls, num.dec = .num.dec)

      # Compute expansion factors, and estimate variables per ha
      .fixedAreaPlot <- .fixed.area.k.tree.calculation(data = .fixedAreaPlot)

      # Select last row for each radius value
      .fixedAreaPlot <- matrix(apply(.fixedAreaPlot, 2, tapply,
                                     .fixedAreaPlot[, "radius"],
                                     .select.last.value),
                               ncol = ncol(.fixedAreaPlot),
                               dimnames = list(NULL, colnames(.fixedAreaPlot)))

      # Save fixed area plot results
      .fixedAreaPlot <- data.frame(id = .i, .fixedAreaPlot, row.names = NULL,
                                   stringsAsFactors = FALSE)
      .fixedAreaPlot <- .fixedAreaPlot[, colnames(.fixed.area.plot),
                                       drop = FALSE]
      .fixed.area.plot <- rbind(.fixed.area.plot, .fixedAreaPlot)

    }


    # Compute k-tree plot simulations for plot .i ----


    if (!is.null(plot.parameters$k.tree.max)) {

      # Define number of decimals places to be considered
      .num.dec <- 0

      # Define maximum number of trees according to 'plot.parameters' argument,
      # and number of trees' database
      .k.tree.max <- nrow(.tree.tls)
      if (.k.tree.max >= plot.parameters$k.tree.max)
        .k.tree.max <- round(plot.parameters$k.tree.max)

      # Compute radius for each k-tree plot, and select trees according to
      # maximum number of trees
      .kTreePlot <- .radius.k.tree.calculation(k.tree.max = .k.tree.max,
                                               tree = .tree.tls)

      # Compute expansion factors, and estimate variables per ha
      .kTreePlot <- .fixed.area.k.tree.calculation(data = .kTreePlot)

      # Save k-tree plot results
      .kTreePlot <- data.frame(id = .i, .kTreePlot, row.names = NULL,
                               stringsAsFactors = FALSE)
      .kTreePlot <- .kTreePlot[, colnames(.k.tree.plot), drop = FALSE]
      .k.tree.plot <- rbind(.k.tree.plot, .kTreePlot)

    }


    # Compute angle-count plot simulations for plot .i ----


    if (!is.null(plot.parameters$BAF.max)) {

      # Define number of decimals places to be considered
      .num.dec <- .decimals(.BAF.increment)

      # Define minimum and maximum value for BAF sequence according to
      # 'plot.parameters' argument, and BAF values in trees' database
      .BAF.rang <- range(.customFloor(.tree.tls[, "BAF"], Decimals = .num.dec))
      .BAF.min <- .BAF.rang[1]
      .BAF.max <- plot.parameters$BAF.max
      if (.BAF.rang[2] < .BAF.max) .BAF.max <- .BAF.rang[2]

      # Compute a BAF sequence for angle-count plots
      .BAF.seq <- seq(from = .BAF.max, to = .BAF.min, by = - .BAF.increment)
      .BAF.seq <- sort(unique(round(.BAF.seq, .num.dec)))
      names(.BAF.seq) <- .format.numb(x = .BAF.seq, dec = .num.dec)

      # Compute estimate variables per ha
      .angleCountPlot <- lapply(.BAF.seq, .angle.count.calculation,
                                data = .tree.tls)
      .angleCountPlot <- do.call(rbind, .angleCountPlot)

      # Save angle-count plot results
      .angleCountPlot <- data.frame(id = .i, .angleCountPlot, row.names = NULL,
                                    stringsAsFactors = FALSE)
      .angleCountPlot <- .angleCountPlot[, colnames(.angle.count.plot),
                                         drop = FALSE]
      .angle.count.plot <- rbind(.angle.count.plot, .angleCountPlot)

    }

  }


  # Obtaining mean values and standard deviation when mean argument
  # is specified in the function ----

  if (average | all.plot.designs) {

    .fixed.area.plot.i <- split(.fixed.area.plot, .fixed.area.plot$stratum)

    .agregate <- function(data){

      .out <- data.frame(n = tapply(data$radius, data$radius, length),
                         radius = tapply(data$radius, data$radius, mean,
                                         na.rm = TRUE),
                         N = tapply(data$N, data$radius, mean, na.rm = TRUE),
                         N.sd = tapply(data$N, data$radius, stats::sd,
                                       na.rm = TRUE),
                         G = tapply(data$G, data$radius, mean, na.rm = TRUE),
                         G.sd = tapply(data$G, data$radius, stats::sd,
                                       na.rm = TRUE))

      .out$stratum <- data$stratum[1]

      return(.out)

    }

    .fixed.area.plot <- do.call(rbind, lapply(.fixed.area.plot.i, .agregate))

    .k.tree.plot.i <- split(.k.tree.plot, .k.tree.plot$stratum)

    .agregate <- function(data){

      .out <- data.frame(n = tapply(data$k, data$k, length),
                         k = tapply(data$k, data$k, mean, na.rm = TRUE),
                         N = tapply(data$N, data$k, mean, na.rm = TRUE),
                         N.sd = tapply(data$N, data$k, stats::sd, na.rm = TRUE),
                         G = tapply(data$G, data$k, mean, na.rm = TRUE),
                         G.sd = tapply(data$G, data$k, stats::sd, na.rm = TRUE))

      .out$stratum <- data$stratum[1]

      return(.out)

    }

    .k.tree.plot <- do.call(rbind, lapply(.k.tree.plot.i, .agregate))


    .angle.count.plot.i <- split(.angle.count.plot, .angle.count.plot$stratum)

    .agregate <- function(data){

      .out <- data.frame(n = tapply(data$BAF, data$BAF, length),
                         BAF = tapply(data$BAF, data$BAF, mean, na.rm = TRUE),
                         N = tapply(data$N, data$BAF, mean, na.rm = TRUE),
                         N.sd = tapply(data$N, data$BAF, stats::sd,
                                       na.rm = TRUE),
                         G = tapply(data$G, data$BAF, mean, na.rm = TRUE),
                         G.sd = tapply(data$G, data$BAF, stats::sd,
                                       na.rm = TRUE))

      .out$stratum <- data$stratum[1]

      return(.out)

    }

    .angle.count.plot <- do.call(rbind, lapply(.angle.count.plot.i,
                                               .agregate))

  }

  # List with plots simulatioins and axis information

  .plots <- list(.fixed.area.plot = list(.fixed.area.plot,
                                         .fixed.area.plot.axis),
                 .k.tree.plot = list(.k.tree.plot, .k.tree.plot.axis),
                 .angle.count.plot = list(.angle.count.plot,
                                          .angle.count.plot.axis))


  ##### PLOTS ####

  # Individual plots ----


  if (!average & !all.plot.designs) {

    # Defining margins

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    graphics::par(mfrow = c(3, 2))
    graphics::par(mar = c(4, 5, 2, 1))



    for (i in c(1:length(.plots))) {

      .data <- .plots[[i]]

      .axis <- .data[[2]]
      .data <- .data[[1]]

      # Density (tree/ha)

      plot(.data[, 3], .data$N, type = "n",
           main = .axis[[1]],
           xlab = .axis[[2]],
           ylab = expression("N.tls" ~ (trees/ha)),
           xlim = c(min(.data[, 3]), max(.data[, 3])),
           ylim = c(min(.data$N), max(.data$N)))

      # The rest of the plots

      for (.i in unique(.data$id)) {

        .dat <- subset(.data, .data$id == .i)

        .x <- .dat[, 3]
        .y <- .dat$N

        graphics::lines(.x, .y, col = scales::alpha(.dat$stratum, 0.5))

      }

      graphics::legend("topright",
                       legend = unique(as.character(.data$stratum)),
                       text.col = unique(.data$stratum),
                       bty = "n", pch = 15,
                       col = unique(.data$stratum))

      # Basal area (m2/ha)

      plot(.data[, 3], .data$G, type = "n",
           main = .axis[[1]],
           xlab = .axis[[2]],
           ylab = expression("G.tls" ~ (m^{2}/ha)),
           xlim = c(min(.data[, 3]), max(.data[, 3])),
           ylim = c(min(.data$G), max(.data$G)))

      # The rest of the plots

      for (.i in unique(.data$id)) {

        .dat <- subset(.data, .data$id == .i)

        .x <- .dat[, 3]
        .y <- .dat$G

        graphics::lines(.x, .y, col = scales::alpha(.dat$stratum, 0.5))

      }

      graphics::legend("topright",
                       legend = unique(as.character(.data$stratum)),
                       text.col = unique(.data$stratum),
                       bty = "n", pch = 15,
                       col = unique(.data$stratum))

    }
  }

  # Mean values ----

  if (average) {

    # Defining margins

    graphics::par(mfrow = c(3, 2))
    graphics::par(mar = c(4, 5, 2, 4))

    for (i in c(1:length(.plots))) {

      .data <- .plots[[i]]

      .axis <- .data[[2]]
      .data <- .data[[1]]

      .data$N.min <- .data$N - (.data$N.sd / 2)
      .data$N.max <- .data$N + (.data$N.sd / 2)

      .data$G.min <- .data$G - (.data$G.sd / 2)
      .data$G.max <- .data$G + (.data$G.sd / 2)


      plot(.data[, 2], .data$N, type = "n",
           main = .axis[[1]],
           xlab = .axis[[2]],
           ylab = expression("N.tls" ~ (trees/ha)),
           xlim = c(min(.data[, 2]), max(.data[, 2])),
           ylim = c(min(.data[, c("N", "N.min", "N.max")], na.rm = TRUE),
                    max(.data[, c("N", "N.min", "N.max")], na.rm = TRUE)))


      for (i in unique(.data$stratum)) {

        .dat <- subset(.data, .data$stratum == i)

        .area <- .dat[which(!is.na(.dat$N.sd)), , drop = FALSE]

        graphics::par(new=T)

        plot(.data[, 2], .data$N, type = "n",
             xlab = "", ylab = "",
             xlim = c(min(.data[, 2]), max(.data[, 2])),
             ylim = c(min(.data[, c("N", "N.min", "N.max")], na.rm = TRUE),
                      max(.data[, c("N", "N.min", "N.max")], na.rm = TRUE)))

        graphics::polygon(c(.area[, 2], rev(.area[, 2])),
                          c(.area$N.max, rev(.area$N.min)),
                          col = scales::alpha(.dat$stratum, 0.5), border = NA)

        graphics::lines(.dat[, 2], .dat$N, col = .dat$stratum)

        graphics::par(new=T)

        plot(.dat[, 2], .dat$n,
             ylim = c(0, max(.data$n)),
             xaxt="n", yaxt="n",
             ylab = "", xlab = "",
             type="l", lty = 3, col = .dat$stratum)

        graphics::axis(4)

        graphics::mtext(side = 4, line = 3, "No plots", col = "black",
                        cex = 0.75)

      }

      graphics::legend("topright",
                       legend = unique(as.character(.data$stratum)),
                       text.col = unique(.data$stratum),
                       bty = "n", pch = 15,
                       col = unique(.data$stratum))

      # Basal area (m2/ha)

      plot(.data[, 2], .data$G, type = "n",
           main = .axis[[1]],
           xlab = .axis[[2]],
           ylab = expression("G.tls" ~ (m^{2}/ha)),
           xlim = c(min(.data[, 2]), max(.data[, 2])),
           ylim = c(min(.data[, c("G", "G.min", "G.max")], na.rm = TRUE),
                    max(.data[, c("G", "G.min", "G.max")], na.rm = TRUE)))


      for (i in unique(.data$stratum)) {

        .dat <- subset(.data, .data$stratum == i)

        .area <- .dat[which(!is.na(.dat$G.sd)), , drop = FALSE]

        graphics::par(new=T)

        plot(.data[, 2], .data$G, type = "n",
             xlab = "", ylab = "",
             xlim = c(min(.data[, 2]), max(.data[, 2])),
             ylim = c(min(.data[, c("G", "G.min", "G.max")], na.rm = TRUE),
                      max(.data[, c("G", "G.min", "G.max")], na.rm = TRUE)))


        graphics::polygon(c(.area[, 2], rev(.area[, 2])),
                          c(.area$G.max, rev(.area$G.min)),
                          col = scales::alpha(.dat$stratum, 0.5), border = NA)

        graphics::lines(.dat[, 2], .dat$G, col = .dat$stratum)

        graphics::par(new=T)

        plot(.dat[, 2], .dat$n,
             ylim = c(0, max(.data$n)),
             xaxt="n", yaxt="n",
             ylab = "", xlab = "",
             type="l", lty = 3, col = .dat$stratum)

        graphics::axis(4)

        graphics::mtext(side = 4, line = 3, "No plots", col = "black",
                        cex = 0.75)

      }

      graphics::legend("topright",
                       legend = unique(as.character(.data$stratum)),
                       text.col = unique(.data$stratum),
                       bty = "n", pch = 15,
                       col = unique(.data$stratum))

    }

  }

  # All plot design together ----

  if (all.plot.designs) {

    if (average) grDevices::dev.new()

    # Defining margins

    graphics::par(mfrow = c(1, 2))
    graphics::par(mar = c(7, 5, 2, 1))

    .N.max <- max(.plots$.fixed.area.plot[[1]]$N,
                  .plots$.k.tree.plot[[1]]$N,
                  .plots$.angle.count.plot[[1]]$N)

    .N.min <- min(.plots$.fixed.area.plot[[1]]$N,
                  .plots$.k.tree.plot[[1]]$N,
                  .plots$.angle.count.plot[[1]]$N)



    for (i in c(1:length(.plots))) {

      .data <- .plots[[i]]

      .axis <- .data[[2]]
      .data <- .data[[1]]

      .data$N.min <- .data$N - (.data$N.sd / 2)
      .data$N.max <- .data$N + (.data$N.sd / 2)

      if (i > 1)
        graphics::par(new=T)

      plot(.data[, 2], .data$N, type = "n",
           main = "Density",
           xlab = "",
           xaxt = "n",
           ylab = expression("N.tls" ~ (trees/ha)),
           xlim = c(0, max(.data[, 2])),
           ylim = c(.N.min, .N.max))

      graphics::par(new = T)

      if(i == 1)
        j <- 1

      if(i == 2)
        j <- 3

      if(i == 3)
        j <- 5

      graphics::axis(1, line = j, col = i, col.ticks = i, col.axis = i)

      graphics::mtext(.axis[[2]], side = 1, line = j, at = 0, adj = 1, col = i)

      .area <- .data[which(!is.na(.data$N.sd)), , drop = FALSE]

      graphics::polygon(c(.area[, 2], rev(.area[, 2])),
                        c(.area$N.max, rev(.area$N.min)),
                        col = scales::alpha(i, 0.25), border = NA)

      graphics::lines(.data[, 2], .data$N, col = scales::alpha(i, 0.75),
                      lwd = 2)

    }


    graphics::legend("topright",
                     legend = c(.plots$.fixed.area.plot[[2]]$main,
                                .plots$.k.tree.plot[[2]]$main,
                                .plots$.angle.count.plot[[2]]$main),
                     text.col = c(1:length(.plots)), cex = 0.75,
                     bty = "n", pch = 15,
                     col = c(1:length(.plots)))


    # Basal area (m2/ha)

    .G.max <- max(.plots$.fixed.area.plot[[1]]$G,
                  .plots$.k.tree.plot[[1]]$G,
                  .plots$.angle.count.plot[[1]]$G)

    .G.min <- min(.plots$.fixed.area.plot[[1]]$G,
                  .plots$.k.tree.plot[[1]]$G,
                  .plots$.angle.count.plot[[1]]$G)


    for (i in c(1:length(.plots))) {

      .data <- .plots[[i]]

      .axis <- .data[[2]]
      .data <- .data[[1]]

      .data$G.min <- .data$G - (.data$G.sd / 2)
      .data$G.max <- .data$G + (.data$G.sd / 2)

      if (i > 1)
        graphics::par(new=T)

      plot(.data[, 2], .data$G, type = "n",
           main = "Basal Area",
           xlab = "",
           xaxt = "n",
           ylab = expression("G.tls" ~ (m^{2}/ha)),
           xlim = c(0, max(.data[, 2])),
           ylim = c(.G.min, .G.max))

      graphics::par(new=T)

      if(i == 1)
        j <- 1

      if(i == 2)
        j <- 3

      if(i == 3)
        j <- 5

      graphics::axis(1, line = j, col = i, col.ticks = i, col.axis = i)

      graphics::mtext(.axis[[2]], side = 1, line = j, adj = 1, at = 0, col = i)

      .area <- .data[which(!is.na(.data$G.sd)), , drop = FALSE]

      graphics::polygon(c(.area[, 2], rev(.area[, 2])),
                        c(.area$G.max, rev(.area$G.min)),
                        col = scales::alpha(i, 0.25), border = NA)

      graphics::lines(.data[, 2], .data$G, col = scales::alpha(i, 0.75),
                      lwd = 2)

    }

    graphics::legend("topright",
                     legend = c(.plots$.fixed.area.plot[[2]]$main,
                                .plots$.k.tree.plot[[2]]$main,
                                .plots$.angle.count.plot[[2]]$main),
                     text.col = c(1:length(.plots)), cex = 0.75,
                     bty = "n", pch = 15,
                     col = c(1:length(.plots)))

  }

  return(invisible())

}
