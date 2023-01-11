
estimation.plot.size <- function(tree.tls,
                                 plot.parameters = data.frame(radius.max = 25,
                                                              k.max = 50,
                                                              BAF.max = 4),
                                 dbh.min = 4, average = FALSE,
                                 all.plot.designs = FALSE) {


  # Define stratum as 1 (by default) when tree.list$stratum is missed
  if(is.null(tree.tls$stratum) | all.plot.designs)
    tree.tls$stratum <- 1

  if(is.factor(tree.tls$stratum))
    tree.tls$stratum <- as.numeric(tree.tls$stratum)

  if(is.character(tree.tls$stratum))
    tree.tls$stratum <- as.numeric(as.factor(tree.tls$stratum))



  # Call internal function
  .est <- .sim.calc(funct = "est", tree.tls = tree.tls,
                    tree.ds = NULL, tree.field = NULL,
                    plot.design = c("fixed.area", "k.tree", "angle.count"),
                    plot.parameters = plot.parameters, scan.approach = "single",
                    var.metr = list(tls = c("N.tls", "G.tls")), v.calc = "parab",
                    dbh.min = dbh.min,  h.min = 1.3, max.dist = Inf,
                    dir.data = NULL, save.result = FALSE, dir.result = NULL)
  .fixed.area.plot <- .est$fixed.area
  .k.tree.plot <- .est$k.tree
  .angle.count.plot <- .est$angle.count

  # Create labels for plot axis
  .fixed.area.plot.axis <- list(main = "Circular fixed area plots", xlab = "Radius (m)")
  .k.tree.plot.axis <- list(main = "k-tree plots", xlab = "k-tree (trees)")
  .angle.count.plot.axis <- list(main = "Angle-count plots",
                                 xlab = expression("BAF" ~ (m^{2}/ha)))



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
