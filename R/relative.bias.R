
relative.bias <- function(simulations,
                          variables = c("N", "G", "V", "d", "dg", "d.0", "h",
                                        "h.0"),
                          save.result = TRUE, dir.result = NULL) {


  # Define a character vector containing index name (radius, k or BAF) for each
  # available plot design
  .plot.design <- c(fixed.area.plot = "radius", k.tree.plot = "k",
                    angle.count.plot = "BAF")

  # Define character vectors containing the implemented field variables and TLS
  # metrics
  .field.names <- c(
                    # Density (trees/ha), basal area (m2/ha) and volume (m3/ha)
                    "N", "G", "V",

                    # Mean diameters (cm), and mean dominant diameters (cm)
                    "d", "dg", "dgeom", "dharm",
                    paste(c("d", "dg", "dgeom", "dharm"), "0", sep = "."),

                    # Mean heights (m), and mean dominant heights (m)
                    "h", "hg", "hgeom", "hharm",
                    paste(c("h", "hg", "hgeom", "hharm"), "0", sep = "."))
  .tls.names <- c(
                  # Density (trees/ha)
                  "N.tls", "N.hn", "N.hr", "N.hn.cov", "N.hr.cov", "N.sh",
                  "N.pam",

                  # Number of points
                  "num.points", "num.points.est", "num.points.hom",
                  "num.points.hom.est",

                  # Basal area (m2/ha)
                  "G.tls", "G.hn", "G.hr", "G.hn.cov", "G.hr.cov", "G.sh",
                  "G.pam",

                  # Volume (m3/ha)
                  "V.tls", "V.hn", "V.hr", "V.hn.cov", "V.hr.cov", "V.sh",
                  "V.pam",

                  # Mean diameters (cm), and mean dominant diameters (cm)
                  paste(c("d", "dg", "dgeom", "dharm"), "tls", sep = "."),
                  paste(c("d", "dg", "dgeom", "dharm"), "0.tls", sep = "."),

                  # Mean heights (m), and mean dominant heights (m)
                  paste(c("h", "hg", "hgeom", "hharm"), "tls", sep = "."),
                  paste(c("h", "hg", "hgeom", "hharm"), "0.tls", sep = "."),

                  # Height percentiles P99 (m)
                  # Remark: only for relative bias regarding heights and
                  # dominant heights
                  sprintf("P%02i", 99))


  # Check arguments


  # 'simulations' must be a list with all or any of the preset elements, and
  # being at least one of them no NULL
  if (!is.list(simulations)) stop("'simulations' must be a list")
  if (is.null(simulations) || all(!names(.plot.design) %in% names(simulations)))
    stop("'simulations' must have at least one of the following elements:",
         "'fixed.area.plot', 'k.tree.plot' or 'angle.count.plot'")
  if (any(!names(simulations) %in% names(.plot.design))) {

    simulations <- simulations[names(simulations) %in% names(.plot.design)]
    warning("There is any element in 'simulations' which do not correspond ",
            "with preset ones. It was not taken into account during the ",
            "execution")

  }
  simulations <- simulations[!sapply(simulations, is.null)]
  if (length(simulations) == 0)
    stop("'simulations' must have at least one of the following elements ",
         "different from 'NULL': 'fixed.area.plot', 'k.tree.plot' or ",
         "'angle.count.plot'")
  for (.i in names(simulations)) {

    # All elements in 'simulations' must be data frames with at least one row,
    # certain mandatory columns and most of them numeric
    if (!is.data.frame(simulations[[.i]]))
      stop("'simulations$", .i, "' must be a data.frame")
    if (nrow(simulations[[.i]]) < 1)
      stop("'simulations$", .i, "' must have at least one row")
    if (!"id" %in% colnames(simulations[[.i]]))
      stop("'simulations$", .i, "' must have a column named 'id'")
    if (!.plot.design[.i] %in% colnames(simulations[[.i]]))
      stop("'simulations$", .i, "' must have a column named '",
           .plot.design[.i], "'")
    if (all(!.field.names %in% colnames(simulations[[.i]])))
      stop("'simulations$", .i, "' must have at least one field estimations ",
           "column")
    if (all(!.tls.names %in% colnames(simulations[[.i]])))
      stop("'simulations$", .i, "' must have at least one TLS metrics column")
    for (.j in colnames(simulations[[.i]])) {

      if (.j != "id" && !is.numeric(simulations[[.i]][, .j]))
        stop("Column '", .j,"' of 'simulations$", .i, "' must be numeric")

    }

  }

  # 'variables' must be c("N", "G", "V", "d", "dg", "d.0", "h", "h.0") (by
  # default) or a character vector with all or any of the preset estimations of
  # variables based on field data. Besides, simulations elements must have at
  # least the columns corresponding to 'variables'
  if (!is.character(variables) || !is.vector(variables))
    stop("'variables' must be a character vector")
  if (length(variables) == 0 || all(!variables %in% .field.names))
    stop("'variables' must have at least one of the following character ",
         "strings: ", paste("'", .field.names, "'", sep = "", collapse = ", "))
  if (any(!variables %in% .field.names)) {

    variables <- variables[variables %in% .field.names]
    warning("There is any character string in 'variables' which do not ",
            "correspond with preset ones. It was not taken into account ",
            "during the execution")

  }
  for (.i in names(simulations)) {

    if (any(!variables %in% colnames(simulations[[.i]])))
      stop("Any column corresponding to 'variables' is missing in ",
           "'simulations$", .i, "'")

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


  # Define values for certain auxiliary objects, and create empty lists where
  # results will be saved


  # Restrict field variables to them specified in 'variables' argument
  .field.names <- .field.names[.field.names %in% variables]

  # Define a character vector containing color palette for relative bias plots,
  # if 'save.result' is TRUE
  if (save.result) {

    .tls.color <- substr(.tls.names[substr(.tls.names, 1, 1) != "P"], 1, 1)
    names(.tls.color) <- .tls.names[substr(.tls.names, 1, 1) != "P"]
    # Density
    .tls.color[.tls.color == "N"] <-
      grDevices::hcl.colors(sum(.tls.color == "N") + 1,
                            "Reds")[-(sum(.tls.color == "N") + 1)]
    # Number of points
    .tls.color[.tls.color == "n"] <-
      grDevices::hcl.colors(sum(.tls.color == "n") + 1,
                            "Greens")[-(sum(.tls.color == "n") + 1)]
    # Basal area
    .tls.color[.tls.color == "G"] <-
      grDevices::hcl.colors(sum(.tls.color == "G") + 1,
                            "Blues")[-(sum(.tls.color == "G") + 1)]
    # Volume
    .tls.color[.tls.color == "V"] <-
      grDevices::hcl.colors(sum(.tls.color == "V") + 1,
                            "Oranges")[-(sum(.tls.color == "V") + 1)]
    # Diameters
    .tls.color[.tls.color == "d"] <-
      grDevices::hcl.colors(sum(.tls.color == "d") + 1,
                            "Purples")[-(sum(.tls.color == "d") + 1)]
    # Heights
    .tls.color[.tls.color == "h"] <-
      grDevices::hcl.colors(sum(.tls.color == "h"), "Magenta")
    # Percentiles
    .perc.names <- outer(.field.names[substr(.field.names, 1, 1) == "h"],
                         .tls.names[substr(.tls.names, 1, 1) == "P"], paste,
                         sep = ".")
    .perc.names <- as.vector(.perc.names)
    .perc.color <-  grDevices::gray.colors(length(.perc.names))
    names(.perc.color) <- .perc.names
    .tls.color <- c(.tls.color, .perc.color)

  }

  # Define an empty list where relative bias will be saved
  RB <- vector("list", length(simulations))
  names(RB) <- names(simulations)


  # Loop for each plot design in 'simulations' argument
  for (.i in names(simulations)) {

    # Define initial time, and print message
    t0 <- Sys.time()
    message("Computing relative bias for ",
        switch(.i, fixed.area.plot = "fixed area ", k.tree.plot = "k-tree ",
               angle.count.plot = "angle-count "), "plots")


    # Rearrange simulated data ----


    # Convert data.frame containing simulated data into a matrix
    simulations[[.i]] <- as.matrix(simulations[[.i]])

    # Create an empty numeric 3-dimensional array where simulated data will be
    # rearranged
    .index.val <- sort(unique(simulations[[.i]][, .plot.design[.i]]))
    .index.dec <- .decimals(.index.val)
    names(.index.val) <- .format.numb(x = .index.val, dec = .index.dec)
    .col.names <- c(.field.names,
                    .tls.names[.tls.names %in% colnames(simulations[[.i]])])
    .id.val <- sort(unique(simulations[[.i]][, "id"]))
    names(.id.val) <- as.character(.id.val)
    .data <- array(numeric(), dim = c(length(.index.val), length(.col.names),
                                      length(.id.val)),
                   dimnames = list(names(.index.val), .col.names,
                                   names(.id.val)))

    # Rearrange simulated data into the previously created array
    for (.j in names(.id.val)) {

      .row.names <- simulations[[.i]][simulations[[.i]][, "id"] == .id.val[.j],
                                      .plot.design[.i]]
      .row.names <- .format.numb(x = .row.names, dec = .index.dec)
      .data[.row.names, , .j] <-
        simulations[[.i]][simulations[[.i]][, "id"] == .id.val[.j],
                          dimnames(.data)[[2]]]

    }


    # Compute relative bias ----


    # Create a character matrix containing all the pairs of field variable and
    # TLS metric to be considered for relative bias calculations
    .RB.pairs <- .tls.names[.tls.names %in% dimnames(.data)[[2]]]
    .RB.pairs <- cbind(field = rep(.field.names, each = length(.RB.pairs)),
                        tls = rep(.RB.pairs, length(.field.names)))
    rownames(.RB.pairs) <- apply(.RB.pairs, 1, paste, collapse = ".")
    # Restrict possible pairs to comparable ones: field and corresponding TLS
    # counterpart; or field heights and TLS percentiles
    .row.names <-
      (substr(.RB.pairs[, "tls"], 1, nchar(.RB.pairs[, "field"])) ==
         .RB.pairs[, "field"]) |
      ((substr(.RB.pairs[, "field"], 1, 1) == "h" &
          substr(.RB.pairs[, "tls"], 1, 1) == "P"))
    .RB.pairs <- .RB.pairs[.row.names, , drop = FALSE]
    # Restrict pairs related to diameters and heights to those corresponding to
    # the same mean function: different from diameters and heights; or
    # diameters/heights corresponding to the same mean function; or field
    # heights and TLS percentiles
    .row.names <-
      (!substr(.RB.pairs[, "field"], 1, 1) %in% c("d", "h")) |
      (substr(.RB.pairs[, "field"], 1, 1) %in% c("d", "h") &
         paste(.RB.pairs[, "field"], "tls", sep = ".") == .RB.pairs[, "tls"]) |
      (substr(.RB.pairs[, "field"], 1, 1) == "h" &
         substr(.RB.pairs[, "tls"], 1, 1) == "P")
    .RB.pairs <- .RB.pairs[.row.names, , drop = FALSE]
    # Check if all field estimations have at least one TLS counterpart
    if (any(!.field.names %in% .RB.pairs[, "field"]))
      stop("'simulations$", .i, "' must have at least one TLS counterpart for ",
           "each element of 'variables' argument")


    # Compute relative bias for each radius, k or BAF


    # Compute relative bias
    .RB.data <-
      lapply(dimnames(.data)[[1]],
             function(index, data, RB.pairs){
               # Compute relative bias for each pair
               RB.data <-
                 sapply(1:nrow(RB.pairs),
                        function(nr, data, RB.pairs) {
                          # Select simulated data
                          data <- matrix(data[RB.pairs[nr, ], ],
                                         nrow = ncol(RB.pairs),
                                         ncol = dim(data)[2],
                                         dimnames = list(RB.pairs[nr,],
                                                         dimnames(data)[[2]]))
                          # Select plots with non NA values
                          data <- data[, apply(!is.na(data), 2, all),
                                       drop = FALSE]
                          # Compute relative bias only if number of plots is
                          # greater than 0
                          if (ncol(data) > 0) {
                            RB.data <-
                              ((mean(data[RB.pairs[nr, "tls"], ]) -
                                  mean(data[RB.pairs[nr, "field"], ])) /
                                 mean(data[RB.pairs[nr, "field"], ])) * 100
                          } else {
                            RB.data <- NA
                          }
                          return(RB.data)
                        },
                        data = matrix(data[index, , ], nrow = dim(data)[[2]],
                                      ncol = dim(data)[[3]],
                                      dimnames = dimnames(data)[2:3]),
                        RB.pairs = RB.pairs)
               return(RB.data)
             },
             data = .data, RB.pairs = .RB.pairs)
    .RB.data <- do.call(rbind, .RB.data)
    rownames(.RB.data) <- dimnames(.data)[[1]]
    colnames(.RB.data) <- rownames(.RB.pairs)

    # Remove results corresponding to first radius, k or BAF values if all
    # relative bias are NA values
    .any.nona <- which(apply(!is.na(.RB.data), 1, any))
    if (length(.any.nona) == 0) {

      warning("All the computed relative bias for ",
              switch(.i, fixed.area.plot = "fixed area", k.tree.plot = "k-tree",
                     angle.count.plot = "angle-count"),
              " plot design case are 'NA' values")

    } else {

      .RB.data <- .RB.data[min(.any.nona):dim(.RB.data)[1], , drop = FALSE]

    }

    # Save index and relative bias values in results list, and write csv file
    # containing them if 'save.result' is TRUE
    RB[[.i]] <- matrix(c(.index.val[rownames(.RB.data)], .RB.data),
                       nrow = nrow(.RB.data), ncol = 1 + ncol(.RB.data),
                       dimnames = list(NULL, c(.plot.design[.i],
                                               colnames(.RB.data))))
    if (save.result)
      utils::write.csv(RB[[.i]],
                       file = file.path(dir.result,
                                        paste("RB", .i, "csv", sep =".")),
                       row.names = FALSE)


    # Loop for plotting relative bias related to each field variable if
    # 'save.result' is TRUE
    # Remark: all field variables related to diameters are plotted in the same
    # interactive graphic, and the same for heights


    if (save.result) {

      .j.seq <- .field.names[!substr(.field.names, 1, 1) %in% c("d", "h")]
      if (any(substr(.field.names, 1, 1) == "d")) .j.seq <- c(.j.seq, "d")
      if (any(substr(.field.names, 1, 1) == "h")) .j.seq <- c(.j.seq, "h")

      for (.j in .j.seq) {

        # Select pairs corresponding to the field variable
        if (!.j %in% c("d", "h")) {

          .RB.data.j <-
            .RB.data[, rownames(.RB.pairs)[.RB.pairs[, "field"] == .j],
                     drop = FALSE]
          colnames(.RB.data.j) <- .RB.pairs[colnames(.RB.data.j), "tls"]
          .color <- .tls.color[colnames(.RB.data.j)]

        } else {

          .col.names <-
            rownames(.RB.pairs)[substr(.RB.pairs[, "field"], 1, 1) == .j &
                                  substr(.RB.pairs[, "tls"], 1, 1) != "P"]
          .RB.data.j <- .RB.data[, .col.names, drop = FALSE]
          .color <- .tls.color[.RB.pairs[.col.names, "tls"]]
          names(.color) <- colnames(.RB.data.j)

          # Add percentiles
          .col.names <-
            rownames(.RB.pairs)[substr(.RB.pairs[, "field"], 1, 1) == .j &
                                  substr(.RB.pairs[, "tls"], 1, 1) == "P"]
          if (length(.col.names) > 0) {

            .RB.data.j <- cbind(.RB.data.j,
                                .RB.data[, .col.names, drop = FALSE])
            .color <- c(.color, .tls.color[.col.names])

          }

        }

        # Remove results corresponding to radius, k or BAF values such that all
        # relative bias are NA values
        .RB.data.j <- .RB.data.j[apply(!is.na(.RB.data.j), 1, any), ,
                                 drop = FALSE]

        # Convert relative bias matrix to data.frame, and add index values
        .RB.data.j <- data.frame(.index.val[rownames(.RB.data.j)], .RB.data.j,
                                 stringsAsFactors = FALSE)
        colnames(.RB.data.j)[1] <- .plot.design[.i]

        # Define title, subtitle, and axis names, and create empty plot
        .title <- switch(.j, N = "Density (N, trees/ha)",
                         G = "Basal area (G, m<sup>2</sup>/ha)",
                         V = "Volume (V, m<sup>3</sup>/ha)",
                         d = "Mean diameters (cm)",
                         h = "Mean heights (m)")
        .subtitle <- paste("<br> <span style='font-size: 12px;'>",
                           switch(.i, fixed.area.plot = "Fixed area",
                                  k.tree.plot = "K-tree",
                                  angle.count.plot = "Angle-count"),
                           "</span>", sep ="")
        .xaxis <- switch(.i, fixed.area.plot = "Radius (m)",
                         k.tree.plot = "K-tree (trees)",
                         angle.count.plot = "BAF (m<sup>2</sup>/ha)")
        .yaxis <- "Relative bias (%)"
        fig <-
          plotly::plot_ly(.RB.data.j, type = 'scatter', mode = 'lines') %>%
          plotly::layout(title = paste(.title, .subtitle, sep = ""),
                         xaxis = list(title = .xaxis),
                         yaxis = list (title = .yaxis), margin = list(t = 50))

        # Add traces
        for (.k in
             colnames(.RB.data.j)[colnames(.RB.data.j) != .plot.design[.i]]) {

          fig <- fig %>%
            plotly::add_trace(x = .RB.data.j[, .plot.design[.i]],
                              y = .RB.data.j[, .k], name = .k,
                              line = list(color = .color[.k], width = 2,
                                          dash = 'dot'))

        }

        # Save relative bias plot related to the field variable
        suppressWarnings(
          htmlwidgets::saveWidget(widget = fig,
                                  file = file.path(dir.result,
                                                   paste("RB", .j, .i, "html",
                                                         sep = ".")),
                                  selfcontained = TRUE,
                                  libdir = file.path(dir.result, "RB_files")))

      }

    }

    # Define final time, and print message
    t1 <- Sys.time()
    message(" (", format(round(difftime(t1, t0, units = "secs"), 2)), ")")

  }


  return(RB)

}
