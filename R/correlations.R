
correlations <- function(simulations,
                         variables = c("N", "G", "V", "d", "dg", "d.0", "h",
                                       "h.0"),
                         method = c("pearson", "spearman"), save.result = TRUE,
                         dir.result = NULL) {

  # Checking if id columns are characters or factors, and converting to numeric in that cases
  if(!is.null(simulations$fixed.area) & is.character(simulations$fixed.area$id) | is.factor(simulations$fixed.area$id))
    simulations$fixed.area$id <- as.numeric(as.factor(simulations$fixed.area$id))

  if(!is.null(simulations$k.tree) & is.character(simulations$k.tree$id) | is.factor(simulations$k.tree$id))
    simulations$k.tree$id <- as.numeric(as.factor(simulations$k.tree$id))

  if(!is.null(simulations$angle.count) & is.character(simulations$angle.count$id) | is.factor(simulations$angle.count$id))
    simulations$angle.count$id <- as.numeric(as.factor(simulations$angle.count$id))


  # Define a character vector containing index name (radius, k or BAF) for each
  # available plot design
  .plot.design <- c(fixed.area = "radius", k.tree = "k", angle.count = "BAF")

  # Define character vectors containing the implemented field variables and TLS
  # metrics
  .field.names <- c(
    # Density (trees/ha), basal area (m2/ha) and volume (m3/ha)
    "N", "G", "V",

    # Volume (m3/ha) and biomass (Mg/ha) provided by user
    # if("W" %in% colnames(simulations$fixed.area))
    "V.user", "W.user",

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
    "n.pts", "n.pts.est", "n.pts.red", "n.pts.red.est",

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

    # Height percentiles (m)
    sprintf("P%02i", c(1, 5, 10, 20, 25, 30, 40, 50, 60, 70, 75,
                       80, 90, 95, 99)),

    # Points metrics
    # Z coordinate
    "mean.z", "mean.q.z", "mean.g.z", "mean.h.z", "median.z", "mode.z",
    "max.z", "min.z", "var.z", "sd.z", "CV.z", "D.z", "ID.z",
    "kurtosis.z", "skewness.z",
    "p.a.mean.z", "p.a.mode.z", "p.a.2m.z",
    "p.b.mean.z", "p.b.mode.z", "p.b.2m.z", "CRR.z",
    "L2.z", "L3.z", "L4.z", "L3.mu.z", "L4.mu.z",
    "L.CV.z",
    "median.a.d.z", "mode.a.d.z",
    "weibull_c.z", "weibull_b.z",

    # Rho coordinate
    "mean.rho", "mean.q.rho", "mean.g.rho", "mean.h.rho", "median.rho", "mode.rho",
    "max.rho", "min.rho", "var.rho", "sd.rho", "CV.rho", "D.rho", "ID.rho",
    "kurtosis.rho", "skewness.rho",
    "p.a.mean.rho", "p.a.mode.rho",
    "p.b.mean.rho", "p.b.mode.rho", "CRR.rho",
    "L2.rho", "L3.rho", "L4.rho", "L3.mu.rho", "L4.mu.rho",
    "L.CV.rho",
    "median.a.d.rho", "mode.a.d.rho",
    "weibull_c.rho", "weibull_b.rho",

    # R coordinate
    "mean.r", "mean.q.r", "mean.g.r", "mean.h.r", "median.r", "mode.r",
    "max.r", "min.r", "var.r", "sd.r", "CV.r", "D.r", "ID.r",
    "kurtosis.r", "skewness.r",
    "p.a.mean.r", "p.a.mode.r",
    "p.b.mean.r", "p.b.mode.r", "CRR.r",
    "L2.r", "L3.r", "L4.r", "L3.mu.r", "L4.mu.r",
    "L.CV.r",
    "median.a.d.r", "mode.a.d.r",
    "weibull_c.r", "weibull_b.r")

  # Define a character vector containing the available correlation measurements
  .cor.method <- c("pearson", "spearman")


  # Check arguments


  # 'simulations' must be a list with all or any of the preset elements, and
  # being at least one of them no NULL
  if (!is.list(simulations)) stop("'simulations' must be a list")
  if (is.null(simulations) || all(!names(.plot.design) %in% names(simulations)))
    stop("'simulations' must have at least one of the following elements:",
         "'fixed.area', 'k.tree' or 'angle.count'")
  if (any(!names(simulations) %in% names(.plot.design))) {

    simulations <- simulations[names(simulations) %in% names(.plot.design)]
    warning("There is any element in 'simulations' which do not correspond ",
            "with preset ones. It was not taken into account during the ",
            "execution")

  }
  simulations <- simulations[!sapply(simulations, is.null)]
  if (length(simulations) == 0)
    stop("'simulations' must have at least one of the following elements ",
         "different from 'NULL': 'fixed.area', 'k.tree' or 'angle.count'")
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

  # 'method' must be c("pearson", "spearman") (by default) or a character vector
  # with all or any of the preset methods
  if (!is.character(method) || !is.vector(method))
    stop("'method' must be a character vector")
  if (length(method) == 0 || all(!method %in% .cor.method))
    stop("'method' must have at least one of the following character strings: ",
         paste("'", .cor.method, "'", sep = "", collapse = ", "))
  if (any(!method %in% .cor.method)) {

    method <- method[method %in% .cor.method]
    warning("There is any character string in 'method' which do not ",
            "correspond with preset ones. It was not taken into account ",
            "during the execution")

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

  # Define a character vector containing color palette for correlations plots,
  # if 'save.result' is TRUE
  if (save.result) {

    .tls.color <- substr(.tls.names, 1, 1)
    names(.tls.color) <- .tls.names
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
    .tls.color[.tls.color == "P"] <-
      grDevices::gray.colors(sum(.tls.color == "P"))

    # Means
    .tls.color[.tls.color == "m"] <-
      grDevices::hcl.colors(sum(.tls.color == "m"), "Reds2")

    # Percentages
    .tls.color[.tls.color == "p"] <-
      grDevices::hcl.colors(sum(.tls.color == "p"), "Greens2")

    # Weibull
    .tls.color[.tls.color == "w"] <-
      grDevices::hcl.colors(sum(.tls.color == "w"), "Blues2")


  }

  # Define character vectors containing names for correlation values and their
  # associated p-values, and names for optimal correlation values and their
  # corresponding TLS metric
  .cor.names <- c("estimate", "p.value")
  .opt.cor.names <- c("cor", "metric")

  # Define empty lists where correlations, p-values, and optimal correlations
  # and corresponding metrics will be saved
  correlations <- vector("list", length(method))
  names(correlations) <- method
  correlations.pval <-correlations
  opt.correlations <- correlations
  for (.i in names(correlations)) {

    correlations[[.i]] <- vector("list", length(simulations))
    names(correlations[[.i]]) <- names(simulations)
    correlations.pval[[.i]] <- correlations[[.i]]
    opt.correlations[[.i]] <- correlations[[.i]]

  }


  # Loop for each plot design in 'simulations' argument
  for (.i in names(simulations)) {

    # Define initial time, and print message
    t0 <- Sys.time()
    message("Computing correlations for ",
            switch(.i, fixed.area = "fixed area ", k.tree = "k-tree ",
                   angle.count = "angle-count "), "plots")


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


    # Compute correlations, and optimal correlations ----


    # Create a character matrix containing all the pairs of field variable and
    # TLS metric to be considered for correlation calculations
    .cor.pairs <- .tls.names[.tls.names %in% dimnames(.data)[[2]]]
    .cor.pairs <- cbind(field = rep(.field.names, each = length(.cor.pairs)),
                        tls = rep(.cor.pairs, length(.field.names)))
    rownames(.cor.pairs) <- apply(.cor.pairs, 1, paste, collapse = ".")


    # Loop for each method
    for (.j in method){


      # Compute correlations and associated p-values for each radius, k or BAF


      # Compute correlations
      .cor.data <-
        sapply(dimnames(.data)[[1]],
               function(index, data, cor.pairs, method, cor.names){
                 # Compute correlations for each pair
                 cor.data <-
                   sapply(1:nrow(cor.pairs),
                          function(nr, data, cor.pairs, method, cor.names) {
                            # Select simulated data
                            data <- matrix(data[cor.pairs[nr, ], ],
                                           nrow = ncol(cor.pairs),
                                           ncol = dim(data)[2],
                                           dimnames = list(cor.pairs[nr,],
                                                           dimnames(data)[[2]]))
                            # Select plots with non NA values
                            data <- data[, apply(!is.na(data), 2, all),
                                         drop = FALSE]
                            # Compute correlations only if number of plots is
                            # greater than 2 and data have non zero standard
                            # deviation
                            if (ncol(data) > 2 &
                                all(apply(data, 1, stats::sd) != 0)) {
                              # Compute correlations for each method
                              cor.data <- get(paste(".cor", method, sep ="."))(
                                x = data[cor.pairs[nr, "field"], ],
                                y = data[cor.pairs[nr, "tls"], ])
                              cor.data <- unlist(cor.data[cor.names])
                            } else {
                              cor.data <- rep(NA, length(cor.names))
                              names(cor.data) <- cor.names
                            }
                            return(cor.data)
                          },
                          data = matrix(data[index, , ], nrow = dim(data)[[2]],
                                        ncol = dim(data)[[3]],
                                        dimnames = dimnames(data)[2:3]),
                          cor.pairs = cor.pairs, method = method,
                          cor.names = cor.names)
                 return(cor.data)
               },
               data = .data, cor.pairs = .cor.pairs, method = .j,
               cor.names = .cor.names)

      # Rearrange correlations results into a numeric 3-dimensional array
      .cor.data <- array(.cor.data,
                         dim = c(length(.cor.names), nrow(.cor.pairs),
                                 dim(.data)[1]),
                         dimnames = list(.cor.names, rownames(.cor.pairs),
                                         dimnames(.data)[[1]]))
      .cor.data <- aperm(.cor.data, 3:1)

      # Remove results corresponding to first radius, k or BAF values if all
      # correlations are NA values
      .any.nona <- which(apply(!is.na(.cor.data), 1, any))
      if (length(.any.nona) == 0) {

        warning("All the computed ",
                switch(.j, pearson = "Pearson correlations",
                       spearman = "Spearman correlations"),
                " for ",
                switch(.i, fixed.area = "fixed area",
                       k.tree = "k-tree",
                       angle.count = "angle-count"),
                " plot design case are 'NA' values")
      } else {

        .cor.data <- .cor.data[min(.any.nona):dim(.cor.data)[1], , ,
                               drop = FALSE]

      }

      # Save index and correlations values in results list, and write a csv file
      # containing them if 'save.result' is TRUE
      correlations[[.j]][[.i]] <-
        matrix(c(.index.val[dimnames(.cor.data)[[1]]],
                 .cor.data[, , .cor.names[1]]),
               nrow = dim(.cor.data)[1], ncol = 1 + dim(.cor.data)[2],
               dimnames = list(NULL,
                               c(.plot.design[.i], dimnames(.cor.data)[[2]])))
      if (save.result)
        utils::write.csv(correlations[[.j]][[.i]],
                         file = file.path(dir.result,
                                          paste("correlations", .i, .j, "csv",
                                                sep =".")),
                         row.names = FALSE)

      # Save index and correlations p-values in results lists
      correlations.pval[[.j]][[.i]] <-
        matrix(c(.index.val[dimnames(.cor.data)[[1]]],
                 .cor.data[, , .cor.names[2]]),
               nrow = dim(.cor.data)[1], ncol = 1 + dim(.cor.data)[2],
               dimnames = list(NULL,
                               c(.plot.design[.i], dimnames(.cor.data)[[2]])))


      # Compute optimal correlations for each radius, k or BAF


      # Create an empty data.frame where optimal correlations for each field
      # variable, and the name of the corresponding TLS metric, will be saved
      # for each radius, k or BAF
      .col.names <- as.vector(t(outer(.field.names, .opt.cor.names, paste,
                                      sep = ".")))
      .opt.cor.data <- matrix(NA, nrow = dim(.cor.data)[1],
                              ncol = length(.col.names),
                              dimnames = list(dimnames(.cor.data)[[1]],
                                              .col.names))
      .opt.cor.data <- data.frame(.opt.cor.data, stringsAsFactors = FALSE)


      # Loop for each field variable
      for (.k in .field.names) {

        # Select pairs corresponding to the field variable
        .col.names <- rownames(.cor.pairs)[.cor.pairs[, "field"] == .k]
        .cor.data.k <- matrix(.cor.data[, .col.names, "estimate"],
                              nrow = dim(.cor.data)[1],
                              ncol = length(.col.names),
                              dimnames = list(dimnames(.cor.data)[[1]],
                                              .col.names))

        # Remove results corresponding to radius, k or BAF values such that all
        # correlations are NA values
        .cor.data.k <- .cor.data.k[apply(!is.na(.cor.data.k), 1, any), ,
                                   drop = FALSE]

        # Obtain and save optimal correlations and corresponding TLS metric
        if (nrow(.cor.data.k) > 0) {

          # Obtain optimal correlations and the corresponding TLS metric
          .opt.cor.data.k <-
            lapply(1:nrow(.cor.data.k),
                   function(nr, data, names) {
                     pos <- colnames(data)[which.max(abs(data[nr, ]))]
                     val <- data[nr, pos]
                     opt.cor <- data.frame(val = val, pos = pos,
                                           row.names = rownames(data)[nr],
                                           stringsAsFactors = FALSE)
                     return(opt.cor)
                   },
                   data = .cor.data.k, names = .opt.cor.names)
          .opt.cor.data.k <- do.call(rbind, .opt.cor.data.k)
          .opt.cor.data.k[, "pos"] <- .cor.pairs[.opt.cor.data.k[, "pos"],
                                                 "tls"]
          colnames(.opt.cor.data.k) <- paste(.k, .opt.cor.names, sep = ".")

          # Save optimal correlations and the corresponding TLS metric
          .opt.cor.data[rownames(.opt.cor.data.k), colnames(.opt.cor.data.k)] <-
            .opt.cor.data.k

        }


        # Plot correlations related to the field variable, if 'save.result' is
        # TRUE


        if (save.result) {

          # Convert correlations matrix to data.frame, and add index values
          .cor.data.k <- data.frame(.index.val[rownames(.cor.data.k)],
                                    .cor.data.k, stringsAsFactors = FALSE)
          colnames(.cor.data.k)[1] <- .plot.design[.i]

          # Define title, subtitle, and axis names, and create empty plot
          .title <-
            switch(.k, N = "Density (N, trees/ha)",
                   G = "Basal area (G, m<sup>2</sup>/ha)",
                   V = "Volume (V, m<sup>3</sup>/ha)",
                   V.user = "Volume (V, m<sup>3</sup>/ha) provided by user",
                   W.user = "Biomass (W, Mg/ha) provided by user",
                   d = "Arithmetic mean diameter (d, cm)",
                   dg = "Quadratic mean diameter (dg, cm)",
                   dgeom = "Geometric mean diameter (dgeom, cm)",
                   dharm = "Harmonic mean diameter (dharm, cm)",
                   d.0 = "Arithmetic mean dominant diameter (d.0, cm)",
                   dg.0 = "Quadratic mean dominant diameter (dg.0, cm)",
                   dgeom.0 = "Geometric mean dominant diameter (dgeom.0, cm)",
                   dharm.0 = "Harmonic mean dominant diameter (dharm.0, cm)",
                   h = "Arithmetic mean height (h, m)",
                   hg = "Quadratic mean height (hg, m)",
                   hgeom = "Geometric mean height (hgeom, m)",
                   hharm = "Harmonic mean height (hharm, m)",
                   h.0 = "Arithmetic mean dominant height (h.0, m)",
                   hg.0 = "Quadratic mean dominant height (hg.0, m)",
                   hgeom.0 = "Geometric mean dominant height (hgeom.0, m)",
                   hharm.0 = "Harmonic mean dominant height (hharm.0, m)")
          .subtitle <- paste("<br> <span style='font-size: 20px;'>",
                             switch(.i, fixed.area = "Circular fixed area plot",
                                    k.tree = "K-tree plot",
                                    angle.count = "Angle-count plot"),
                             "</span>", sep ="")
          .xaxis <- switch(.i, fixed.area = "Radius (m)",
                           k.tree = "K-tree (trees)",
                           angle.count = "BAF (m<sup>2</sup>/ha)")
          .yaxis <- switch(.j, pearson = "Pearson correlation",
                           spearman = "Spearman correlation")
          fig <-
            plotly::plot_ly(.cor.data.k, type = 'scatter', mode = 'lines') %>%
            plotly::layout(title = paste(.title, .subtitle, sep = ""), font = list(size = 25),
                           xaxis = list(title = .xaxis),
                           yaxis = list (title = .yaxis),
                           margin = list(t = 100))

          # Add traces
          for (.l in colnames(.cor.data.k)[colnames(.cor.data.k) !=
                                           .plot.design[.i]]) {

            .color <- .tls.color[.cor.pairs[.l, "tls"]]
            fig <- fig %>%
              plotly::add_trace(x = .cor.data.k[, .plot.design[.i]],
                                y = .cor.data.k[, .l],
                                name = .cor.pairs[.l, "tls"],
                                line = list(color = .color, width = 2,
                                            dash = 'dot'))

          }

          # Save correlations plot related to the field variable
          suppressWarnings(
            htmlwidgets::saveWidget(widget = fig,
                                    file = file.path(dir.result,
                                                     paste("correlations", .k,
                                                           .i, .j, "html",
                                                           sep = ".")),
                                    selfcontained = TRUE,
                                    libdir = file.path(dir.result,
                                                       "correlations_files")))

        }

      }

      # Save index, optimal correlations, and names of the corresponding TLS
      # metrics in results list, and write csv file containing them
      opt.correlations[[.j]][[.i]] <- cbind(.index.val[rownames(.opt.cor.data)],
                                            .opt.cor.data)
      rownames(opt.correlations[[.j]][[.i]]) <- NULL
      colnames(opt.correlations[[.j]][[.i]])[1] <- .plot.design[.i]
      if (save.result)
        utils::write.csv(opt.correlations[[.j]][[.i]],
                         file = file.path(dir.result, paste("opt.correlations",
                                                            .i, .j, "csv",
                                                            sep = ".")),
                         row.names = FALSE)

    }

    # Define final time, and print message
    t1 <- Sys.time()
    message(" (", format(round(difftime(t1, t0, units = "secs"), 2)), ")")

  }


  return(list(correlations = correlations,
              correlations.pval =  correlations.pval,
              opt.correlations = opt.correlations))

}
