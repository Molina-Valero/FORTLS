
optimize.plot.design <- function(correlations,
                                 variables = c("N", "G", "V", "d", "dg", "d.0",
                                               "h", "h.0"),
                                 dir.result = NULL) {

  opt.cor <- correlations

  if("W" %in% colnames(opt.cor$pearson$fixed.area.plot) | "W" %in% colnames(opt.cor$spearman$fixed.area.plot))
    variables = c(variables, "W")

  # Define a character vector containing index name (radius, k or BAF) for each
  # available plot design
  .plot.design <- c(fixed.area.plot = "radius", k.tree.plot = "k",
                    angle.count.plot = "BAF")

  # Define character vectors containing the implemented field variables and TLS
  # metrics
  .field.names <- c(
                    # Density (trees/ha), basal area (m2/ha) and volume (m3/ha)
                    "N", "G", "V",

                    # Biomass (Mg/ha)
                    if("W" %in% colnames(opt.cor$pearson$fixed.area.plot) | "W" %in% colnames(opt.cor$spearman$fixed.area.plot))
                      "W",

                    # Mean diameters (cm), and mean dominant diameters (cm)
                    "d", "dg", "dgeom", "dharm",
                    paste(c("d", "dg", "dgeom", "dharm"), "0", sep = "."),

                    # Mean heights (m), and mean dominant heights (m)
                    "h", "hg", "hgeom", "hharm",
                    paste(c("h", "hg", "hgeom", "hharm"), "0", sep = "."))

  # Define a character vector containing colors for heatmaps
  .color <- c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#f7f7f7", "#fddbc7",
              "#ef8a62", "#b2182b")

  # Define a character vector containing the available correlation measurements
  .cor.method <- c("pearson", "spearman")

  # Define character matrix containing column names for optimal correlation
  # values and their corresponding TLS metric
  .opt.cor.names <- c("cor", "metric")
  .opt.cor.names <- matrix(outer(.field.names, .opt.cor.names, paste,
                                 sep = "."),
                           nrow = length(.field.names),
                           ncol = length(.opt.cor.names),
                           dimnames = list(.field.names, .opt.cor.names))


  # Check arguments


  # 'correlations' must be a list with all or any of the preset elements, and
  # being at least one of them no NULL
  if (!is.list(correlations)) stop("'correlations' must be a list")
  if (is.null(correlations) || all(!.cor.method %in% names(correlations)))
    stop("'correlations' must have at least one of the following elements:",
         "'pearson' or 'spearman'")
  if (any(!names(correlations) %in% .cor.method)) {

    correlations <- correlations[names(correlations) %in% .cor.method]
    warning("There is any element in 'correlations' which do not correspond ",
            "with preset ones. It was not taken into account during the ",
            "execution")

  }
  correlations <- correlations[!sapply(correlations, is.null)]
  if (length(correlations) == 0)
    stop("'correlations' must have at least one of the following elements ",
         "different from 'NULL': 'pearson' or 'spearman'")
  for (.i in names(correlations)) {

    # Each 'correlations' element must be a list with all or any of the preset
    # elements, and being at least one of them no NULL
    if (!is.list(correlations[[.i]]))
      stop("'correlations$", .i, "' must be a list")
    if (is.null(correlations[[.i]]) ||
        all(!names(.plot.design) %in% names(correlations[[.i]])))
      stop("'correlations$", .i, "' must have at least one of the following ",
           "elements: 'fixed.area.plot', 'k.tree.plot' or 'angle.count.plot'")
    if (any(!names(correlations[[.i]]) %in% names(.plot.design))) {

      correlations[[.i]] <- correlations[[.i]][names(correlations[[.i]]) %in%
                                                 names(.plot.design)]
      warning("There is any element in 'correlations$", .i, "' which do not ",
              "correspond with preset ones. It was not taken into account ",
              "during the execution")

    }
    correlations[[.i]] <- correlations[[.i]][!sapply(correlations[[.i]],
                                                     is.null)]
    if (length(correlations[[.i]]) == 0)
      stop("'correlations$", .i, "' must have at least one of the following ",
           "elements different from 'NULL': 'fixed.area.plot', 'k.tree.plot' ",
           "or 'angle.count.plot'")
    for (.j in names(correlations[[.i]])) {

      # All elements in each 'correlations' element must be data frames with at
      # least one row, certain mandatory columns and all of them numeric
      if (!is.data.frame(correlations[[.i]][[.j]]))
        stop("'correlations$", .i, "$", .j, "' must be a data.frame")
      if (nrow(correlations[[.i]][[.j]]) < 1)
        stop("'correlations$", .i, "$", .j, "' must have at least one row")
      if (!.plot.design[.j] %in% colnames(correlations[[.i]][[.j]]))
        stop("'correlations$", .i, "$", .j, "' must have a column named '",
             .plot.design[.j], "'")
      if (sum(apply(apply(.opt.cor.names, 1:2,
                          function(x) {x %in%
                              colnames(correlations[[.i]][[.j]])}),
                    1, all)) == 0)
        stop("'correlations$", .i, "$", .j, "' must have at least a ",
             "('<x>.cor', '<x>.metric') pair for the same field estimation")
      for (.k in colnames(correlations[[.i]][[.j]])) {

        if (!.k %in% .opt.cor.names[, "metric"] &&
            !is.numeric(correlations[[.i]][[.j]][, .k]))
          stop("Column '", .k,"' of 'correlations$", .i, "$", .j, "' must be ",
               "numeric")

      }

    }

  }

  # 'variables' must be c("N", "G", "V", "d", "dg", "d.0", "h", "h.0") (by
  # default) or a character vector with all or any of the preset estimations of
  # variables based on field data. Besides, 'correlations' elements must have at
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
  for (.i in names(correlations)) {

    for (.j in names(correlations[[.i]])) {

      if (any(!as.vector(.opt.cor.names[variables, ]) %in%
              colnames(correlations[[.i]][[.j]])))
        stop("Any column corresponding to 'variables' is missing in ",
             "'correlations$", .i, "$", .j, "'")

    }

  }

  # 'dir.result' must be NULL (by default) or a character string containing the
  # absolute path to a existing directory
  if (!is.null(dir.result)) {

    if (!is.character(dir.result))
      stop("'dir.result' must be a character string")
    if (length(dir.result) > 1) {

      dir.result <- dir.result[1]
      warning("Only first value in 'dir.result' was taken into account during ",
              "the execution")

    }
    if (!dir.exists(dir.result)) stop("'dir.result' directory must exist")

  } else {

    # Define working directory as directory by default for saving files
    dir.result <- getwd()

  }


  # Loop for each correlation measure
  for (.i in names(correlations)) {

    # Define initial time, and print message
    t0 <- Sys.time()
    message("Plotting heatmap(s) for ",
        switch(.i, pearson = "optimal Pearson's correlations ",
               spearman = "optimal Spearman's correlations "))


    # Loop for each plot design
    for (.j in names(correlations[[.i]])) {


      # Plot optimal correlations for each radius, k or BAF


      # Select data corresponding to the field variables, and add index values
      .col.names <- .opt.cor.names[rownames(.opt.cor.names) %in% variables, ,
                                   drop = FALSE]
      .index <- correlations[[.i]][[.j]][, .plot.design[.j]]
      .opt.cor <- correlations[[.i]][[.j]][, .col.names[, "cor"], drop = FALSE]
      .opt.metric <- correlations[[.i]][[.j]][, .col.names[, "metric"],
                                              drop = FALSE]

      # Define title, subtitle, and axis names
      .title <- switch(.i, pearson = "Pearson correlation",
                       spearman = "Spearman correlation")
      .subtitle <- paste("<br> <span style='font-size: 12px;'>",
                         switch(.j, fixed.area.plot = "Fixed area",
                                k.tree.plot = "K-tree",
                                angle.count.plot = "Angle-count"),
                         "</span>", sep ="")
      .xaxis <- switch(.j, fixed.area.plot = "Radius (m)",
                       k.tree.plot = "K-tree (trees)",
                       angle.count.plot = "BAF (m<sup>2</sup>/ha)")
      .yaxis <- "Variables"

      # Create heatmap
      fig <- plotly::plot_ly(x = .index, y = rownames(.col.names),
                             z = t(.opt.cor), type = "heatmap",
                             zmin = -1, zmax = 1, zmid = 0, colors = .color,
                             hoverinfo = "x+y+text",
                             text = matrix(paste("Cor:", round(t(.opt.cor), 3),
                                                 "<br>Metric:", t(.opt.metric)),
                                           nrow = ncol(.opt.cor),
                                           ncol = nrow(.opt.cor))) %>%
        plotly::layout(title = paste(.title, .subtitle, sep = ""),
                       xaxis = list(title = .xaxis),
                       yaxis = list (title = .yaxis), margin = list(t = 50))

      # Save heatmap
      suppressWarnings(
        htmlwidgets::saveWidget(widget = fig,
                                file = file.path(dir.result,
                                                 paste("opt.correlations", .j,
                                                       .i, "html", sep = ".")),
                                selfcontained = TRUE,
                                libdir = file.path(dir.result,
                                                   "opt.correlations_files")))

    }

    # Define final time, and print message
    t1 <- Sys.time()
    message(" (", format(round(difftime(t1, t0, units = "secs"), 2)), ")")

  }


  return(invisible())

}
