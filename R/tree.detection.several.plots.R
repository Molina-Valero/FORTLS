
tree.detection.several.plots <- function(las.list, id.list = NULL, file = NULL,

                                         scan.approach = "single",

                                         pcd.red = NULL, normalized = NULL,
                                         center.coord = NULL,
                                         x.side = NULL, y.side = NULL,
                                         max.dist = NULL, min.height = NULL, max.height = 50,
                                         algorithm.dtm = "knnidw", res.dtm = 0.2, csf = list(cloth_resolution = 0.5),
                                         intensity = NULL, RGB = NULL,

                                         single.tree = NULL,
                                         dbh.min = 4, dbh.max = 200, h.min = 1.3,
                                         ncr.threshold = 0.1,
                                         tls.resolution = NULL, tls.precision = NULL,
                                         density.reduction = 2,
                                         stem.section = c(0.7, 3.5), stem.range = NULL, breaks = NULL,
                                         slice = 0.1, understory = NULL, bark.roughness = 1,
                                         den.type = 1, d.top = NULL,
                                         plot.attributes = NULL, plot = NULL,

                                         dir.data = NULL, save.result = TRUE, dir.result = NULL){

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()


  for (i in 1:length(las.list)) {

    # Assign id

    if(!is.null(id.list)){

      .id <- id.list[i]

    } else {

      .id <- as.integer(i)

    }

    # Assign coordinates

    if(!is.null(center.coord) & !is.null(center.coord$id)){

      x.center <- center.coord[center.coord$id == .id, "x"]
      y.center <- center.coord[center.coord$id == .id, "y"]

    } else if(!is.null(center.coord)) {

      x.center <- center.coord$x[i]
      y.center <- center.coord$y[i]

    } else {

      x.center = NULL
      y.center = NULL

    }



    # Messages

    message("Computing plot: ", .id)

    message("Normalizing")

    # File name

    if(!is.null(file)){

      .file <- paste(file[i], ".txt", sep = "")

    } else {

      .file <- paste(.id, ".txt", sep = "")

    }

    .data <- normalize(las = las.list[i], normalized = normalized,

                       x.center = x.center, y.center = y.center,

                       max.dist = max.dist, min.height = min.height, max.height = max.height,

                       algorithm.dtm = algorithm.dtm, res.dtm = res.dtm, csf = csf,

                       intensity = intensity, RGB = RGB,

                       scan.approach = scan.approach,

                       id = .id, file = .file, plot = NULL,

                       dir.data = dir.data, save.result = save.result, dir.result = dir.result)

    if(!is.null(pcd.red))
      .data <- .data[.data$prob.selec == 1, ]

    message("Detecting trees")

    if(scan.approach == "single"){

    if(i > 1)
      .stem.curve <- read.csv("stem.curve.csv")

    .tree.tls.i <- tree.detection.single.scan(data = .data, single.tree = single.tree,

                                              dbh.min = dbh.min, dbh.max = dbh.max, h.min = h.min,

                                              ncr.threshold = ncr.threshold,

                                              tls.resolution = tls.resolution,

                                              density.reduction = density.reduction,

                                              stem.section = stem.section, stem.range = stem.range, breaks = breaks,

                                              slice = slice, understory = understory, bark.roughness = bark.roughness,

                                              den.type = den.type, d.top = d.top,

                                              plot.attributes = plot.attributes, plot = plot,

                                              save.result = FALSE, dir.result = dir.result)

    rm(.data)

    }

    if(scan.approach == "multi"){

    if(i > 1)
      .stem.curve <- read.csv("stem.curve.csv")

    .tree.tls.i <- tree.detection.multi.scan(data = .data, single.tree = single.tree,

                                             dbh.min = dbh.min, dbh.max = dbh.max, h.min = h.min,

                                             ncr.threshold = ncr.threshold,

                                             tls.precision = tls.precision,

                                             density.reduction = density.reduction,

                                             stem.section = stem.section, stem.range = stem.range, breaks = breaks,

                                             slice = slice, understory = understory, bark.roughness = bark.roughness,

                                             den.type = den.type, d.top = d.top,

                                             plot.attributes = plot.attributes, plot = plot,

                                             save.result = FALSE, dir.result = dir.result)

    rm(.data)

    }

    if(i < 2){

      .tree.tls <- .tree.tls.i

    } else {

      .tree.tls <- rbind(.tree.tls, .tree.tls.i)

      .stem.curve <- rbind(.stem.curve, read.csv("stem.curve.csv"))

      utils::write.csv(.stem.curve,
                       file = file.path(dir.result, "stem.curve.csv"),
                       row.names = FALSE)

    }

    if(isTRUE(save.result)){

      utils::write.csv(.tree.tls,
                       file = file.path(dir.result, "tree.tls.csv"),
                       row.names = FALSE)

    }


  }


  #####
  return(.tree.tls)

}

