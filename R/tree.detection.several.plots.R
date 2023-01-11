
tree.detection.several.plots <- function(las.list, id = NULL, file = NULL,

                                         scan.approach = "single",

                                         pcd.red = NULL, normalized = NULL,
                                         x.center = NULL, y.center = NULL,
                                         x.side = NULL, y.side = NULL,
                                         max.dist = NULL, min.height = NULL, max.height = 50,
                                         algorithm.dtm = "knnidw", res.dtm = 0.2, csf = list(cloth_resolution = 0.5),
                                         RGB = NULL,

                                         single.tree = NULL,
                                         dbh.min = 4, dbh.max = 200, h.min = 1.3,
                                         ncr.threshold = 0.1,
                                         tls.resolution = NULL, tls.precision = NULL,
                                         stem.section = NULL, breaks = NULL,
                                         slice = 0.1, understory = NULL, bark.roughness = 2,
                                         den.type = 1, d.top = NULL,
                                         plot.attributes = NULL,

                                         dir.data = NULL, save.result = TRUE, dir.result = NULL){

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()


  for (i in 1:length(las.list)) {

    # Assign id

    if(!is.null(id)){

      .id <- id[i]

    } else {

      .id <- i

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

                       RGB = RGB,

                       scan.approach = scan.approach,

                       id = .id, file = .file,

                       dir.data = dir.data, save.result = save.result, dir.result = dir.result)

    if(!is.null(pcd.red))
      .data <- .data[.data$prob.selec == 1, ]

    message("Detecting trees")

    if(scan.approach == "single"){

    .tree.tls.i <- tree.detection.single.scan(data = .data,

                                              dbh.min = dbh.min, dbh.max = dbh.max, h.min = h.min,

                                              ncr.threshold = ncr.threshold,

                                              tls.resolution = tls.resolution,

                                              stem.section = stem.section, breaks = breaks,

                                              slice = slice, understory = understory, bark.roughness = bark.roughness,

                                              den.type = den.type, d.top = d.top,

                                              plot.attributes = plot.attributes,

                                              save.result = FALSE, dir.result = dir.result)}

    if(scan.approach == "multi"){

    .tree.tls.i <- tree.detection.multi.scan(data = .data, single.tree = single.tree,

                                             dbh.min = dbh.min, dbh.max = dbh.max, h.min = h.min,

                                             ncr.threshold = ncr.threshold,

                                             tls.precision = tls.precision,

                                             stem.section = stem.section, breaks = breaks,

                                             slice = slice, understory = understory, bark.roughness = bark.roughness,

                                             den.type = den.type, d.top = d.top,

                                             plot.attributes = plot.attributes,

                                             save.result = FALSE, dir.result = dir.result)

    }

    if(i < 2){

      .tree.tls <- .tree.tls.i

    } else {

      .tree.tls <- rbind(.tree.tls, .tree.tls.i)

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
