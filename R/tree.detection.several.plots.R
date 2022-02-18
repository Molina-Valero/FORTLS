
tree.detection.several.plots <- function(las.list, scan.approach = "single",

                                         id = NULL, file = NULL,


                                         x.center = NULL, y.center = NULL,

                                         max.dist = NULL, min.height = NULL, max.height = NULL,

                                         algorithm.dtm = "tin", res.dtm = 0.2,


                                         dbh.min = 7.5, dbh.max = 200, h.min = 1.3,

                                         ncr.threshold = 0.1,

                                         tls.resolution = NULL, tls.precision = NULL,

                                         breaks = c(1.0, 1.3, 1.6),


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

    .data <- normalize(las = las.list[[i]],

                       x.center = x.center, y.center = y.center,

                       max.dist = max.dist,

                       min.height = min.height, max.height = max.height,

                       algorithm.dtm = algorithm.dtm,

                       res.dtm = res.dtm,

                       scan.approach = scan.approach,

                       id = .id, file = .file,

                       dir.data = dir.data, save.result = save.result, dir.result = dir.result)

    message("Detecting trees")

    if(scan.approach == "single"){

    .tree.tls.i <- tree.detection.single.scan(data = .data,

                                              dbh.min = dbh.min,

                                              dbh.max = dbh.max,

                                              h.min = h.min,

                                              breaks = breaks,

                                              ncr.threshold = ncr.threshold,

                                              tls.resolution = tls.resolution,

                                              plot.attributes = plot.attributes,

                                              save.result = FALSE, dir.result = dir.result)}

    if(scan.approach == "multi"){

    .tree.tls.i <- tree.detection.multi.scans(data = .data,

                                              dbh.min = dbh.min,

                                              dbh.max = dbh.max,

                                              h.min = h.min,

                                              breaks = breaks,

                                              ncr.threshold = ncr.threshold,

                                              tls.precision = tls.precision,

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
