
tree.detection.several.plots <- function(las.list, scan.approach = "single",

                                         id = NULL, file = NULL,

                                    normalize.arguments =

                                      list(x.center = NULL, y.center = NULL,

                                           max.dist = NULL, min.height = NULL, max.height = NULL,

                                           algorithm.dtm = "tin", res.dtm = 0.2),

                                    tree.detection.arguments =

                                      list(dbh.min = 7.5, dbh.max = 200, h.min = 1.3,

                                           ncr.threshold = 0.1,

                                           tls.resolution = list(),

                                           tls.precision = NULL,

                                           breaks = c(1.0, 1.3, 1.6),

                                           plot.attributes = NULL),

                                    dir.data = NULL, save.result = TRUE, dir.result = NULL){


  # DTM algorhitm

  if(is.null(normalize.arguments$algorithm.dtm)){

    .algorithm.dtm <- "tin"

  } else {

    .algorithm.dtm <- normalize.arguments$algorithm.dtm

  }


  if(is.null(normalize.arguments$res.dtm)){

    .res.dtm <- 0.2

  } else {

    .res.dtm <- normalize.arguments$res.dtm

  }

  # Tree detection arguments

  if(is.null(tree.detection.arguments$tls.precision)){

    .tls.precision <- NULL

  } else {

    .tls.precision <- tree.detection.arguments$tls.precision

  }

  if(is.null(tree.detection.arguments$dbh.min)){

    .dbh.min <- 7.5

  } else {

    .dbh.min <- tree.detection.arguments$dbh.min

  }

  if(is.null(tree.detection.arguments$dbh.max)){

    .dbh.max <- 200

  } else {

    .dbh.max <- tree.detection.arguments$dbh.max

  }

  if(is.null(tree.detection.arguments$h.min)){

    .h.min <- 1.3

  } else {

    h.min <- tree.detection.arguments$h.min

  }


  if(is.null(tree.detection.arguments$ncr.threshold)){

    .ncr.threshold <- 0.1

  } else {

    .ncr.threshold <- tree.detection.arguments$ncr.threshold

  }

  if(is.null(tree.detection.arguments$breaks)){

    .breaks <- c(1.0, 1.3, 1.6)

  } else {

    .breaks <- tree.detection.arguments$breaks

  }


  for (i in (1:length(las.list))) {

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

                       x.center = normalize.arguments$x.center,

                       y.center = normalize.arguments$y.center,

                       max.dist = normalize.arguments$max.dist,

                       min.height = normalize.arguments$min.height,

                       max.height = normalize.arguments$max.height,

                       algorithm.dtm = .algorithm.dtm,

                       res.dtm = .res.dtm,

                       scan.approach = scan.approach,

                       id = .id, file = .file,

                       dir.data = dir.data, save.result = save.result, dir.result = dir.result)

    message("Detecting trees")

    if(scan.approach == "single"){

    .tree.list.tls.i <- tree.detection.single.scan(data = .data,

                                       dbh.min = .dbh.min,

                                       dbh.max = .dbh.max,

                                       h.min = .h.min,

                                       breaks = .breaks,

                                       ncr.threshold = .ncr.threshold,

                                       tls.resolution = tree.detection.arguments$tls.resolution,

                                       plot.attributes = tree.detection.arguments$plot.attributes,

                                       save.result = FALSE, dir.result = dir.result)}

    if(scan.approach == "multi"){

    .tree.list.tls.i <- tree.detection.multi.scans(data = .data,

                                                   dbh.min = .dbh.min,

                                                   dbh.max = .dbh.max,

                                                   h.min = .h.min,

                                                   breaks = .breaks,

                                                   ncr.threshold = .ncr.threshold,

                                                   tls.precision = .tls.precision,

                                                   plot.attributes = tree.detection.arguments$plot.attributes,

                                                   save.result = FALSE, dir.result = dir.result)

    }

    if(i < 2){

      .tree.list.tls <- .tree.list.tls.i

    } else {

      .tree.list.tls <- rbind(.tree.list.tls, .tree.list.tls.i)

    }

    if(isTRUE(save.result)){

      utils::write.csv(.tree.list.tls,
                       file = file.path(dir.result, "tree.list.tls.csv"),
                       row.names = FALSE)

    }


  }

  #####
  return(.tree.list.tls)

}
