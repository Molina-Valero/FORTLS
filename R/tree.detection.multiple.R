
tree.detection.multiple <- function(las.list, id = NULL, file = NULL,

                                    normalize.arguments =

                                      list(max.dist = NULL, min.height = NULL, max.height = NULL),

                                    tree.detection.arguments =

                                      list(dbh.min = 7.5, dbh.max = 200,

                                           ncr.threshold = 0.1,

                                           tls.resolution = list(),

                                           breaks = c(1.0, 1.3, 1.6),

                                           plot.attributes = NULL),

                                    dir.data = NULL, save.result = TRUE, dir.result = NULL){

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

    message("Computing plot: ", i)

    message("Normalizing")

    # Assign id

    if(!is.null(id)){

      .id <- id[i]

    } else {

      .id <- i

    }


    # File name

    if(!is.null(file)){

      .file <- paste(file[i], ".txt", sep = "")

    } else {

      .file <- paste(.id, ".txt", sep = "")

    }

    .data <- normalize(las = las.list[[i]],

                       max.dist = normalize.arguments$max.dist,

                       min.height = normalize.arguments$min.height,

                       max.height = normalize.arguments$max.height,

                       id = .id, file = .file,

                       dir.data = dir.data, save.result = save.result, dir.result = dir.result)

    message("Detecting trees")

    .tree.list.tls.i <- tree.detection(data = .data,

                                       dbh.min = tree.detection.arguments$dbh.min,

                                       dbh.max = tree.detection.arguments$dbh.max,

                                       breaks = .breaks,

                                       ncr.threshold = .ncr.threshold,

                                       tls.resolution = tree.detection.arguments$tls.resolution,

                                       plot.attributes = tree.detection.arguments$plot.attributes,

                                       save.result = FALSE, dir.result = dir.result)

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
