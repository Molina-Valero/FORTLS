

#................................................................................ secondary function to use and sort the sublists before using as input to the Rcpp function

sort_sublists = function(lst_unsorted,
                         MIN_PNTS_VX = 2) {              # keep sublists with at least that many observations

  ROWS = lapply(lst_unsorted, nrow) |>
    unlist() |>
    as.vector()

  idx_keep = which(ROWS > MIN_PNTS_VX)

  lst_unsorted = lst_unsorted[idx_keep]
  ROWS = ROWS[idx_keep]

  # rank the rows to process the smaller objects first
  ROWS_rank = rank(ROWS, ties.method = 'first')
  rank_lst = list(rank = ROWS_rank,
                  item = 1:length(ROWS_rank)) |>
    data.table::setDT()
  rank_lst = rank_lst[order(rank_lst$rank, decreasing = F), ]

  lst_unsorted = lst_unsorted[rank_lst$item]
  return(lst_unsorted)
}



geometric.features <- function(data,
                               dist,
                               approximate_KNN = FALSE,
                               features = c("PCA1", "PCA2", "anisotropy", "planarity",
                                            "linearity", "surface_variation", "sphericity", "verticality",
                                            "number_neighbors", "omnivariance", "eigenentropy", "surface_density",
                                            "volume_density"),
                               threads = 1,
                               keep_NaN = FALSE,
                               verbose = FALSE,
                               solver_threshold = 50000) {

  valid_features = c("PCA1", "PCA2", "anisotropy", "planarity", "linearity",          # The "first_eigenvalue", "second_eigenvalue", "third_eigenvalue", "sum_eigenvalues" are not features because the other features depend on these for the computation thus they have to be returned by default
                     "surface_variation", "sphericity", "verticality",
                     "number_neighbors", "omnivariance", "eigenentropy",
                     "surface_density", "volume_density")

  verify = (features %in% valid_features)
  if (length(verify) == 0) stop(glue::glue("Invalid 'features' input! Select one of the {paste(valid_features, collapse = ', ')}"))

  if (!all(verify)) {
    stop(glue::glue("Valid features are: {paste(valid_features, collapse = ', ')}"), '\n')
  }

  # Select necessary fields from original txt file of point cloud
  if (!inherits(data, 'data.table')) {
    data = data.table::as.data.table(data[, c("point", "x", "y", "z", 'id')])
  }

  # Create x and y coordinates for grid
  x <- seq(min(data$x), max(data$x)+0.5, by = 0.5)
  y <- seq(min(data$y), max(data$y)+0.5, by = 0.5)


  if (length(x) < 2 | length(y) < 2) {

    if (verbose) cat(glue::glue("The dbscan::frNN() will be used to find the nearest points per row using a distance of '{dist}' ..."), '\n')

    mt_knn = as.matrix(dat[, c("x", "y", "z")])

    res_ann = dbscan::frNN(
      x = mt_knn,
      eps = dist,
      query = NULL,
      sort = FALSE,
      search = "kdtree",
      bucketSize = 10,
      splitRule = "suggest",
      approx = ifelse(approximate_KNN, 1, 0)
    )

    knn_idxs = res_ann$id


    geo.fea = subset_matrix_by_indices(mat = as.matrix(dat),
                                       index_lists = knn_idxs,
                                       dist = dist,
                                       First_eigenvalue = TRUE,
                                       Second_eigenvalue = TRUE,
                                       Third_eigenvalue = TRUE,
                                       Sum_of_eigenvalues = TRUE,
                                       PCA_1 = ifelse("PCA1" %in% features, TRUE, FALSE),
                                       PCA_2 = ifelse("PCA2" %in% features, TRUE, FALSE),
                                       Anisotropy = ifelse("anisotropy" %in% features, TRUE, FALSE),
                                       Planarity = ifelse("planarity" %in% features, TRUE, FALSE),
                                       Linearity = ifelse("linearity" %in% features, TRUE, FALSE),
                                       Surface_variation = ifelse("surface_variation" %in% features, TRUE, FALSE),
                                       Normal_change_rate = ifelse("sphericity" %in% features, TRUE, FALSE),
                                       Verticality = ifelse("verticality" %in% features, TRUE, FALSE),
                                       Number_of_points = ifelse("number_neighbors" %in% features, TRUE, FALSE),
                                       omnivariance = ifelse("omnivariance" %in% features, TRUE, FALSE),
                                       eigenentropy = ifelse("eigenentropy" %in% features, TRUE, FALSE),
                                       surface_density = ifelse("surface_density" %in% features, TRUE, FALSE),
                                       volume_density = ifelse("volume_density" %in% features, TRUE, FALSE),
                                       threads = threads,
                                       solver_thresh = solver_threshold)

    geo.fea = do.call(rbind, geo.fea) |>
      data.table::as.data.table()

    # required because we have different names in the cpp code and we use different names in the R code [ we might have to adjust the names in the R code ]
    for (NAME in features) {
      if (NAME == 'anisotropy') colnames(geo.fea)[which(colnames(geo.fea) == 'Anisotropy')] = 'anisotropy'
      if (NAME == 'planarity') colnames(geo.fea)[which(colnames(geo.fea) == 'Planarity')] = 'planarity'
      if (NAME == 'linearity') colnames(geo.fea)[which(colnames(geo.fea) == 'Linearity')] = 'linearity'
      if (NAME == 'surface_variation') colnames(geo.fea)[which(colnames(geo.fea) == 'Surface_variation')] = 'surface_variation'
      if (NAME == 'sphericity') colnames(geo.fea)[which(colnames(geo.fea) == 'Normal_change_rate')] = 'sphericity'
      if (NAME == 'verticality') colnames(geo.fea)[which(colnames(geo.fea) == 'Verticality')] = 'verticality'
      if (NAME == 'number_neighbors') colnames(geo.fea)[which(colnames(geo.fea) == 'Number_of_points')] = 'number_neighbors'
    }

    cols_keep = c('point', features)

    geo.fea = geo.fea[, ..cols_keep]

    if (!keep_NaN) {
      geo.fea = geo.fea[stats::complete.cases(geo.fea), , drop = F]
    }

  } else{

    if (verbose) cat(glue::glue("The dbscan::frNN() will be used to find the nearest points per row using a distance of '{dist}' ..."), '\n')

    # if (verbose) cat('A grid will be created ...\n')
    # Empty data frame where coordinates neccesaries for creating grid will be saved
    grid <- data.frame(x = rep(x, each = length(y)),
                       y = rep(y, times = length(x)))

    grid <- sf::st_as_sf(grid, coords = c("x","y"))

    # if (verbose) cat('A buffer will be used for overlapping grid-cells ...\n')
    grid <- sf::st_buffer(grid, dist = 0.25 + dist, endCapStyle = "SQUARE")
    grid <- sf::st_cast(grid, "POLYGON")

    # if (verbose) cat('Intersection between the initial xyz coordinates and the grid-cells ...\n')
    grid.2 <- sf::st_intersects(grid, sf::st_as_sf(data, coords = c("x","y")))
    grid.2 <- as.data.frame(grid.2)
    colnames(grid.2) <- c("id", "code")
    data$code <- as.numeric(row.names(data))

    data <- merge(data[, c('code', 'point', 'x', 'y', 'z')], grid.2, by = "code")
    data <- data[, 2:ncol(data)]
    rm(grid, grid.2)

    dat <- split(data[, c('point', 'x', 'y', 'z', 'id')], by = 'id')
    dat = sort_sublists(lst_unsorted = dat, MIN_PNTS_VX = 2)
    # if (verbose) cat(glue::glue("{length(dat)} grid-cells will be processed ..."), '\n')

    LEN = length(dat)

    t_start = proc.time()
    geo.fea = lapply(1:LEN, function(y) {
      if (verbose) {
        if (y %in% as.integer(seq(from = 1, to = LEN, length.out = 10))) {
          cat(glue::glue("Geometric properties processing: {round(y / LEN, 4) * 100} %"), "\n")
        }
      }

      x = dat[[y]]

      if ('id' %in% colnames(x)) {
        x$id = NULL

      }


    mt_knn = as.matrix(x[, c("x", "y", "z")])

    res_ann = dbscan::frNN(
      x = mt_knn,
      eps = dist,
      query = NULL,
      sort = FALSE,
      search = "kdtree",
      bucketSize = 10,
      splitRule = "suggest",
      approx = ifelse(approximate_KNN, 1, 0)
    )

    knn_idxs = res_ann$id

    geo.fea = subset_matrix_by_indices(mat = as.matrix(x[, c("point", "x", "y", "z")]),
                                       index_lists = knn_idxs,
                                       dist = dist,
                                       First_eigenvalue = TRUE,
                                       Second_eigenvalue = TRUE,
                                       Third_eigenvalue = TRUE,
                                       Sum_of_eigenvalues = TRUE,
                                       PCA_1 = ifelse("PCA1" %in% features, TRUE, FALSE),
                                       PCA_2 = ifelse("PCA2" %in% features, TRUE, FALSE),
                                       Anisotropy = ifelse("anisotropy" %in% features, TRUE, FALSE),
                                       Planarity = ifelse("planarity" %in% features, TRUE, FALSE),
                                       Linearity = ifelse("linearity" %in% features, TRUE, FALSE),
                                       Surface_variation = ifelse("surface_variation" %in% features, TRUE, FALSE),
                                       Normal_change_rate = ifelse("sphericity" %in% features, TRUE, FALSE),
                                       Verticality = ifelse("verticality" %in% features, TRUE, FALSE),
                                       Number_of_points = ifelse("number_neighbors" %in% features, TRUE, FALSE),
                                       omnivariance = ifelse("omnivariance" %in% features, TRUE, FALSE),
                                       eigenentropy = ifelse("eigenentropy" %in% features, TRUE, FALSE),
                                       surface_density = ifelse("surface_density" %in% features, TRUE, FALSE),
                                       volume_density = ifelse("volume_density" %in% features, TRUE, FALSE),
                                       threads = threads,
                                       solver_thresh = solver_threshold)

    geo.fea = do.call(rbind, geo.fea) |>
      data.table::as.data.table()

    # required because we have different names in the cpp code and we use different names in the R code [ we might have to adjust the names in the R code ]
    for (NAME in features) {
      if (NAME == 'anisotropy') colnames(geo.fea)[which(colnames(geo.fea) == 'Anisotropy')] = 'anisotropy'
      if (NAME == 'planarity') colnames(geo.fea)[which(colnames(geo.fea) == 'Planarity')] = 'planarity'
      if (NAME == 'linearity') colnames(geo.fea)[which(colnames(geo.fea) == 'Linearity')] = 'linearity'
      if (NAME == 'surface_variation') colnames(geo.fea)[which(colnames(geo.fea) == 'Surface_variation')] = 'surface_variation'
      if (NAME == 'sphericity') colnames(geo.fea)[which(colnames(geo.fea) == 'Normal_change_rate')] = 'sphericity'
      if (NAME == 'verticality') colnames(geo.fea)[which(colnames(geo.fea) == 'Verticality')] = 'verticality'
      if (NAME == 'number_neighbors') colnames(geo.fea)[which(colnames(geo.fea) == 'Number_of_points')] = 'number_neighbors'
    }

    cols_keep = c('point', features)

    geo.fea = geo.fea[, ..cols_keep]

    if (!keep_NaN) {
      geo.fea = geo.fea[stats::complete.cases(geo.fea), , drop = F]
    }

    })
  }

  t_end = proc.time()
  if (verbose) {
    cat(glue::glue("Elapsed time for Geometric Features computation (in seconds): {round((t_end - t_start)['elapsed'], digits = 1)}"), '\n')
  }

  geo.fea = data.table::rbindlist(geo.fea)

  geo.fea <- geo.fea[!duplicated(geo.fea$point), ]

  geo.fea <- merge(data[, c("point", "x", "y", "z")], geo.fea, by = "point", all.x = TRUE)

  geo.fea <- geo.fea[!duplicated(geo.fea$point), ]

  return(geo.fea)

}


