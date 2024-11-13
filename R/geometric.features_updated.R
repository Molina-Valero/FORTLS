
# required R packages
# library(lidR)
# library(data.table)
# library(VoxR)
# library(glue)
# library(sf)
# library(Rfast)
# library(dplyr)
# library(mapview)

#................................................................................. The updated geometric.features function

geometric.features <- function(data,
                               grid_method,                    # either 'sf_grid', 'voxel_grid', 'data_table'
                               voxel_resolution = 6.0,         # stands for resolution in meters if the 'voxel_grid' method is selected
                               data_table_grid_cells = 20,     # grid cell matrix of maximum size (data_table_grid_cells, data_table_grid_cells)
                               features = c("PCA1", "PCA2", "anisotropy", "planarity",
                                            "linearity", "surface_variation", "sphericity", "verticality",
                                            "number_neighbors", "omnivariance", "eigenentropy", "surface_density",
                                            "volume_density"),
                               dist = 0.05,
                               threads = 1,
                               keep_NaN = FALSE,
                               verbose = FALSE,
                               solver_threshold = 50000) {

  # if (grid_method == 'sf_grid' & !('code' %in% colnames(data))) stop("The 'sf_grid' method requires that the 'code' column name exists in the data!")

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
  x <- seq(min(data$x), max(data$x)+1)
  y <- seq(min(data$y), max(data$y)+1)

  if (length(x) < 2 | length(y) < 2) {

    dat <- data[, c("point", "x", "y", "z")]

    geo.fea <- geometric_features(m = as.matrix(dat),
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

    geo.fea = geo.fea[, c('point', features)]
    geo.fea = data.table::as.data.table(geo.fea)

  } else{

    if (verbose) cat(glue::glue("The '{grid_method}' method will be used ..."), '\n')

    if (grid_method == 'voxel_grid') {

      vx = VoxR_vox_update(data = data[, c( 'x', 'y', 'z')],
                           res = voxel_resolution,
                           message = T,
                           full.grid = F)       # don't set 'full.grid' to TRUE because for some reason it doesn't return the 'data_out' sublist

      vx_voxels = vx$data_out
      vx_voxels$id = glue::glue("{vx_voxels$x}_{vx_voxels$y}_{vx_voxels$z}")

      vx_init = vx$all_data
      vx_init$id = glue::glue("{vx_init$x}_{vx_init$y}_{vx_init$z}")

      if (length(setdiff(x = unique(vx_voxels$id), y = vx_init$id)) > 0) {
        stop("We don't expect a difference in the voxels 'id' column between the data.tables!")
      }

      dat_3d_compute = data[, c('point', 'x', 'y', 'z')]

      if (nrow(dat_3d_compute) != nrow(vx_init)) {
        stop("We expect the same number of rows between computed voxel data.tables!")
      }

      dat_3d_compute$id = vx_init$id

      if (verbose) cat('Split input data by voxel-id ...\n')
      dat = split(dat_3d_compute, by = 'id')

      if (verbose) cat(glue::glue("{length(dat)} voxels will be processed ..."), '\n')

      dat = sort_sublists(lst_unsorted = dat, MIN_PNTS_VX = 2)

    } else if (grid_method == 'sf_grid') {

      if (verbose) cat('A grid will be created ...\n')
      # Empty data frame where coordinates neccesaries for creating grid will be saved
      grid <- data.frame(x = rep(x, each = length(y)),
                         y = rep(y, times = length(x)))

      grid <- sf::st_as_sf(grid, coords = c("x","y"))

      if (verbose) cat('A buffer will be used for overlapping grid-cells ...\n')
      grid <- sf::st_buffer(grid, dist = 0.5 + dist, endCapStyle = "SQUARE")
      grid <- sf::st_cast(grid, "POLYGON")
      # mapview::mapview(grid)

      if (verbose) cat('Intersection between the initial xyz coordinates and the grid-cells ...\n')
      grid.2 <- sf::st_intersects(grid, sf::st_as_sf(data, coords = c("x","y")))
      grid.2 <- as.data.frame(grid.2)
      colnames(grid.2) <- c("id", "code")
      data$code <- as.numeric(row.names(data))

      data <- merge(data[, c('code', 'point', 'x', 'y', 'z')], grid.2, by = "code")
      data <- data[, 2:ncol(data)]
      rm(grid, grid.2)

      dat <- split(data[, c('point', 'x', 'y', 'z', 'id')], by = 'id')
      dat = sort_sublists(lst_unsorted = dat, MIN_PNTS_VX = 2)
      if (verbose) cat(glue::glue("{length(dat)} grid-cells will be processed ..."), '\n')

    } else if (grid_method == 'data_table') {

      dat = datatable_grid(data_inp = data,
                               dist = dist,
                               n_grid = data_table_grid_cells,
                               verify_visually = FALSE)

      dat = dat$data[, c('I', 'x', 'y', 'z', 'x_bin', 'y_bin')]
      dat$id = glue::glue("{dat$x_bin}_{dat$y_bin}")
      dat$point = data$point[dat$I]                     # the 'point' column exists and it will be re-ordered based on 'I'
      dat = dat[, !c('x_bin', 'y_bin')]

      dat <- split(dat, by = 'id')
      dat = sort_sublists(lst_unsorted = dat, MIN_PNTS_VX = 2)
      if (verbose) cat(glue::glue("{length(dat)} grid-cells will be processed ..."), '\n')

    } else {
      stop("The grid method must be either 'sf_grid' or 'voxel_grid'!")
    }

    LEN = length(dat)

    if (verbose) cat(glue::glue("The Geometric Features computation starts using {threads} threads ..."), '\n')

    t_start = proc.time()
    geo.fea = lapply(1:LEN, function(y) {
      if (verbose) {
        if (y %in% as.integer(seq(from = 1, to = LEN, length.out = 10))) {
          cat(glue::glue("Processed chunks: {round(y / LEN, 4) * 100} %"), "\n")
        }
      }
      x = dat[[y]]

      if ('id' %in% colnames(x)) {
        x$id = NULL
      }

      x = geometric_features(m = as.matrix(x),
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

      x = x[, c('point', features)]

      if (!keep_NaN) {
        x = x[complete.cases(x), , drop = F]
      }

      x = x |>
        data.table::as.data.table()
      x
    })
    t_end = proc.time()
    if (verbose) {
      cat(glue::glue("Elapsed time for Geometric Features computation (in seconds): {round((t_end - t_start)['elapsed'], digits = 1)}"), '\n')
    }

    geo.fea = data.table::rbindlist(geo.fea)

    if (verbose) {
      cat(glue::glue("Processed rows (valid computations) for {round(nrow(geo.fea) / nrow(data), 4) * 100} % of the input data!"), "\n")
    }
  }

  if (nrow(geo.fea) > 0) {
    data <- merge(data[, c("point", "x", "y", "z")], geo.fea, by = "point", all.x = TRUE)
  }

  data <- data[!duplicated(data), ]

  return(data)
}

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

#................................................................................. # This function was updated to return also the input data with the computed voxels
# reference: https://github.com/Blecigne/VoxR/blob/Blecigne/VoxR/R/vox.R

VoxR_vox_update = function (data,
                            res,
                            full.grid,
                            message) {

  x = y = z = npts = .N = . = ":=" = NULL
  check = VoxR::ck_conv_dat(data, message = message)
  if (missing(res)) {
    stop("No voxel resolution (res) provided")
  }
  else {
    if (!is.vector(res))
      stop("res must be a vector of length 1")
    if (!is.numeric(res))
      stop("res must be numeric")
    if (res <= 0)
      stop("res must be positive")
    if (length(res) > 1) {
      res = res[1]
      warning("res contains more than 1 element. Only the first was used")
    }
  }
  if (missing(full.grid))
    full.grid = FALSE
  data = check$data

  data[, `:=`(x = Rfast::Round(x/res) * res, y = Rfast::Round(y/res) * res, z = Rfast::Round(z/res) * res)]
  all_copy = data.table::copy(data)
  data = unique(data[, `:=`(npts, .N), by = .(x, y, z)])
  if (full.grid) {
    x_seq = seq(min(data$x), max(data$x), res)
    y_seq = seq(min(data$y), max(data$y), res)
    z_seq = seq(min(data$z), max(data$z), res)
    est_weight = round(length(x_seq) * length(y_seq) * length(z_seq) *
                         4 * 8/2^{
                           20
                         }/1024, 1)
    if (est_weight > 2) {
      cat(paste("Final data is estimated to be ", est_weight,
                "GB in size. Continue execution ? y or n"))
      test = readline()
      if (test != "y") {
        stop("Execution stoped.")
      }
    }
    empty = data.table::data.table(expand.grid(x_seq, y_seq,
                                               z_seq))
    data.table::setnames(empty, c("x", "y", "z"))
    empty[, `:=`(npts, 0)]
    data = dplyr::bind_rows(data, empty)
    data = data[, `:=`(npts, sum(npts)), keyby = .(x, y,
                                                   z)]
  }
  if (check$dfr)
    data = as.data.frame(data)
  return(list(data_out = data,
              all_data = all_copy))
}


#...............................................................................  compute the chunk size from the area of the las-catalog

chunk_size_from_area = function(las_catalog,
                                num_parts) {
  # Get the extent of the file
  EXT = terra::ext(las_catalog) |> sf::st_bbox()

  # # Calculate the total area
  total_area <- (as.numeric(EXT['xmax']) - as.numeric(EXT['xmin'])) * (as.numeric(EXT['ymax']) - as.numeric(EXT['ymin']))

  # Calculate the area of each part
  area_m_sq <- total_area / num_parts

  # # Calculate the chunk size (assuming square chunks)
  chunk_size <- sqrt(area_m_sq)

  return(chunk_size)
}

# chunk_size_from_area = function(las_catalog,
#                                 num_parts) {
#
#   # Get the extent of the file
#   ext <- lidR::extent(las_catalog)
#
#   # # Calculate the total area
#   total_area <- (ext@xmax - ext@xmin) * (ext@ymax - ext@ymin)
#
#   # Calculate the area of each part
#   area_m_sq <- total_area / num_parts
#
#   # # Calculate the chunk size (assuming square chunks)
#   chunk_size <- sqrt(area_m_sq)
#
#   return(chunk_size)
# }

#...............................................................................  save to tiles and receive the metadata

save_to_tiles = function(las_catalog,
                         chunk_size_input,
                         dir_save_tiles,
                         buffer_size = 2,
                         increment_chunk_size_factor = 1.4) {

  # Set the chunk size option
  opt_chunk_size(las_catalog) <- chunk_size_input * increment_chunk_size_factor                  # unfortunately this does not give the exact number of tiles

  chunk_buffer <- buffer_size                                                                    # in the same units as your coordinate system
  opt_chunk_buffer(las_catalog) <- chunk_buffer
  opt_output_files(las_catalog) <- file.path(dir_save_tiles, "/retile_{XLEFT}_{YBOTTOM}")        # see:  https://gis.stackexchange.com/a/441222/147911

  # Create a new catalog with overlapping chunks
  ctg_chunked <- catalog_retile(las_catalog)

  return(ctg_chunked@data)
}

#............................................................................... save a file as a .laz file

save_file_as_laz = function(input_file,
                            output_laz_file,
                            threads,
                            EPSG = 32633) {

  data = data.table::fread(input_file, stringsAsFactors = F, header = T, nThread = threads)
  data = data[, c("x", "y", "z")]
  colnames(data) = c('X', 'Y', 'Z')

  #................................................................................................ verify visually
  # rgl::plot3d(data$X, data$Y, data$Z, col = lidR::height.colors(50)[cut(data$Z, breaks = 50)],
  #             xlab = "X", ylab = "Y", zlab = "Z",
  #             main = "3D LiDAR Point Cloud")
  #................................................................................................

  # convert to a LAS object by using a header without the coordinate reference
  HEADER = lidR::LASheader(data)
  ctg = lidR::LAS(data = data, header = HEADER)
  projection(ctg) <- EPSG

  #................................................. verify visually
  # EXT = sf::st_bbox(ext(ctg), crs = EPSG) |>
  #   sf::st_as_sfc()
  # mapview::mapview(EXT)
  #.................................................

  lidR::writeLAS(las = ctg, file = output_laz_file, index = TRUE)
}

#............................................................................... create overlapping polygons using the data.table R package

datatable_grid = function(data_inp,
                          dist = 0.05,
                          n_grid = 20,
                          verify_visually = FALSE) {

  # Calculate the ranges of x and y
  x_range <- range(data_inp$x, na.rm = TRUE)
  y_range <- range(data_inp$y, na.rm = TRUE)

  # Grid cell width and height
  x_step <- (x_range[2] - x_range[1]) / n_grid
  y_step <- (y_range[2] - y_range[1]) / n_grid

  # Create a grid with 15x15 cells
  grid <- CJ(
    x_bin = seq(1, n_grid, by = 1),
    y_bin = seq(1, n_grid, by = 1)
  )

  # Add bounding box information to each grid cell (with a buffer of 'dist')
  grid[, `:=`(
    xmin = x_range[1] + (x_bin - 1) * x_step - dist,
    xmax = x_range[1] + x_bin * x_step + dist,
    ymin = y_range[1] + (y_bin - 1) * y_step - dist,
    ymax = y_range[1] + y_bin * y_step + dist
  )]

  # Assign each point to the nearest grid cell using findInterval
  data_inp[, `:=`(
    x_bin = findInterval(x, seq(x_range[1], x_range[2], length.out = n_grid + 1), rightmost.closed = TRUE),
    y_bin = findInterval(y, seq(y_range[1], y_range[2], length.out = n_grid + 1), rightmost.closed = TRUE)
  )]

  # Create data_overlaps
  data_overlaps <- data_inp[, .(x_bin = rep(x_bin, each = 9),
                                y_bin = rep(y_bin, each = 9),
                                x = rep(x, each = 9),
                                y = rep(y, each = 9),
                                z = rep(z, each = 9)),
                            by = .I][
                              , `:=`(
                                x_bin_neighbor = rep(c(-1,0,1), times = .N / 9, each = 3) + x_bin,
                                y_bin_neighbor = rep(c(-1,0,1), times = .N / 9) + y_bin
                              )][
                                x_bin_neighbor >= 1 & x_bin_neighbor <= n_grid &
                                  y_bin_neighbor >= 1 & y_bin_neighbor <= n_grid
                              ][
                                grid, on = .(x_bin_neighbor = x_bin, y_bin_neighbor = y_bin),
                                nomatch = 0
                              ][
                                x >= xmin & x <= xmax & y >= ymin & y <= ymax,
                                .(I, x, y, z, x_bin = x_bin_neighbor, y_bin = y_bin_neighbor, xmin, xmax, ymin, ymax)
                              ]

  # remove duplicates
  data_overlaps_xy = unique(data_overlaps, use.key = FALSE)

  # verify visually
  if (verify_visually) {

    idx_keep = which(!duplicated(data_overlaps[, c('xmin', 'xmax', 'ymin', 'ymax')]))
    grid_viz = data_overlaps_xy[idx_keep, , drop = F]

    res_sf = lapply(1:nrow(grid_viz), function(y) {
      x = grid_viz[y, , drop = F]
      bbx_x = sf::st_bbox(c(xmin = x$xmin, xmax = x$xmax, ymax = x$ymax, ymin = x$ymin)) |>
        sf::st_as_sfc() |>
        sf::st_as_sf() |>
        data.table::as.data.table()

      bbx_x
    }) |>
      data.table::rbindlist() |>
      sf::st_as_sf()

    print(mapview::mapview(res_sf))
  }

  # create frequency table for the computed grid
  N_per_grid_cell = data_overlaps_xy[, .(count = .N), by = c('x_bin', 'y_bin')]
  N_per_grid_cell = N_per_grid_cell[order(N_per_grid_cell$count, decreasing = T), ]

  return(list(data = data_overlaps_xy,
              freq_tbl = N_per_grid_cell,
              total_points = sum(N_per_grid_cell$count)))
}


