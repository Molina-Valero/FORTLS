


species.classification <- function(data, tree.tls, dist = 0.1, type = 3){

  data <- data[, c("x", "y", "z")]

  sp <- data.frame(tree = as.integer(), x = as.numeric(), y = as.numeric(), z = as.numeric())

  if(type == 3){

  voro <- sf::st_as_sf(data, coords = c("x", "y", "z"))
  voronoi <- sf::st_as_sf(tree.tls, coords = c("x", "y"))
  voronoi <- sf::st_collection_extract(sf::st_voronoi(do.call(c, sf::st_geometry(voronoi))))
  voro$tree <- unlist(sf::st_intersects(voro, voronoi))

  }


  for (i in tree.tls$tree) {

  tree <- tree.tls[tree.tls$tree == i, ]

  if(type != 3){
  dat <- data[data$x < tree$x + tree$dbh / 100 &
              data$x > tree$x - tree$dbh / 100 &
              data$y < tree$y + tree$dbh / 100 &
              data$y > tree$y - tree$dbh / 100, ]
  }


  if(nrow(dat) < 1){next}

  if(type == 1){dat <- dat[dat$z <= 1.65 & dat$z >= 0.95, ]}

  if(type == 2){dat <- dat[dat$z >= 0.95, ]}

  if(type == 3){

    dat <- as.data.frame(sf::st_coordinates(voro[voro$tree == i, ]))
    if(nrow(dat) < 1){next}

    colnames(dat) <- c("x", "y", "z")


    kk <- quantile(dat$z, probs = c(0.25,0.5,0.75))

    dat <- dat[dat$z >= kk[1] - 1 - dist &
               dat$z <= kk[1] + 1 + dist |
               dat$z >= kk[2] - 1 - dist &
               dat$z <= kk[2] + 1 + dist |
               dat$z >= kk[3] - 1 - dist &
               dat$z <= kk[3] + 1 + dist, ]
    }

  if(nrow(dat) < 1){next}

  dat$tree <- i

  sp <- rbind(sp, dat[, c("tree", "x", "y", "z")])

  }

  sp$point <- 1:nrow(sp)

  variables <- geometric.features(sp, dist)

  if(type == 1){
  variables <- variables[variables$z < 1.6 & variables$z > 1, ]}

  if(type == 2){
  variables <- variables[variables$z > 1, ]}

  if(type == 3){
    variables <- variables[variables$z > kk[1] - 1 &
                           variables$z < kk[1] + 1 |
                           variables$z > kk[2] - 1 &
                           variables$z < kk[2] + 1 |
                           variables$z > kk[3] - 1 &
                           variables$z < kk[3] + 1, ]
    }


  variables <- merge(variables, sp[, c("point", "tree")], by = "point", all = FALSE)

  variables <- variables[variables$number_neighbors > 2, ]

  variables <- variables[!is.na(variables$point), ]

  tree.variables <- data.frame(tree = tapply(variables$tree, variables$tree, mean),
                               first_eigenvalue = tapply(variables$first_eigenvalue, variables$tree, mean, na.rm = TRUE),
                               second_eigenvalue = tapply(variables$second_eigenvalue, variables$tree, mean, na.rm = TRUE),
                               third_eigenvalue = tapply(variables$third_eigenvalue, variables$tree, mean, na.rm = TRUE),
                               sum_eigenvalues = tapply(variables$sum_eigenvalues, variables$tree, mean, na.rm = TRUE),
                               omnivariance = tapply(variables$omnivariance, variables$tree, mean, na.rm = TRUE),
                               eigenentropy = tapply(variables$eigenentropy, variables$tree, mean, na.rm = TRUE),
                               PCA1 = tapply(variables$PCA1, variables$tree, mean, na.rm = TRUE),
                               PCA2 = tapply(variables$PCA2, variables$tree, mean, na.rm = TRUE),
                               anisotropy = tapply(variables$anisotropy, variables$tree, mean, na.rm = TRUE),
                               planarity = tapply(variables$planarity, variables$tree, mean, na.rm = TRUE),
                               linearity = tapply(variables$linearity, variables$tree, mean, na.rm = TRUE),
                               surface_variation = tapply(variables$surface_variation, variables$tree, mean, na.rm = TRUE),
                               sphericity = tapply(variables$sphericity, variables$tree, mean, na.rm = TRUE),
                               verticality = tapply(variables$verticality, variables$tree, mean, na.rm = TRUE),
                               number_neighbors = tapply(variables$number_neighbors, variables$tree, mean, na.rm = TRUE),
                               surface_density = tapply(variables$surface_density, variables$tree, mean, na.rm = TRUE),
                               volume_density = tapply(variables$volume_density, variables$tree, mean, na.rm = TRUE))

  return(tree.variables)

}
