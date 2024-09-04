


species.classification <- function(data, tree.tls, dist = 0.1, type = 3,
                                   save.result = TRUE, dir.result = NULL){

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()

  slice <- 0.5

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

  if(nrow(dat) < 1){next}

  }


  if(type == 1){dat <- dat[dat$z <= 1.65 & dat$z >= 0.95, ]}

  if(type == 2){dat <- dat[dat$z >= 0.95, ]}

  if(type == 3){

    dat <- as.data.frame(sf::st_coordinates(voro[voro$tree == i, ]))
    if(nrow(dat) < 1){next}

    colnames(dat) <- c("x", "y", "z")


    kk <- quantile(dat$z, probs = c(0.25,0.5,0.75))

    dat <- dat[dat$z >= kk[1] - slice - dist &
               dat$z <= kk[1] + slice + dist |
               dat$z >= kk[2] - slice - dist &
               dat$z <= kk[2] + slice + dist |
               dat$z >= kk[3] - slice - dist &
               dat$z <= kk[3] + slice + dist, ]

    if(nrow(dat) < 1){next}

    }


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

    variables_25 <- variables[variables$z > kk[1] - slice &
                              variables$z < kk[1] + slice, ]

    variables_50 <- variables[variables$z > kk[2] - slice &
                              variables$z < kk[2] + slice, ]

    variables_75 <- variables[variables$z > kk[3] - slice &
                              variables$z < kk[3] + slice, ]

    }


  if(type != 3){

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
  }

  if(type == 3){

    variables_25 <- merge(variables_25, sp[, c("point", "tree")], by = "point", all = FALSE)

    variables_25 <- variables_25[variables_25$number_neighbors > 2, ]

    variables_25 <- variables_25[!is.na(variables_25$point), ]

    tree.variables_25 <- data.frame(tree = tapply(variables_25$tree, variables_25$tree, mean),
                                    first_eigenvalue = tapply(variables_25$first_eigenvalue, variables_25$tree, mean, na.rm = TRUE),
                                    second_eigenvalue = tapply(variables_25$second_eigenvalue, variables_25$tree, mean, na.rm = TRUE),
                                    third_eigenvalue = tapply(variables_25$third_eigenvalue, variables_25$tree, mean, na.rm = TRUE),
                                    sum_eigenvalues = tapply(variables_25$sum_eigenvalues, variables_25$tree, mean, na.rm = TRUE),
                                    omnivariance = tapply(variables_25$omnivariance, variables_25$tree, mean, na.rm = TRUE),
                                    eigenentropy = tapply(variables_25$eigenentropy, variables_25$tree, mean, na.rm = TRUE),
                                    PCA1 = tapply(variables_25$PCA1, variables_25$tree, mean, na.rm = TRUE),
                                    PCA2 = tapply(variables_25$PCA2, variables_25$tree, mean, na.rm = TRUE),
                                    anisotropy = tapply(variables_25$anisotropy, variables_25$tree, mean, na.rm = TRUE),
                                    planarity = tapply(variables_25$planarity, variables_25$tree, mean, na.rm = TRUE),
                                    linearity = tapply(variables_25$linearity, variables_25$tree, mean, na.rm = TRUE),
                                    surface_variation = tapply(variables_25$surface_variation, variables_25$tree, mean, na.rm = TRUE),
                                    sphericity = tapply(variables_25$sphericity, variables_25$tree, mean, na.rm = TRUE),
                                    verticality = tapply(variables_25$verticality, variables_25$tree, mean, na.rm = TRUE),
                                    number_neighbors = tapply(variables_25$number_neighbors, variables_25$tree, mean, na.rm = TRUE),
                                    surface_density = tapply(variables_25$surface_density, variables_25$tree, mean, na.rm = TRUE),
                                    volume_density = tapply(variables_25$volume_density, variables_25$tree, mean, na.rm = TRUE))



    variables_50 <- merge(variables_50, sp[, c("point", "tree")], by = "point", all = FALSE)

    variables_50 <- variables_50[variables_50$number_neighbors > 2, ]

    variables_50 <- variables_50[!is.na(variables_50$point), ]

    tree.variables_50 <- data.frame(tree = tapply(variables_50$tree, variables_50$tree, mean),
                                    first_eigenvalue = tapply(variables_50$first_eigenvalue, variables_50$tree, mean, na.rm = TRUE),
                                    second_eigenvalue = tapply(variables_50$second_eigenvalue, variables_50$tree, mean, na.rm = TRUE),
                                    third_eigenvalue = tapply(variables_50$third_eigenvalue, variables_50$tree, mean, na.rm = TRUE),
                                    sum_eigenvalues = tapply(variables_50$sum_eigenvalues, variables_50$tree, mean, na.rm = TRUE),
                                    omnivariance = tapply(variables_50$omnivariance, variables_50$tree, mean, na.rm = TRUE),
                                    eigenentropy = tapply(variables_50$eigenentropy, variables_50$tree, mean, na.rm = TRUE),
                                    PCA1 = tapply(variables_50$PCA1, variables_50$tree, mean, na.rm = TRUE),
                                    PCA2 = tapply(variables_50$PCA2, variables_50$tree, mean, na.rm = TRUE),
                                    anisotropy = tapply(variables_50$anisotropy, variables_50$tree, mean, na.rm = TRUE),
                                    planarity = tapply(variables_50$planarity, variables_50$tree, mean, na.rm = TRUE),
                                    linearity = tapply(variables_50$linearity, variables_50$tree, mean, na.rm = TRUE),
                                    surface_variation = tapply(variables_50$surface_variation, variables_50$tree, mean, na.rm = TRUE),
                                    sphericity = tapply(variables_50$sphericity, variables_50$tree, mean, na.rm = TRUE),
                                    verticality = tapply(variables_50$verticality, variables_50$tree, mean, na.rm = TRUE),
                                    number_neighbors = tapply(variables_50$number_neighbors, variables_50$tree, mean, na.rm = TRUE),
                                    surface_density = tapply(variables_50$surface_density, variables_50$tree, mean, na.rm = TRUE),
                                    volume_density = tapply(variables_50$volume_density, variables_50$tree, mean, na.rm = TRUE))



    variables_75 <- merge(variables_75, sp[, c("point", "tree")], by = "point", all = FALSE)

    variables_75 <- variables_75[variables_75$number_neighbors > 2, ]

    variables_75 <- variables_75[!is.na(variables_75$point), ]

    tree.variables_75 <- data.frame(tree = tapply(variables_75$tree, variables_75$tree, mean),
                                    first_eigenvalue = tapply(variables_75$first_eigenvalue, variables_75$tree, mean, na.rm = TRUE),
                                    second_eigenvalue = tapply(variables_75$second_eigenvalue, variables_75$tree, mean, na.rm = TRUE),
                                    third_eigenvalue = tapply(variables_75$third_eigenvalue, variables_75$tree, mean, na.rm = TRUE),
                                    sum_eigenvalues = tapply(variables_75$sum_eigenvalues, variables_75$tree, mean, na.rm = TRUE),
                                    omnivariance = tapply(variables_75$omnivariance, variables_75$tree, mean, na.rm = TRUE),
                                    eigenentropy = tapply(variables_75$eigenentropy, variables_75$tree, mean, na.rm = TRUE),
                                    PCA1 = tapply(variables_75$PCA1, variables_75$tree, mean, na.rm = TRUE),
                                    PCA2 = tapply(variables_75$PCA2, variables_75$tree, mean, na.rm = TRUE),
                                    anisotropy = tapply(variables_75$anisotropy, variables_75$tree, mean, na.rm = TRUE),
                                    planarity = tapply(variables_75$planarity, variables_75$tree, mean, na.rm = TRUE),
                                    linearity = tapply(variables_75$linearity, variables_75$tree, mean, na.rm = TRUE),
                                    surface_variation = tapply(variables_75$surface_variation, variables_75$tree, mean, na.rm = TRUE),
                                    sphericity = tapply(variables_75$sphericity, variables_75$tree, mean, na.rm = TRUE),
                                    verticality = tapply(variables_75$verticality, variables_75$tree, mean, na.rm = TRUE),
                                    number_neighbors = tapply(variables_75$number_neighbors, variables_75$tree, mean, na.rm = TRUE),
                                    surface_density = tapply(variables_75$surface_density, variables_75$tree, mean, na.rm = TRUE),
                                    volume_density = tapply(variables_75$volume_density, variables_75$tree, mean, na.rm = TRUE))



    tree.variables <- merge(tree.variables_25, tree.variables_50, by = "tree", all = FALSE)
    tree.variables <- merge(tree.variables, tree.variables_75, by = "tree", all = FALSE)

    tree.variables[is.na(tree.variables)] <- 0

    tree.variables$first_eigenvalue <- (tree.variables$first_eigenvalue + tree.variables$first_eigenvalue.x + tree.variables$first_eigenvalue.y) / 3
    tree.variables$second_eigenvalue <- (tree.variables$second_eigenvalue + tree.variables$second_eigenvalue.x + tree.variables$second_eigenvalue.y) / 3
    tree.variables$third_eigenvalue <- (tree.variables$third_eigenvalue + tree.variables$third_eigenvalue.x + tree.variables$third_eigenvalue.y) / 3
    tree.variables$sum_eigenvalues <- (tree.variables$sum_eigenvalues + tree.variables$sum_eigenvalues.x + tree.variables$sum_eigenvalues.y) / 3
    tree.variables$omnivariance <- (tree.variables$omnivariance + tree.variables$omnivariance.x + tree.variables$omnivariance.y) / 3
    tree.variables$eigenentropy <- (tree.variables$eigenentropy + tree.variables$eigenentropy.x + tree.variables$eigenentropy.y) / 3
    tree.variables$PCA1 <- (tree.variables$PCA1 + tree.variables$PCA1.x + tree.variables$PCA1.y) / 3
    tree.variables$PCA2 <- (tree.variables$PCA2 + tree.variables$PCA2.x + tree.variables$PCA2.y) / 3
    tree.variables$anisotropy <- (tree.variables$anisotropy + tree.variables$anisotropy.x + tree.variables$anisotropy.y) / 3
    tree.variables$planarity <- (tree.variables$planarity + tree.variables$planarity.x + tree.variables$planarity.y) / 3
    tree.variables$linearity <- (tree.variables$linearity + tree.variables$linearity.x + tree.variables$linearity.y) / 3
    tree.variables$surface_variation <- (tree.variables$surface_variation + tree.variables$surface_variation.x + tree.variables$surface_variation.y) / 3
    tree.variables$sphericity <- (tree.variables$sphericity + tree.variables$sphericity.x + tree.variables$sphericity.y) / 3
    tree.variables$verticality <- (tree.variables$verticality + tree.variables$verticality.x + tree.variables$verticality.y) / 3
    tree.variables$number_neighbors <- (tree.variables$number_neighbors + tree.variables$number_neighbors.x + tree.variables$number_neighbors.y) / 3
    tree.variables$surface_density <- (tree.variables$surface_density + tree.variables$surface_density.x + tree.variables$surface_density.y) / 3
    tree.variables$volume_density <- (tree.variables$volume_density + tree.variables$volume_density.x + tree.variables$volume_density.y) / 3

    tree.variables <- tree.variables[, names(tree.variables_25)]

    }

  # return(tree.variables)

  sp <- data.frame(sp.id = c(0, 1, 2, 3),
                   sp = c("Eucaliptus nitens", "Fagus sylvatica", "Pinus nigra", "Pinus pinea"))

  tree.variables$sp.id <- random.forest.sp(tree.variables)
  tree.tls <- merge(tree.tls, tree.variables[, c("tree", "sp.id")], by = "tree")
  tree.tls <- merge(tree.tls, sp, by = "sp.id")
  tree.tls <- tree.tls[, !(names(tree.tls) %in% c("sp.id"))]


  if(isTRUE(save.result)){

    utils::write.csv(tree.tls,
                     file = file.path(dir.result, "tree.tls.csv"),
                     row.names = FALSE)
  }


  return(tree.tls)


}
