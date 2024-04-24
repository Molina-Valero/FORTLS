


species.classification <- function(data, tree.tls, dist = 0.05, median = NULL){

  data <- data[, c("x", "y", "z")]

  sp <- data.frame(tree = as.integer(), x = as.numeric(), y = as.numeric(), z = as.numeric())

  for (i in tree.tls$tree) {

  tree <- tree.tls[tree.tls$tree == i, ]

  dat <- data[data$x < tree$x + tree$dbh / 100 &
              data$x > tree$x - tree$dbh / 100 &
              data$y < tree$y + tree$dbh / 100 &
              data$y > tree$y - tree$dbh / 100, ]

  # dat <- dat[dat$z <= 1.65 & dat$z >= 0.95, ]
  dat <- dat[dat$z >= 0.95, ]

  # dat$x <- dat$x - tree$x
  # dat$y <- dat$y - tree$y

  # dat$phi_2<- atan2(dat$y, dat$x)
  # dat$phi_2 <- ifelse(dat$phi_2 < 0, dat$phi_2 + (2 * pi), dat$phi_2)

  # dat$col <- ifelse(dat$phi_2 < pi/2 & dat$phi_2 > 0, 1,
  #                    ifelse(dat$phi_2 < pi & dat$phi_2 > pi/2, 2,
  #                           ifelse(dat$phi_2 < pi*3/2 & dat$phi_2 > pi, 3, 4)))

  # plot(data$x, data$y, asp = 1, col = data$col)

  dat$tree <- i

  sp <- rbind(sp, dat[, c("tree", "x", "y", "z")])

  }

  sp$point <- 1:nrow(sp)

  variables <- geometric.features(sp, dist)
  # variables <- variables[dat$z <= 1.6 & dat$z >= 1, ]
  variables <- variables[dat$z > 1, ]
  variables <- merge(variables, sp[, c("point", "tree")])

  variables <- variables[variables$number_neighbors > 2, ]

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

  if(!is.null(median)){
  tree.variables <- data.frame(tree = tapply(variables$tree, variables$tree, median),
                               first_eigenvalue = tapply(variables$first_eigenvalue, variables$tree, median, na.rm = TRUE),
                               second_eigenvalue = tapply(variables$second_eigenvalue, variables$tree, median, na.rm = TRUE),
                               third_eigenvalue = tapply(variables$third_eigenvalue, variables$tree, median, na.rm = TRUE),
                               sum_eigenvalues = tapply(variables$sum_eigenvalues, variables$tree, median, na.rm = TRUE),
                               omnivariance = tapply(variables$omnivariance, variables$tree, median, na.rm = TRUE),
                               eigenentropy = tapply(variables$eigenentropy, variables$tree, median, na.rm = TRUE),
                               PCA1 = tapply(variables$PCA1, variables$tree, median, na.rm = TRUE),
                               PCA2 = tapply(variables$PCA2, variables$tree, median, na.rm = TRUE),
                               anisotropy = tapply(variables$anisotropy, variables$tree, median, na.rm = TRUE),
                               planarity = tapply(variables$planarity, variables$tree, median, na.rm = TRUE),
                               linearity = tapply(variables$linearity, variables$tree, median, na.rm = TRUE),
                               surface_variation = tapply(variables$surface_variation, variables$tree, median, na.rm = TRUE),
                               sphericity = tapply(variables$sphericity, variables$tree, median, na.rm = TRUE),
                               verticality = tapply(variables$verticality, variables$tree, median, na.rm = TRUE),
                               number_neighbors = tapply(variables$number_neighbors, variables$tree, median, na.rm = TRUE),
                               surface_density = tapply(variables$surface_density, variables$tree, median, na.rm = TRUE),
                               volume_density = tapply(variables$volume_density, variables$tree, median, na.rm = TRUE))
  }

  return(tree.variables)

}
