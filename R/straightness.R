
# Calculate overlap area between two circles
.circle_overlap_area <- function(x1, y1, r1, x2, y2, r2) {
  # Distance between centers
  d <- sqrt((x2 - x1)^2 + (y2 - y1)^2)

  # Case 1: No overlap
  if (d >= r1 + r2) return(0)

  # Case 2: One circle contains the other
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)

  # Case 3: Partial overlap (lens intersection)
  angle1 <- 2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  angle2 <- 2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))

  area1 <- 0.5 * r1^2 * (angle1 - sin(angle1))
  area2 <- 0.5 * r2^2 * (angle2 - sin(angle2))

  return(area1 + area2)
}


# Calculate union area of two circles
.circle_union_area <- function(x1, y1, r1, x2, y2, r2) {
  area1 <- pi * r1^2
  area2 <- pi * r2^2
  overlap <- .circle_overlap_area(x1, y1, r1, x2, y2, r2)
  return(area1 + area2 - overlap)
}


# Helper Functions ----


# Calculate overlap ratio between two consecutive measurements
.calculate_overlap_ratio <- function(row1, row2) {
  overlap <- .circle_overlap_area(
    row1$x, row1$y, row1$dhi / 2,
    row2$x, row2$y, row2$dhi / 2
  )
  union <- .circle_union_area(
    row1$x, row1$y, row1$dhi / 2,
    row2$x, row2$y, row2$dhi / 2
  )
  return(overlap / union)
}


# Calculate combined relative RMSE for diameter and position models
.calculate_combined_rmse <- function(data) {
  fit_d <- lm(dhi ~ hi, data = data)
  fit_x <- lm(x ~ hi, data = data)
  fit_y <- lm(y ~ hi, data = data)

  rel_rmse_d <- summary(fit_d)$sigma / mean(predict(fit_d))
  rel_rmse_x <- summary(fit_x)$sigma / mean(predict(fit_x))
  rel_rmse_y <- summary(fit_y)$sigma / mean(predict(fit_y))

  # rel_rmse_d <- summary(fit_d)$sigma / mean(data$dhi, na.rm = TRUE)
  # rel_rmse_x <- summary(fit_x)$sigma / mean(data$x, na.rm = TRUE)
  # rel_rmse_y <- summary(fit_y)$sigma / mean(data$y, na.rm = TRUE)

  return(mean(c(rel_rmse_d, rel_rmse_x, rel_rmse_y), na.rm = TRUE))
}


# Predict diameter and position at given height
.predict_at_height <- function(data, target_height) {
  fit_d <- lm(dhi ~ hi, data = data)
  fit_x <- lm(x ~ hi, data = data)
  fit_y <- lm(y ~ hi, data = data)

  return(c(
    dhi = coef(fit_d)[1] + coef(fit_d)[2] * target_height,
    x = coef(fit_x)[1] + coef(fit_x)[2] * target_height,
    y = coef(fit_y)[1] + coef(fit_y)[2] * target_height
  ))
}


# Main Function ----

# Refine tree diameter measurements
#
# Three-step process:
# 1. Find the most reliable segment (lowest RMSE)
# 2. Correct measurements above the segment
# 3. Correct measurements below the segment
#
# @param data Data frame with: tree, hi (height), dhi (diameter), x, y
# @param window_size Number of consecutive points for fitting (default: 10)
# @param overlap_threshold Min overlap ratio for valid measurements (default: 0.75)
# @param make_plots Create diagnostic plots (default: TRUE)
# @return Corrected data frame

.refine_diameters <- function(data,
                              window_size = 10,
                              overlap_threshold = 0.75,
                              make_plots = TRUE) {

  # Validation
  required_cols <- c("tree", "hi", "dhi", "x", "y")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain: ", paste(required_cols, collapse = ", "))
  }
  if (nrow(data) < window_size) {
    warning("Insufficient data. Returning original.")
    return(data)
  }

  corrected <- data

  # Step 1: Find best segment
  if (make_plots) {
    plot(data$dhi ~ data$hi,
         main = paste("Tree:", data$tree[1]),
         xlab = "Height", ylab = "Diameter",
         pch = 1, col = "gray60")
  }

  n_windows <- nrow(data) - window_size + 1
  rmse_vals <- sapply(1:n_windows, function(i) {
    .calculate_combined_rmse(data[i:(i + window_size - 1), ])
  })

  best_idx <- which.min(rmse_vals)
  best_segment <- data[best_idx:(best_idx + window_size - 1), ]

  if (make_plots) {
    points(best_segment$dhi ~ best_segment$hi, pch = 19, col = "darkgreen")
  }

  # Step 2: Forward pass (upward correction)
  last_good <- best_idx + window_size - 1
  if (last_good < nrow(data)) {
    for (i in last_good:(nrow(data) - 1)) {
      ratio <- .calculate_overlap_ratio(corrected[i, ], corrected[i + 1, ])

      if (ratio < overlap_threshold) {
        start <- max(1, i - window_size + 1)
        pred <- .predict_at_height(corrected[start:i, ], corrected[i + 1, "hi"])
        corrected[i + 1, c("dhi", "x", "y")] <- pred
      }
    }
  }

  # Step 3: Backward pass (downward correction)
  if (best_idx > 1) {
    for (i in best_idx:2) {
      ratio <- .calculate_overlap_ratio(corrected[i - 1, ], corrected[i, ])

      if (ratio < overlap_threshold) {
        end <- min(nrow(data), i + window_size - 1)
        pred <- .predict_at_height(corrected[i:end, ], corrected[i - 1, "hi"])
        corrected[i - 1, c("dhi", "x", "y")] <- pred
      }
    }
  }

  if (make_plots) {
    points(corrected$dhi ~ corrected$hi, col = "blue", pch = 19, cex = 0.5)
    legend("topright",
           legend = c("Original", "Best segment", "Corrected"),
           pch = c(1, 19, 19),
           col = c("gray60", "darkgreen", "blue"))
  }

  return(corrected)
}





# setwd("C:/pruebas")
#
# tree <- read.csv("stem.curve.csv")
#
# # Example
# for (i in unique(tree$tree)) {
#
# data <- tree[tree$tree == i, ]
#
# kk <- .refine_diameters(data)
#
# }





.straightness <- function(data, stem.range = NULL){


  n <- data$tree[1]


  # Retaining sections within the stem range defined in the arguments

  if(!is.null(stem.range) & nrow(data) > 2){

    data <- data[data$hi >= stem.range[1] & data$hi <= stem.range[2], ]

  }


  # SS.max, sinuosity and lean are not estimated if there are less than 3 sections

  if(nrow(data) < 3){

    SS.max <- NA
    sinuosity <- NA
    lean <- NA

  } else {

    # Filtering outliers

    data$dif.sec <- c(abs(diff(data$hi)), 0)
    data$dif <- c(diff(data$dhi),
                  data$dhi[nrow(data)] - data$dhi[nrow(data) - 1])
    data$dif <- ifelse(data$dif.sec > 1, data$dif / data$dif.sec, data$dif)

    while (suppressWarnings(max(abs(data$dif))) > 10) {

      data$filter <- ifelse(data$dif > 10 | data$dif < -10, 0, 1)
      data <- data[data$filter == 1, ]
      if(nrow(data) < 2) {next}

      data$dif.sec <- c(abs(diff(data$hi)), 0)
      data$dif <- c(diff(data$dhi),
                    data$dhi[nrow(data)] - data$dhi[nrow(data) - 1])
      data$dif <- ifelse(data$dif.sec > 1, data$dif / data$dif.sec, data$dif)

    }


  # SS.max, sinuosity and lean are not estimated if there are less than 3 sections

  if(nrow(data) < 3){

  SS.max <- NA
  sinuosity <- NA
  lean <- NA

  out <- data.frame(tree = n, SS.max = SS.max, sinuosity = sinuosity, lean = lean)

  return(out)

  }


  mod.x <- lm(data$x[c(1,nrow(data))]~data$hi[c(1,nrow(data))])
  mod.y <- lm(data$y[c(1,nrow(data))]~data$hi[c(1,nrow(data))])

  data$sagita.x <- coef(mod.x)[1] + coef(mod.x)[2] * data$hi
  data$sagita.y <- coef(mod.y)[1] + coef(mod.y)[2] * data$hi
  data$sagita <- sqrt((data$sagita.x-data$x)^2+(data$sagita.y-data$y)^2)

  S.max <- max(data$sagita)

  # h.range <- max(data$hi)-min(data$hi)
  # h.range <- sqrt((data$x[1]-data$x[nrow(data)]) ^ 2 + (data$y[1]-data$y[nrow(data)]) ^ 2 + (data$hi[1]-data$hi[nrow(data)]) ^ 2)
  h.range <- sqrt((data$x[nrow(data)]-data$x[1]) ^ 2 + (data$y[nrow(data)]-data$y[1]) ^ 2 + (data$hi[nrow(data)]-data$hi[1]) ^ 2)

  R <- (S.max ^ 2 + (h.range / 2) ^ 2) / (2 * S.max)

  if((R ^ 2 - (max(data$hi) / 2) ^ 2) > 0){
    SS.max <- R - sqrt(R ^ 2 - (max(data$hi) / 2) ^ 2)
    } else {SS.max <- S.max}


  # Sinuosity

  data$x.2 <- c(data$x[2:nrow(data)], data$x[nrow(data)])
  data$y.2 <- c(data$y[2:nrow(data)], data$y[nrow(data)])
  data$hi.2 <- c(data$hi[2:nrow(data)], data$hi[nrow(data)])

  data$l <- sqrt((data$x.2-data$x) ^ 2 + (data$y.2-data$y) ^ 2 + (data$hi.2-data$hi) ^ 2)

  sinuosity <- sum(data$l) / h.range

  # Lean

  pca <- stats::princomp(data[, c("x", "y", "hi")])
  XY <- sqrt(pca$loadings[, 1][1] ^ 2 + pca$loadings[, 1][2] ^ 2)

  if(XY == 0){XY <- 2e-16}
  lean <- atan(abs(pca$loadings[, 1][3]) / XY) * (360 / (2 * pi))

  }

  out <- data.frame(tree = data$tree[1], SS.max = SS.max, sinuosity = sinuosity, lean = lean)

  return(out)

}



# Example
# for (i in unique(tree$tree)) {
#
#   data <- tree[tree$tree == i, ]
#
#   kk <- .refine_diameters(data)
#
#   kk1 <- .straightness(kk)
#   kk2 <- .straightness(data)
#
# }


