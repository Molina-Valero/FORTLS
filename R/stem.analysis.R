

.straightness <- function(data, stem.range = NULL){

  n <- data$tree[1]

  if(!is.null(stem.range) & nrow(data) > 2){

    data <- data[data$hi >= stem.range[1] & data$hi <= stem.range[2], ]

  }


  if(nrow(data) < 3){

    SS.max <- NA
    sinuosity <- NA
    lean <- NA

  } else {


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
