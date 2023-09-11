

.straightness <- function(data, stem.range = NULL){

  if(!is.null(stem.range) & nrow(data) > 2){

    data <- data[data$hi >= stem.range[1] & data$hi <= stem.range[2], ]

  }


  if(nrow(data) < 3){

    SS.max <- NA
    sinuosity <- NA

  } else {


  mod.x <- lm(data$x[c(1,nrow(data))]~data$hi[c(1,nrow(data))])
  mod.y <- lm(data$y[c(1,nrow(data))]~data$hi[c(1,nrow(data))])


  data$sagita.x <- coef(mod.x)[1] + coef(mod.x)[2] * data$hi
  data$sagita.y <- coef(mod.y)[1] + coef(mod.y)[2] * data$hi
  data$sagita <- sqrt((data$sagita.x-data$x)^2+(data$sagita.y-data$y)^2)

  S.max <- max(data$sagita)

  # h.range <- max(data$hi)-min(data$hi)
  h.range <- sqrt((data$x[1]-data$x[nrow(data)]) ^ 2 + (data$y[1]-data$y[nrow(data)]) ^ 2 + (data$hi[1]-data$hi[nrow(data)]) ^ 2)

  R <- (S.max ^ 2 + (h.range / 2) ^ 2) / (2 * S.max)

  if((R ^ 2 - (max(data$hi) / 2) ^ 2) > 0){
    SS.max <- R - sqrt(R ^ 2 - (max(data$hi) / 2) ^ 2)
    } else {SS.max <- S.max}


  data$x.2 <- c(data$x[2:nrow(data)], data$x[nrow(data)])
  data$y.2 <- c(data$y[2:nrow(data)], data$y[nrow(data)])
  data$hi.2 <- c(data$hi[2:nrow(data)], data$hi[nrow(data)])

  data$l <- sqrt((data$x-data$x.2) ^ 2 + (data$y-data$y.2) ^ 2 + (data$hi-data$hi.2) ^ 2)

  sinuosity <- sum(data$l) / h.range

  }

  out <- data.frame(tree = data$tree[1], SS.max = SS.max, sinuosity = sinuosity)

  return(out)

}
