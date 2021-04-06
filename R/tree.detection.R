tree.detection <- function(data, dbh.min = 7.5, dbh.max = 200, ncr.threshold = 0.1, tls.resolution = list(), breaks=c(1.0, 1.3, 1.6), plot.attributes = NULL,
                           dir.result = NULL, save.result = TRUE){

  # Obtaining working directory for saving files
  if(is.null(dir.result))
    dir.result <- getwd()

  # Converting arguments to International System of Units
  # Arguments of the forest inventory

  .dbh.min <- dbh.min / 100
  .dbh.max <- dbh.max / 100

  # Arguments of the TLS precision
  # If resolution is defined by points distance at a certain distance from TLS (mm/m):
  .point.dist <- tls.resolution$point.dist / 1000
  .tls.dist <- tls.resolution$tls.dist

  # If resolution is defined by angles (?):
  .vertical.angle <- tls.resolution$vertical.angle * 2 * pi / 360
  .horizontal.angle <- tls.resolution$horizontal.angle * 2 * pi / 360

  # Ahora en funci?n de los par?metros de la resoluci?n del TLS que tengamos, hay que definir
  # la apertura angular:

  # Angular resolution:
  if(is.null(.point.dist)){

    .alpha.v <- .vertical.angle
    .alpha.h <- .horizontal.angle

  } else {

    .alpha.v <- atan((.point.dist / 2) / .tls.dist) * 2
    .alpha.h <- .alpha.v

  }

  # Generation of homogenized point cloud

  # .ncr.threshold <- .ncr.threshold.double(data)

  .filteraux <- data.frame(cluster = as.numeric(),
                           center.x = as.numeric(), center.y = as.numeric(),
                           center.phi = as.numeric(), center.rho = as.numeric(),
                           center.r = as.numeric(), center.theta = as.numeric(),
                           radius = as.numeric(),
                           num.points = as.numeric(), num.points.hom = as.numeric(),
                           phi.left = as.numeric(), phi.right = as.numeric(),
                           arc.circ = as.numeric(), sec = as.numeric())


  for(cuts in breaks){

    message("Computing section: ", cuts, " m")

    .cut <- data[which(data$z > (cuts-0.1) & data$z < (cuts+0.1)), ]

    .cut <- .ncr.remove.slice.double(.cut)

    .cut <- .cut[which(.cut$ncr < ncr.threshold | is.na(.cut$ncr)), ]

    # Ahora ya si me quedo con la rebanada de 1m +/- 5 cm
    .cut <- .cut[which(.cut$z > (cuts-0.05) & .cut$z < (cuts+0.05)), ]

    # Dbscan parameters
    .eps <- .dbh.min / 2

    # Clustering
    .dbscan <- dbscan::dbscan(.cut[, c("x", "y")], eps = .eps)
    .cut$cluster <- .dbscan$cluster
    .cut <- .cut[which(.cut$cluster > 0), ]

    # Checking if there are clusters
    if(nrow(.cut) < 1)
      next

    .cut$sec <- cuts

    .pb <- progress::progress_bar$new(total = length(unique(.cut$cluster)))

    .filter <- data.frame(cluster = as.numeric(),

                          # Center coordinates
                          center.x = as.numeric(), center.y = as.numeric(),
                          center.phi = as.numeric(), center.rho = as.numeric(),
                          center.r = as.numeric(), center.theta = as.numeric(),

                          # Radius
                          radius = as.numeric(),

                          # Number of points belowing to cluster (craw and after point crooping)
                          num.points = as.numeric(), num.points.hom = as.numeric(),

                          # Phi coordinates of left and right
                          phi.left = as.numeric(), phi.right = as.numeric(),

                          # Circunference arc
                          arc.circ = as.numeric(),

                          # Partial occulion
                          occlusion = as.numeric())


    for(.i in unique(.cut$cluster)){

      .pb$tick()

      # Selecionamos el cluster i
      .dat <- .cut[which(.cut$cluster == .i), ]

      # Generacion de la malla.

      .x.rang <- max(.dat$x) - min(.dat$x)
      .y.rang <- max(.dat$y) - min(.dat$y)
      .phi.rang <- max(.dat$phi) - min(.dat$phi)
      .rho.rang <- max(.dat$rho) - min(.dat$rho)

      # Ahora se calculan las coordendas de los centroides con respecto al TLS:
      .x.cent <- (.x.rang / 2) + min(.dat$x)
      .y.cent <- (.y.rang / 2) + min(.dat$y)
      .phi.cent <- (.phi.rang / 2) + min(.dat$phi)
      .rho.cent <- (.rho.rang / 2) + min(.dat$rho)

      # Se obtiene el ancho de la malla cuadrada que se aplicar? sobre el cl?ster
      .ancho.malla <- (max(.x.rang, .y.rang) / 2) * 1.5

      # Coordenadas m?ximas y m?nimas de la malla
      .xmin <- .x.cent - .ancho.malla
      .ymin <- .y.cent - .ancho.malla
      .xmax <- .x.cent + .ancho.malla
      .ymax <- .y.cent + .ancho.malla

      # Ahora se obtiene una segunda malla basada en las coornadas cil?ndricas phi y rho.
      # Esta se utilizar? para calcular la densidad media de puntos por celda
      .ancho.malla.2 <- (max(.phi.rang, .rho.rang) / 2)

      .phimin <- .phi.cent - .ancho.malla.2
      .rhomin <- .rho.cent - .ancho.malla.2
      .phimax <- .phi.cent + .ancho.malla.2
      .rhomax <- .rho.cent + .ancho.malla.2

      # Filtro
      .x.values <- seq(from = .xmin, to = .xmax, by = 0.03)
      .y.values <- seq(from = .ymin, to = .ymax, by = 0.03)

      .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

      for(.i in 1:length(.x.values)){
        for(.j in 1:length(.y.values)){

          .den <- .dat[which(.dat$x < ((.x.values[.i]) + 0.015) &
                               .dat$x > ((.x.values[.i]) - 0.015) &
                               .dat$y < ((.y.values[.j]) + 0.015) &
                               .dat$y > ((.y.values[.j]) - 0.015)), ]

          # Aquellas celdas con menos de 2 puntos no las tengo en cuenta
          # para luego m?s tarde calcular la densidad media por celda
          .density[.j, .i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

        }

      }


      # Estmiaci?n de la densidad media por celda
      .threeshold <- mean(.density, na.rm = T)

      if(is.nan(.threeshold)){next}

      .density <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))
      .remove <- data.frame(point = as.numeric())

      for(.i in 1:length(.x.values)){
        for(.j in 1:length(.y.values)){

          .den <- .dat[which(.dat$x < ((.x.values[.i]) + 0.015) &
                             .dat$x > ((.x.values[.i]) - 0.015) &
                             .dat$y < ((.y.values[.j]) + 0.015) &
                             .dat$y > ((.y.values[.j]) - 0.015)), ]

      # Aquellas celdas con menos de 2 puntos no las tengo en cuenta
      # para luego m?s tarde calcular la densidad media por celda
      .density[.j, .i] <- ifelse(nrow(.den) < 1, NA, nrow(.den))

      if(nrow(.den) > .threeshold){

      .rem <- data.frame(point = .den$point)
      .remove <- rbind(.remove, .rem)

          }

        }

      }

      .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)

      if(nrow(.dat) < 1){next}

      .n <- 0.7 * (0.1 / (tan(.alpha.h / 2) * mean(.dat$rho) * 2))

      # Valores de las coordenadas phi y rho correspondientes a las
      # intersecciones de la malla:
      .x2.values <- seq(from = .phimin, to = .phimax, by = .alpha.h)
      .y2.values <- seq(from = .rhomin, to = .rhomax, by = 0.04)

      # Matriz donde se va a almacenar el n?mero de puntos por celda
      .density <- matrix(0, ncol = length(.x2.values), nrow = length(.y2.values))

      .remove <- data.frame(point = as.numeric())

      for(.i in 1:length(.x2.values)){
        for(.j in 1:length(.y2.values)){

          .den <- .dat[which(.dat$phi <= ((.x2.values[.i]) + (.alpha.h/2)) &
                             .dat$phi >= ((.x2.values[.i]) - (.alpha.h/2)) &
                             .dat$rho <= ((.y2.values[.j]) + 0.02) &
                             .dat$rho <= ((.y2.values[.j]) - 0.02)), ]

      # Aquellas celdas con menos de 2 puntos no las tengo en cuenta
      # para luego m?s tarde calcular la densidad media por celda
      .density[.j, .i] <- ifelse(nrow(.den) < 2, NA, nrow(.den))

      if(nrow(.den) > 1){

          .rem <- data.frame(point = .den$point)
          .remove <- rbind(.remove, .rem)

        }

      }

    }


     # Estmiaci?n de la densidad media por celda
    .density <- ifelse(is.nan(.density), NA, .density)

    if(is.nan(mean(.density, na.rm = TRUE))){next}

    if(max(.density[which(!is.na(.density))], na.rm = T) < .n){next}

      # Eliminaci?n de las celdas con solo un punto
      .dat <- merge(.dat, .remove, by = "point", all.y = TRUE)

      # Si despu?s de esta operaci?n el .dat no tiene puntos se pasa a la siguiente iteraci?n
      if(nrow(.dat) < 1){next}

      # Estimaci?n del n?mero de puntos tanto de la nube original (.num) como despu?s del
      # proceso de point crooping (.numHom)
      .num.points <- nrow(.dat)
      .num.points.hom <- nrow(.dat[which(.dat$prob.selec == 1), ])

      # Ahora que se ha hecho un peque?o filtrado se calcula el centroide del cl?ster

      .x.values <- seq(from = .xmin, to = .xmax, by = 0.01)
      .y.values <- seq(from = .ymin, to = .ymax, by = 0.01)

      # Se crea una matriz vac?a donde se va a guardar para cada intersecci?n de la malla,
      # la varianza de las distancias entre los puntos y la correspondiente intersecci?n
      .matriz <- matrix(0, ncol = length(.x.values), nrow = length(.y.values))

      for(.i in 1:length(.x.values)){
        for(.j in 1:length(.y.values)){

          .variance <- stats::var(raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.i], .y.values[.j]), lonlat=FALSE))
          .matriz[.j, .i] <- .variance

        }
      }

      # Aquella intersecci?n donde la variznza sea m?nima ser? considerada como el centro
      # de la secci?n
      .a <- which(.matriz == min(.matriz), arr.ind = TRUE)

      .center.x <- .x.values[.a[2]]
      .center.y <- .y.values[.a[1]]

      .center.phi <- atan2(.center.y, .center.x)
      .center.phi <- ifelse(.center.phi < 0, .center.phi + (2 * pi), .center.phi)
      .center.rho <- sqrt(.center.x ^ 2 + .center.y ^ 2)
      .center.r <- sqrt(.dat$sec[1] ^ 2 + .center.rho ^ 2)
      .center.theta <- atan2(.dat$sec[1], .center.rho)

      # Radius value as the mean distance
      .radio <- mean(raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE))


      # Distancias de los puntos al centro
      .dat$dist <- raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE)


      # Center behind tree surface
      if(stats::quantile(.dat$rho, prob = 0.05) > .center.r) {next}

      # Al menos el 95 % de las distancias deber?a tener un valor superior a radio/2
      .dat <- .dat[order(.dat$dist, decreasing = FALSE), ]
      .circle <- if(stats::quantile(.dat$dist, prob = 0.05) < (.radio / 2)) {next}


      # Calculo de la coordenada rho de los extremos de la secci?n

      if((max(.dat$phi) - min(.dat$phi)) < pi){

        .dat.2 <- .dat[order(.dat$phi, decreasing = F), ]

       } else {

         .dat.2 <- .dat
         .dat.2$phi <- ifelse(.dat.2$phi < 1, .dat.2$phi + (2 * pi), .dat.2$phi)
         .dat.2 <- .dat.2[order(.dat.2$phi, decreasing = F), ]

       }

       # Me quedo con el percentil 1 por si hubiese alg?n punto raro
       # Aqu? me faltar?a por ver qu? pasa si el cl?ster est? ubicado en 0 +/- phi
       .pto.left <- stats::quantile(.dat.2$phi, prob = 0.01)
       .rho.left <- mean(.dat.2$rho[which(.dat.2$phi <= .pto.left)])
       .phi.left <- mean(.dat.2$phi[which(.dat.2$phi <= .pto.left)])

        # Me quedo con el percentil 99 por si hubiese alg?n punto raro
        .pto.right <- stats::quantile(.dat.2$phi, prob = 0.99)
        .rho.right <- mean(.dat.2$rho[which(.dat.2$phi >= .pto.right)])
        .phi.right <- mean(.dat.2$phi[which(.dat.2$phi >= .pto.right)])

        # Para los puntos del centro de la secci?n me quedo con aquellos que est?n en la
        # mitad de la apertura angular phi +/- la apertura del TLS .alpha
        .phi.cent <- max(.dat.2$phi) - ((max(.dat.2$phi) - min(.dat.2$phi)) / 2)
        .rho.cent <- mean(.dat.2$rho[which(.dat.2$phi > (.phi.cent - .alpha.h) & .dat.2$phi < (.phi.cent + .alpha.h))])

        if(is.nan(.rho.cent)){next}

        # Ahora se comprueba que las coordenadas rho de los extremos sean mayores que la central
        .arc.circ <- ifelse(.rho.left > .rho.cent & .rho.right > .rho.cent, 1, 0)

        # Conversi?n a las coordendas originales en caso de que el cl?ster estuviera situado
        # en 0 +/- phi
        .phi.left <- ifelse(.phi.left > (2 * pi), .phi.left - (2 * pi), .phi.left)
        .phi.right <- ifelse(.phi.right > (2 * pi), .phi.right - (2 * pi), .phi.right)

        # En caso de no haber detectado el arco de circunferencia completo y por lo
        # no se cumple el criterio anterior, comprobamos si ha habido una oclusi?n parcial.
        # En caso de pertenecer a una secci?n de un arbol, la regularidad sistem?tica con
        # la que se encuentran debe dar correlaciones alt?simas entre su nemeraci?n correlativa
        # y phi cuando estos han sido ordenados con respecto a phi.
        .dat.2$n <- c(1:nrow(.dat.2))
        .cor <- try(stats::cor.test(x = .dat.2$n, y = .dat.2$phi, method = 'pearson'))#podr?a utilizarse el comando cor.
        # .cor <- try(stats::cor.test(x = .dat2$n, y = .dat2$phi, method = 'spearman'))

        # En caso de error se salta a la siguiente iteraci?n
        if(class(.cor) == "try-error"){next} else{

          .occlusion <- .cor[[4]]

        }


        # Coefficient of variation for distances among cluster points and the estimated center
        .cv <- stats::sd(raster::pointDistance(cbind(.dat$x,.dat$y), c(.x.values[.a[2]], .y.values[.a[1]]), lonlat = FALSE)) / .radio

        if(.cv > 0.1){next}

        # Zhang et al., (2019)
        .n.w.ratio <- stats::sd(.dat$z) / sqrt(stats::sd(.dat$x) ^ 2 + stats::sd(.dat$y) ^ 2)

        if(.n.w.ratio > 1.5){next}

        # Results
        .salida <- data.frame(cluster = .dat$cluster[1],

                              center.x = .center.x, center.y = .center.y,
                              center.phi = .center.phi, center.rho = .center.rho,
                              center.r = .center.r, center.theta = .center.theta,

                              radius = .radio,

                              num.points = .num.points, num.points.hom = .num.points.hom,

                              phi.left = .phi.left, phi.right = .phi.right,

                              arc.circ = .arc.circ, occlusion = .occlusion)

              .filter <- rbind(.filter, .salida)

    }


    # Arch of circle or partial arch of circle?
    .filter$tree <- ifelse(.filter$arc.circ == 1, 1,
                           ifelse(.filter$arc.circ == 0 & .filter$occlusion > 0.995, 1, 0))
    .filter <- .filter[which(.filter$tree == 1), ]

    # Dbh maximum and minimum
    .filter$tree <- ifelse(.filter$radius > (.dbh.max / 2) | .filter$radius < (.dbh.min / 2), 0, 1)
    .filter <- subset(.filter, .filter$tree == 1)


    if(nrow(.filter) < 1){

      .filter1.0 <- data.frame(cluster = as.numeric(),
                               center.x = as.numeric(), center.y = as.numeric(),
                               center.phi = as.numeric(), center.rho = as.numeric(),
                               center.r = as.numeric(), center.theta = as.numeric(),
                               radius = as.numeric(),
                               num.points = as.numeric(), num.points.hom = as.numeric(),
                               phi.left = as.numeric(), phi.right = as.numeric(),
                               arc.cir = as.numeric(), sec = as.numeric())

    } else{

      .filter1.0 <- .filter[, c("cluster",
                                "center.x", "center.y", "center.phi", "center.rho", "center.r", "center.theta",
                                "radius", "num.points", "num.points.hom", "phi.left", "phi.right", "arc.circ")]
      .filter1.0$sec <- cuts

    }

    .filteraux<-rbind(.filteraux,.filter1.0)

  }#fin del bucle de los cortes

  .filter<-.filteraux

  # Esto se repetir?a para el resto de secciones!!!
  # Merge of all sections ----

  #.filter <- rbind(.filter1.0, .filter1.3, .filter1.6)

  if(nrow(.filter) < 1) stop("No tree was detected")

  .dbscan <- dbscan::dbscan(.filter[, c("center.x", "center.y")], eps = max(.filter$radius), minPts = 1)
  .filter$cluster <- .dbscan$cluster
  .filter <- .filter[order(.filter$cluster, .filter$sec), ]

  # Calculate of taper coeficient as the slope coeficient of linear regresion:

  .taper <- .filter[, c("cluster", "sec", "radius")]


  .lm <- stats::lm(radius ~ sec, data = .taper)
  .slope <- stats::coef(.lm)[2]


  # Creando una nueva columna ("dif") con la diferencia entre el di?metro normal y la
  # secci?n desde la cual se estima el radio, no importar?a el n?mero de secciones (o cortes)
  # empleados. Siempre habr?a un radio estimado
  .filter$dif <- 1.3 - .filter$sec
  .filter$radio.est <- ifelse(.filter$dif == 0, .filter$radius,
                              .filter$radius + .slope * .filter$dif)


  .radio.est <- tapply(.filter$radio.est, .filter$cluster, mean)

  # Dendrometric variables
  .tree <- data.frame(tree = tapply(.filter$cluster, .filter$cluster, mean, na.rm = TRUE),
                      center.x = tapply(.filter$center.x, .filter$cluster, mean, na.rm = TRUE),
                      center.y = tapply(.filter$center.y, .filter$cluster, mean, na.rm = TRUE),
                      center.phi = tapply(.filter$center.phi, .filter$cluster, mean, na.rm = TRUE),
                      center.rho = tapply(.filter$center.rho, .filter$cluster, mean, na.rm = TRUE),
                      center.r = tapply(.filter$center.r, .filter$cluster, mean, na.rm = TRUE),
                      center.theta = tapply(.filter$center.theta, .filter$cluster, mean,na.rm = TRUE),

                      horizontal.distance = tapply(.filter$center.rho, .filter$cluster, mean, na.rm = TRUE),#l?nea repetida
                      radius = .radio.est,

                      phi.left = tapply(.filter$phi.left, .filter$cluster, mean, na.rm = TRUE),
                      phi.right = tapply(.filter$phi.right, .filter$cluster, mean,na.rm = TRUE),

                      partial.occlusion = tapply(.filter$arc.circ, .filter$cluster, mean, na.rm = TRUE),

                      num.points = tapply(.filter$num.points, .filter$cluster, mean, na.rm = TRUE),
                      num.points.hom = tapply(.filter$num.points.hom, .filter$cluster, mean, na.rm = TRUE))

  # Indicamos si hay ?rboles con oclusiones parciales, que ser?an aquellos en los que
  # ninguna de las secciones ha podido ser identificada como arco de circunferencia (ArcCirc)
  .tree$partial.occlusion <- ifelse(.tree$partial.occlusion == 0, 1, 0)

  .tree$dbh <- .tree$radius * 2

  # Calculate of points belowing to radio unit
  # Ya que voy a hacer una estimaci?n, me quedo con aquellas secciones que han sido
  # completamente visibles (ArcCirc == 1) en la secci?n de 1.3 m, que es el di?metro normal
  .filter$filter <- ifelse(.filter$sec == 1.3 & .filter$arc.circ == 1, 1, 0)
  .filter2 <- subset(.filter, .filter$filter == 1)

  # Aqu? estimo los puntos por cl?ster, tanto sin como aplicando el proceso de point crooping,
  # que le corresponden a 1 m de radio
  .filter2$points.radio <- .filter2$num.points / .filter2$radio
  .filter2$points.radio.hom <- .filter2$num.points.hom / .filter2$radio

  # Average points after point crooping by m of radius
  # Ahora hago el promedio
  .tree$points.m <- mean(.filter2$points.radio)
  .tree$points.m.hom <- mean(.filter2$points.radio.hom)

  # Finalmente calculo los puntos estimados que le corresponder?an a cada ?rbol
  # en funci?n del radio
  .tree$num.points.est <- .tree$points.m * .tree$radius
  .tree$num.points.hom.est <- .tree$points.m.hom * .tree$radius

  # En caso de que no haya un campo que contenda la identificaci?n de la parcela (id)
  if(is.null(data$id)){

    .tree <- .tree[, c("tree", "center.x", "center.y", "center.phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "num.points", "num.points.hom", "num.points.est", "num.points.hom.est", "partial.occlusion")]
    colnames(.tree) <- c("tree", "x", "y", "phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "num.points", "num.points.hom", "num.points.est", "num.points.hom.est", "partial.occlusion")

  } else{

    # En caso de que exista un campo que contenga la identificaci?n de la parcela (id)

    .tree$id <- data$id[1]
    .tree$file <- data$file[1]

    .tree <- .tree[, c("id", "file", "tree", "center.x", "center.y", "center.phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "num.points", "num.points.hom", "num.points.est", "num.points.hom.est", "partial.occlusion")]
    colnames(.tree) <- c("id", "file", "tree", "x", "y", "phi", "phi.left", "phi.right", "horizontal.distance", "dbh", "num.points", "num.points.hom", "num.points.est", "num.points.hom.est", "partial.occlusion")

  }

  # .tree$id <- as.integer(.tree$id)

  # Por último se agrega la tabla de atributos
  if(!is.null(plot.attributes))
    .tree <- merge(.tree, plot.attributes, by = "id", all = FALSE)


  if(isTRUE(save.result)){

    utils::write.csv(.tree,
                     file = file.path(dir.result, "tree.list.tls.csv"),
                     row.names = FALSE)
  }


  #####
  return(.tree)

}








