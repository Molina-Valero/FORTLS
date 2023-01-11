
distance.sampling <- function(tree.tls,
                              id.plots = NULL,
                              strata.attributes = NULL){

  # Funciona con el data frame que se obtiene de la funci?n tree_detection()
  # El argumento plot.radius se refiere al radio de la parcela considerado
  .data <- tree.tls

  # Convert dbh (cm) to International System of Units (m)
  .data$dbh <- .data$dbh / 100


  # Selecting plots especifiesd in the argument id.plots (NULL by default)
  if(is.null(id.plots)) {

    .data <- .data

  } else {

    .data <- .data[is.element(.data$id, id.plots), , drop = FALSE]

  }


  if(is.null(.data$id)){

    # Si las el data de entrada no tiene una columna con la identificai?n de las parcelas,
    # y por tanto contendr?a una solo parcela, se le asigna el valor de 1 por defecto
    .data$id <- 1

  } else {

    # En caso contrario se mantiene la codificaci?n para distinguir entre las
    # distintas parcelas
    .data$id <- .data$id

  }

  # Por último se agrega la tabla de atributos
  if(!is.null(strata.attributes))
    .data <- merge(.data, strata.attributes, by = "stratum")

  # Stratum
  if(is.null(.data$stratum)){

    .data$stratum <- 1

  } else {

    .data$stratum <- .data$stratum

  }


  # According to ddf() funtion form "mrds" package:

  # Details

  # The fitting code has certain expectations about data.
  # It should be a dataframe with at least the following fields named and
  # defined as follows:

  # Study.Area:    name of the study area (by default "forest")
  # Region.Label: strata (if there are several stratum)
  # Sample.Label: point transect identifier (plot id)
  # Effort:       survey effort (1 for all points because they are visited a signle time)
  # object:       unique identifier for each detected tree
  # distance:     radial distance of detection from observer (TLS)
  # OBS:          initials of observer; TLS in this case
  # DBH:          diameter estimated at breast height (1.3 m)


  # Because of that, I am defining the next fields:

  # "Effort"
  .data$Effort <- 1

  # "OBS"
  .data$OBS <- "TLS"

  # "object"
  .data$object <- c(1:nrow(.data))

  # Naming the rest of fields according to ddf() function
  .data <- .data[, c("stratum", "id", "tree", "Effort", "h.dist", "OBS", "object", "dbh"), drop = FALSE]
  colnames(.data) <- c("Region.Label", "Sample.Label", "tree", "Effort", "distance", "OBS", "object", "dbh")


  # Deleting possible missing distances
  .data <- .data[!is.na(.data$distance), , drop = FALSE]


  # Defining empty data frames where results will be added

  # Trees detection probability

  .tree <- data.frame(stratum = as.numeric(), id = as.numeric(), tree = as.numeric(),
                      P.hn = as.numeric(), P.hn.cov = as.numeric(), P.hn.cov = as.numeric(), P.hr.cov = as.numeric())


  # Parameters of detection probability functions

  .par <- data.frame(P.hn.scale = as.numeric(),
                     P.hn.cov.scale.intercept = as.numeric(), P.hn.cov.dbh = as.numeric(),
                     P.hr.scale = as.numeric(), P.hr.shape = as.numeric(),
                     P.hr.cov.scale.intercept = as.numeric(), P.hr.cov.dbh = as.numeric(), P.hr.cov.shape = as.numeric())

  # AIC criterion of detection probability functions fitted
  .treeAIC <- data.frame(P.hn = as.numeric(), P.hn.cov = as.numeric(), P.hr = as.numeric(), P.hr.cov = as.numeric())



  ## Fitting detection probability functions by stratum

  for (i in unique(.data$Region.Label)) {

    # Defining plots location
    # graphics::par(mfrow = c(2, 2))


    # Selecting the correstponding stratum
    .dat <- .data[which(.data$Region.Label == i), , drop = FALSE]

    # Defining maximum distance for fitting detection probability functions
    if(is.null(strata.attributes$plot.radius)){

      .max.dist <- max(.dat$distance)

    } else {

      .max.dist <- strata.attributes$plot.radius[i]

    }


    .dat <- .dat[which(.dat$distance <= .max.dist), , drop = FALSE]


    # Ahora se utiliza la funci?n ddf() del paquete mrds,
    # el cual ha sido desarrollado para este tipo de muestreos

    # Half-Normal fucntion

    .error_hn <- try(suppressMessages(Distance::ds(.dat, truncation = list(left = 1, right = .max.dist), transect = "point", key = "hn", adjustment = NULL)))

    if(class(.error_hn)[1] == "try-error"){

      .hn <- data.frame(fitted = NA, criterion = NA)
      .hn <- list(ddf = .hn)

    } else {

      .hn <- suppressMessages(Distance::ds(.dat, truncation = list(left = 1, right = .max.dist), transect = "point", key = "hn", adjustment = NULL))

      # plot(.hn, main = paste("Half-normal - stratum", i, sep = " "))

    }

    # Half-Normal fucntion with dbh as coovariate

    .error_hn_cov <- try(suppressMessages(Distance::ds(.dat, truncation = list(left = 1, right = .max.dist), transect = "point", key = "hn", formula=~dbh, adjustment = NULL)))

    if(class(.error_hn_cov)[1] == "try-error"){

      .hn_cov <- data.frame(fitted = NA, criterion = NA)
      .hn_cov <- list(ddf = .hn_cov)

    } else {

      .hn_cov <- suppressMessages(Distance::ds(.dat, truncation = list(left = 1, right = .max.dist), transect = "point", key = "hn", formula=~dbh, adjustment = NULL))

      # plot(.hn_cov, main = paste("Half-normal with dbh as coovariate - Stratum", i, sep = " "))

    }

    # Hazard rate

    .error_hr <- try(suppressMessages(Distance::ds(.dat, truncation = list(left = 1, right = .max.dist), transect = "point", key = "hn", adjustment = NULL)))

    if(class(.error_hr)[1] == "try-error"){

      .hr <- data.frame(fitted = NA, criterion = NA)
      .hr <- list(ddf = .hr)

    } else {

      .hr <- suppressMessages(Distance::ds(.dat, truncation = list(left = 1, right = .max.dist), transect = "point", key = "hr", adjustment = NULL))

      # plot(.hr, main = paste("Hazard rate - stratum", i, sep = " "))

    }

    # Hazard rate with dbh as coovariate

    .error_hr_cov <- try(suppressMessages(Distance::ds(.dat, truncation = list(left = 1, right = .max.dist), transect = "point", key = "hr", formula=~dbh, adjustment = NULL)))

    if(class(.error_hr_cov)[1] == "try-error"){

      .hr_cov <- data.frame(fitted = NA, criterion = NA)
      .hr_cov <- list(ddf = .hr_cov)

    } else {

      .hr_cov <- suppressMessages(Distance::ds(.dat, truncation = list(left = 1, right = .max.dist), transect = "point", key = "hr", formula=~dbh, adjustment = NULL))

      # plot(.hr_cov,  main = paste("Hazard rate with dbh as coovariate - stratum", i, sep = " "))

    }

    # Final merge
    # Se obtiene un data frame con los valores de la probabilidad de detecci?n de cada arbol
    .salida <- data.frame(# object = .dat$object,
      P.hn = .hn$ddf$fitted, P.hn.cov = .hn_cov$ddf$fitted,
      P.hr = .hr$ddf$fitted, P.hr.cov = .hr_cov$ddf$fitted)
    .salida <- merge(data.frame(object = .dat$object), .salida, by = "row.names", all = TRUE)
    .salida$P.hn <- ifelse(is.na(.salida$P.hn),  mean(.salida$P.hn, na.rm = TRUE), .salida$P.hn)
    .salida$P.hn.cov <- ifelse(is.na(.salida$P.hn.cov),  mean(.salida$P.hn.cov, na.rm = TRUE), .salida$P.hn.cov)
    .salida$P.hr <- ifelse(is.na(.salida$P.hr),  mean(.salida$P.hr, na.rm = TRUE), .salida$P.hr)
    .salida$P.hr.cov <- ifelse(is.na(.salida$P.hr.cov),  mean(.salida$P.hr.cov, na.rm = TRUE), .salida$P.hr.cov)

    # Se asocian los valores obtenidos a sus respectivos ?rboles
    .tree_n <- merge(.dat, .salida, by = "object")

    .tree_n <- .tree_n[, c("Region.Label", "Sample.Label", "tree", "P.hn", "P.hn.cov", "P.hr", "P.hr.cov"), drop = FALSE]
    colnames(.tree_n) <- c("stratum", "id", "tree", "P.hn", "P.hn.cov", "P.hr", "P.hr.cov")

    .tree <- rbind(.tree, .tree_n)

    # Obtaining values of detection probability functions fitted
    .par_n <- data.frame(P.hn.scale = exp(.hn$ddf$par[1]),
                         P.hn.cov.scale.intercept = exp(.hn_cov$ddf$par[1]), P.hn.cov.dbh = .hn_cov$ddf$par[2],
                         P.hr.scale = exp(.hr$ddf$par[2]), P.hr.shape = exp(.hr$ddf$par[1]),
                         P.hr.cov.scale.intercept = exp(.hr_cov$ddf$par[2]), P.hr.cov.dbh = .hr_cov$ddf$par[3], P.hr.cov.shape = exp(.hr_cov$ddf$par[1]))

    .par <- rbind(.par, .par_n)

    # Data frame con el valor del AIC de los modelos ajustados
    .treeAIC_n <- data.frame(P.hn = .hn$ddf$criterion, P.hn.cov = .hn_cov$ddf$criterion, P.hr = .hr$ddf$criterion, P.hr.cov = .hr_cov$ddf$criterion)

    .treeAIC <- rbind(.treeAIC, .treeAIC_n)

  }

  # Final list
  .tree <- list(tree = .tree, par = .par, AIC = .treeAIC)

  ###
  return(.tree)

}
