
# Add quotes and paste a vector

.quot.past <- function(x, quot = "'", sep = "", collapse = ", ") {
  
  paste(quot, trimws(x), quot, sep = sep, collapse = collapse)
  
}


# Check object class and dimension, and values validity

.check.class <- function(x, x.class, name, col.of = NULL, n = NULL,
                         val = NULL, order.val = TRUE) {
  
  # Check object class
  if (x.class == "integer")
    .check <- is.numeric(x) && all(round(x[!is.na(x)], 0) == x[!is.na(x)])
  else
    .check <- get(paste("is", x.class, sep = "."))(x)
  if (!.check)
    stop(ifelse(is.null(col.of), .quot.past(name),
                paste("Column", .quot.past(name), "of", .quot.past(col.of))),
         " must be a ", x.class,
         ifelse(x.class %in% c("integer", "numeric", "character", "logical") &
                  !is.null(n) & n > 1, " vector", ""),
         ".")
  
  # Check object dimension
  if (!is.null(n)) {
    
    if (!x.class %in% c("data.frame")) {
      
      if (length(x) == 0)
        stop(ifelse(is.null(col.of), .quot.past(name),
                    paste("Column", .quot.past(name), "of",
                          .quot.past(col.of))),
             " must have length different from zero.")
      if (n == 1) {
        
        x <- unique(x)
        if (length(x) > 1) {
          
          warning("Only first value in ",
                  ifelse(is.null(col.of), .quot.past(name),
                         paste("column", .quot.past(name), "of",
                               .quot.past(col.of))),
                  " was taken into account during the execution.",
                  immediate. = TRUE)
          x <- x[1]
          
        }
        
      }
      
    } else {
      
      if (nrow(x) == 0) stop(.quot.past(name), " must have at least one row.")
      if (n == 1) {
        
        x <- unique(x)
        if (nrow(x) > 1) {
          
          warning("Only first row in ", .quot.past(name), " was taken into ",
                  "account during the execution.", immediate. = TRUE)
          x <- x[1, , drop = FALSE]
          
        }
        
      }
      
    }
    
  }
  
  # Check values and order
  if (!is.null(val)) {
    
    .inval <- x[!x %in% val]
    if (length(.inval) > 0)
      stop("Value(s) ", .quot.past(.inval), " in ",
           ifelse(is.null(col.of), .quot.past(name),
                  paste("column", .quot.past(name), "of",
                        .quot.past(col.of))),
           " is(are) not valid.\nValid value(s) is(are) ", .quot.past(val), ".")
    if (order.val)
      x <- val[val %in% x]
    
    
  }
  
  return(x)
  
}


# Check mandatory columns and their classes for a data.frame

.check.col <- function(x, x.mand, x.class, name, def = NULL,
                       na.action = c("stop", "warning")[1]) {
  
  # Check mandatory columns existence, assign values by default (if necessary),
  # and remove non necessary columns
  .miss <- colnames(x.mand)[apply(x.mand, 2, any)]
  .miss <- .miss[!.miss %in% colnames(x)]
  if (length(.miss) > 0)
    stop("According to specified arguments, ", .quot.past(name), " must have ",
         "column(s) named: ", .quot.past(.miss), ".")
  if (!is.null(def)) {
    
    .miss <- colnames(def)[!colnames(def) %in% colnames(x)]
    if (length(.miss) > 0)
      x <- cbind(x, do.call(cbind, lapply(def[, .miss], rep, times = nrow(x))))
    
  }
  x <- x[, colnames(x.mand), drop = FALSE]
  
  # Check mandatory columns classes
  for (.i in colnames(x)[colnames(x) %in% names(x.class)])
    x[, .i] <- .check.class(x = x[, .i], x.class = x.class[.i], name = .i,
                            col.of = name, n = nrow(x))
  
  # Check all values are no-NA
  .col.na <- colnames(x)[apply(is.na(x), 2, any)]
  if (length(.col.na) > 0) {
    
    if (na.action %in% "stop")
      stop("Column(s) ", .quot.past(.col.na), " of ", .quot.past(name),
           " has(have) missing values.")
    else if (na.action %in% "warning")
      warning("Column(s) ", .quot.past(.col.na),  "of ", .quot.past(name),
              " has(have) missing values.", immediate. = TRUE)
    
  }
  
  return(x)
  
}


# Compute number of decimals places

.decimals <- function(x){
  
  x <- format(x, scient = FALSE)
  ifelse(!base::grepl(".", x, fix = TRUE), 0,
         nchar(base::strsplit(x, ".", fix = TRUE)[[1]][2]))
  
}


# Compute radius, k and BAF, and tree variables according to plot design(s) and
# 'tree.var'. Currently available tree variables: basal area (g) and volume (v)

.tree.calc <- function(tree, plot.design, tree.var,
                       v.calc = c("coeff", "parab")[2]) {
  
  # Create data.frame where results will be saved
  .col.names <- plot.design
  if ("k.tree" %in% names(plot.design))
    .col.names <- c(.col.names,
                    paste("radius", plot.design["k.tree"], sep = "."))
  .col.names <- c(.col.names, tree.var)
  tree <- cbind(tree, matrix(nrow = nrow(tree), ncol = length(.col.names),
                             dimnames = list(NULL, .col.names)))
  
  # Compute radius (m) for fixed area plot design
  if (plot.design["fixed.area"] %in% colnames(tree))
    tree[, plot.design["fixed.area"]] <- tree[, "h.dist"]
  
  # Compute k (trees) and associated radius (m) for k-tree design
  if (plot.design["k.tree"] %in% colnames(tree)) {
    
    # Order trees by horizontal distance
    .ord <- order(tree[, "h.dist"], decreasing = F)
    # Save order as k
    tree[.ord, plot.design["k.tree"]] <- 1:nrow(tree)
    # Compute associated radius
    .dist <- c(tree[.ord[-1], "h.dist"],
               utils::tail(tree[.ord, "h.dist"], n = 1))
    tree[.ord, paste("radius", plot.design["k.tree"], sep = ".")] <-
      (tree[.ord, "h.dist"] + .dist) / 2
    
  }
  
  # Compute BAF (m2/ha) threshold for angle-count design
  if (plot.design["angle.count"] %in% colnames(tree))
    tree[, plot.design["angle.count"]] <- 2500 / (tree[, "h.dist"] /
                                                    tree[, "dbh"]) ^ 2
  
  # Compute basal area (m2)
  if ("g" %in% colnames(tree)) tree[, "g"] <- (pi / 4) * tree[, "dbh"] ^ 2
  
  # Compute volume (m3)
  if ("v" %in% colnames(tree)) {
    
    if (v.calc == "coeff") {
      
      # Coefficient of 0.45
      tree[, "v"] <- (pi / 4) * tree[, "dbh"] ^ 2 * tree[, "h"] * 0.45
      
    } else if (v.calc == "parab") {
      
      # Paraboloid
      tree[, "v"] <- pi * (tree[, "h"] ^ 2 / 2) *
        ((tree[, "dbh"] / 2) ^ 2 / (tree[, "h"] - 1.3) ^ 2)
      
    } else
      stop("Argument for tree volume calculation must be 'coeff' or 'parab'.")
    
  }
  
  return(tree)
  
}


# Custom rounding of numbers

.customCeiling <- function(x, Decimals = 1) {
  
  ceiling(x * 10 ^ Decimals) / 10 ^ Decimals
  
}

.customFloor <- function(x, Decimals = 1) {
  
  floor(x * 10 ^ Decimals) / 10 ^ Decimals
  
}


# Custom format for numbers

.format.numb <- function(x, dec) {
  
  format(x, trim = TRUE, nsmall = dec)
  
}


# Compute several weighted mean functions for a numeric vector

.wmean.calculation <- function(data, w, mean.names) {
  
  sapply(mean.names,
         function(x, data, w) {
           get(paste("weighted_mean", x, sep = "_"))(data, w)
         },
         data = data, w = w)
  
}


# Compute expansion factors with occlusion corrections and estimate stand
# variables per ha, and compute mean diameters and heights

.stand.calc <- function(val, plot.design, tree, var.metr, ds.meth,
                        var.field.user, mean.d, mean.h) {
  
  # Select data according to specified radius/k/BAF value, and define radius
  if (names(plot.design) %in% c("fixed.area", "k.tree")) {
    
    tree <- tree[tree[, plot.design] <= val, , drop = FALSE]
    .radius <- switch(names(plot.design), fixed.area = val,
                      k.tree = tree[nrow(tree), paste("radius", plot.design,
                                                      sep = ".")])
    
  } else if (names(plot.design) %in% "angle.count")
    tree <- tree[tree[, plot.design] >= val, , drop = FALSE]
  
  # Add auxiliary column for density calculations
  .col.names <- c(paste("N", c("tls", names(ds.meth), "sh", "pam"), sep = "."),
                  "N")
  .col.names <- .col.names[.col.names %in% var.metr]
  if (length(.col.names) > 0) tree <- cbind(tree, n = 1)
  
  # Initialize data.frame where results will be saved: id, radius/k/BAF, and
  # variables and/or metrics according to 'var.metr'
  .col.names <- c("id", plot.design, var.metr)
  if ("stratum" %in% colnames(tree)) .col.names <- c("stratum", .col.names)
  stand <- data.frame(matrix(nrow = length(val), ncol = length(.col.names),
                             dimnames = list(names(val), .col.names)),
                      stringsAsFactors = FALSE)
  
  # Plot ID and radius/k/BAF
  stand[, "id"] <- unique(tree[, "id"])
  stand[, plot.design] <- val
  
  
  # Compute expansion factors to convert units per plot into units per ha, ----
  # with occlusion corrections, and estimate stand variables per ha: density
  # (trees/ha), basal area (m2/ha), volume (m3/ha), and optionally volume
  # (m3/ha) and biomass (Mg/ha) derived from trees attributes provided by the
  # user
  
  # Identify stand variables to be computed, and associated tree variable and
  # expansion factor
  .col.names <- c(sapply(c("N", "G", "V"), paste,
                         c("tls", names(ds.meth), "sh", "pam"), sep = "."),
                  "N", "G", "V", names(var.field.user))
  .col.names <- matrix("", nrow = length(.col.names), ncol = 2,
                       dimnames = list(.col.names, c("var", "ef")))
  for (.i in c("N", "G", "V"))
    .col.names[c(sapply(.i, paste, c("tls", names(ds.meth), "sh", "pam"),
                        sep = "."), .i), "var"] <-
    switch(.i, N = "n", G = "g", V = "v")
  for (.i in names(var.field.user))
    .col.names[.i, "var"] <- var.field.user[[.i]]
  .col.names[c(paste(c("N", "G", "V"), "tls", sep = "."), "N", "G", "V",
               names(var.field.user)), "ef"] <- "EF"
  for (.i in c(names(ds.meth), "sh", "pam"))
    .col.names[paste(c("N", "G", "V"), .i, sep = "."), "ef"] <- paste("EF", .i,
                                                                      sep =".")
  .col.names <- .col.names[rownames(.col.names) %in% colnames(stand), ,
                           drop = FALSE]
  if (nrow(.col.names) > 0) {
    
    # Compute expansion factor according to plot design (all plot designs)
    if (names(plot.design) %in% c("fixed.area", "k.tree"))
      .EF <- 10000 / (pi * .radius ^ 2)
    else if (names(plot.design) %in% "angle.count")
      .EF <- val / (pi * (tree[, "dbh"] / 2) ^ 2)
    .EF <- matrix(.EF, nrow = nrow(tree), dimnames = list(NULL, "EF"))
    
    # Compute expansion factor using distance sampling based correction (fixed
    # area and k-tree plot designs)
    .col.names.2 <- ds.meth[paste("EF", names(ds.meth), sep = ".") %in%
                              .col.names[, "ef"]]
    if (length(.col.names.2) > 0) {
      
      .P <- tree[nrow(tree), .col.names.2, drop = FALSE]
      colnames(.P) <- paste("EF", names(.col.names.2), sep = ".")
      if (nrow(tree) > 1) .P <- apply(.P, 2, rep, times = nrow(tree))
      .EF <- cbind(.EF, .EF[, "EF"] / .P)
      
    }
    
    # Compute expansion factor with correction of the shadowing effect (fixed
    # area and k-tree plot designs)
    if ("EF.sh" %in% .col.names[, "ef"]) {
      
      # Compute shadow area for each tree
      .sh <- ifelse(tree[, "partial.occlusion"] == 0,
                    
                    # Non-occluded trees
                    ((((pi * .radius ^ 2) - (pi * tree[, "h.dist"] ^ 2)) /
                        (2 * pi)) * atan2(tree[, "dbh"], tree[, "h.dist"])) -
                      ((pi * (tree[, "dbh"] / 2) ^ 2) / 2),
                    
                    # Occluded trees
                    ((((pi * .radius ^ 2) - (pi * tree[, "h.dist"] ^ 2)) /
                        (2 * pi)) * (tree[, "wide"])) -
                      (((pi * (tree[, "dbh"] / 2) ^ 2) * tree[, "wide"]) /
                         atan2(tree[, "dbh"], tree[, "h.dist"])))
      
      # Correction: sh = 0 if tree is not completely inside the plot
      .sh[tree[, "h.dist"] + tree[, "dbh"] / 2 > .radius] <- 0
      
      # Compute expansion factor
      .sh <- 1 + (sum(.sh, na.rm = TRUE) / (pi * .radius ^ 2))
      .EF <- cbind(.EF, EF.sh = .EF[, "EF"] * .sh)
      
    }
    
    # Compute expansion factor with gap probability attenuation correction
    # (angle-count plot design)
    if ("EF.pam" %in% .col.names[, "ef"]) {
      
      .De <- mean(tree[, "dbh"]) *
        (1 + (stats::sd(tree[, "dbh"], na.rm = TRUE) /
                mean(tree[, "dbh"], na.rm = TRUE)) ^ 2) ^ 0.5
      .t <- (((val / (pi * (tree[, "dbh"] / 2) ^ 2)) / 10000) * .De *
               tree[, "dbh"]) / (2 * sqrt(val / 10000))
      .Ft <- (2 / .t ^ 2) * (1 - exp(-.t) * (1 + .t))
      .EF <- cbind(.EF, EF.pam = .EF[, "EF"] / .Ft)
      
    }
    
    # Estimate variables per ha
    stand[, rownames(.col.names)] <-
      apply(tree[, .col.names[, "var"], drop = FALSE] *
              .EF[, .col.names[, "ef"], drop = FALSE], 2, sum)
    
  }
  
  
  # Compute mean diameters and heights ----
  
  # Compute mean diameters
  .col.names <- sapply(names(mean.d), paste, c("", ".tls"), sep = "")
  .col.names <- .col.names[apply(apply(.col.names, 1:2,  "%in%",
                                       colnames(stand)), 1, any), ]
  if (length(.col.names) > 0)
    stand[, .col.names] <- .wmean.calculation(data = tree[, "dbh"],
                                              w = rep(1, nrow(tree)),
                                              mean.names = mean.d)
  
  # Compute mean heights
  .col.names <- sapply(names(mean.h), paste, c("", ".tls"), sep = "")
  .col.names <- .col.names[apply(apply(.col.names, 1:2,  "%in%",
                                       colnames(stand)), 1, any), ]
  if (length(.col.names) > 0)
    stand[, .col.names] <- .wmean.calculation(data = tree[, "h"],
                                              w = rep(1, nrow(tree)),
                                              mean.names = mean.h)
  
  
  # Compute number of points
  .col.names <- paste("num.points", c("", ".est", ".hom", ".hom.est"), sep = "")
  .col.names <- .col.names[.col.names %in% colnames(stand)]
  if (length(.col.names) > 0)
    stand[, .col.names] <- apply(tree[, .col.names, drop = FALSE], 2, sum)
  
  return(stand)
  
}


# LiDAR metrics: 'mean', 'max', 'min', 'sd', 'var', 'mode', 'kurtosis',
# 'skewness', 'perc_on_mean', 'perc_on_mode', 'weibull_b', 'weibull_c'

.getmode <- function(v) {
  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
  
}

.c_function <- function(c, media, varianza){
  
  varianza - (media^2) * (gamma(1 + 2 / c) -
                            (gamma(1 + 1 / c))^2) / (gamma(1 + 1 / c))^2
  
}

.points.metrics <- function(rho_seq, data, metr) {
  
  # Restrict data
  data <- data[data[, "z"] > 0.1, , drop = FALSE]
  
  # Compute metrics
  .metr <- lapply(rho_seq,
                  function(rho, data, metr) {
                    
                    .sub <- data[data[, "rho"] <= rho, "z"]
                    
                    .metr <- rep(NA, length(metr))
                    names(.metr) <- metr
                    
                    if ("mean" %in% names(.metr)) .metr["mean"] <- mean(.sub)
                    if ("max" %in% names(.metr)) .metr["max"] <- max(.sub)
                    if ("min" %in% names(.metr)) .metr["min"] <- min(.sub)
                    if ("sd" %in% names(.metr)) .metr["sd"] <- sd(.sub)
                    if ("var" %in% names(.metr)) .metr["var"] <- var(.sub)
                    if ("mode" %in% names(.metr))
                      .metr["mode"] <- .getmode(.sub)
                    if ("kurtosis" %in% names(.metr))
                      .metr["kurtosis"] <- moments::kurtosis(.sub)
                    if ("skewness" %in% names(.metr))
                      .metr["skewness"] <- moments::skewness(.sub)
                    
                    if ("perc_on_mode" %in% names(.metr))
                      .metr["perc_on_mode"] <- mean(.sub > .metr["mode"]) * 100
                    if ("perc_on_mean" %in% names(.metr))
                      .metr["perc_on_mean"] <- mean(.sub > .metr["mean"]) * 100
                    
                    if ("weibull_c" %in% names(.metr)) {
                      
                      .metr["weibull_c"] <-
                        stats::uniroot(.c_function, media = .metr["mean"],
                                       varianza = .metr["var"],
                                       interval = c(.metr["min"],
                                                    .metr["max"]))$root
                      
                    }
                    if ("weibull_b" %in% names(.metr)) {
                      
                      .metr["weibull_b"] <-
                        .metr["mean"] / gamma(1 + 1 / .metr["weibull_c"])
                      
                    }
                    
                    return(.metr)
                    
                  },
                  data = data, metr = metr)
  
  .metr <- do.call(rbind, .metr)
  
  return(.metr)
  
}


# Simulate plots and compute stand variables and metrics

.sim.calc <- function(funct, tree.list.tls, tree.list.ds, tree.list.field,
                      plot.design, plot.parameters, scan.approach, var.metr,
                      v.calc, dbh.min, h.min, max.dist, dir.data, save.result,
                      dir.result) {
  
  
  # Define available values for arguments and other auxiliary objects ----
  
  # Define a list containing mandatory trees' database(s) according to 'function
  # case'. Currently available 'function cases': 'sim', 'metr' and 'est'
  .funct <- list(sim = c(field = "field", tls = "TLS"), metr = c(tls = "TLS"),
                 est = c(tls = "TLS"))
  
  # Define a character vector containing index name (radius, k or BAF) for each
  # available plot design. Currently available plot designs: 'fixed.area',
  # 'k.tree' and 'angle.count'
  .plot.design <- c(fixed.area = "radius", k.tree = "k", angle.count = "BAF")
  
  # Define a list containing mandatory columns in 'tree.list.tls' argument
  # associated to each available scan approach. Currently available scan
  # approaches: 'single' and 'multi'
  .scan.approach <- list(single = c("phi.left", "phi.right"),
                         multi = character())
  
  # Define a character vector containing mandatory columns in 'tree.list.ds'
  # argument associated to available distance sampling methodologies. Currently
  # available methodologies: 'hn' (half normal function), 'hn.cov' (half normal
  # function with dbh as covariate), 'hr' (half rate function) and 'hr.cov'
  # (half rate function with dbh as covariate)
  .ds.meth <- c(hn = "P.hn", hr = "P.hr", hn.cov = "P.hn.cov",
                hr.cov = "P.hr.cov")
  
  # Define character vectors containing the mean function to be used for each
  # mean diameter/height computation. Currently available means: arithmetic,
  # quadratic, geometric and harmonic
  .mean.d <- c(d = "arit", dg = "sqrt", dgeom = "geom", dharm = "harm")
  .mean.h <- c(h = "arit", hg = "sqrt", hgeom = "geom", hharm = "harm")
  
  # Define a numeric vector containing probabilities for height percentiles
  # calculation
  .prob = c(1, 5, 10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90, 95, 99)
  
  # Define a list containing mandatory columns in 'tree.list.field' argument
  # associated to each available field variables to be computed from optional
  # trees attributes provided by the user. Currently available field variables:
  # 'V.user' and 'W.user'
  .var.field.user <- list(V.user = "v.user", W.user = "w.user")
  
  # Define a list containing the available TLS variables and metrics according
  # to 'function case', 'tree.list.ds' availability, scan approach and plot
  # design(s), and the available field variables according to 'function case'
  # and trees variables provided by user in 'tree.list.field'.
  # Currently available:
  # - TLS variables and metrics (1): 'N.tls', 'N.hn' (2), 'N.hr' (2),
  #   'N.hn.cov' (2), 'N.hr.cov' (2), 'N.sh' (3), 'N.pam' (4); 'G.tls',
  #   'G.hn' (2), 'G.hr' (2), 'G.hn.cov' (2), 'G.hr.cov' (2), 'G.sh' (3),
  #   'G.pam' (4); 'V.tls', 'V.hn' (2), 'V.hr' (2), 'V.hn.cov' (2),
  #   'V.hr.cov'(2), 'V.sh' (3), 'V.pam' (4); 'd.tls', 'dg.tls', 'dgeom.tls',
  #   'dharm.tls'; 'h.tls', 'hg.tls', 'hgeom.tls', 'hharm.tls'; 'd.0.tls',
  #   'dg.0.tls', 'dgeom.0.tls', 'dharm.0.tls'; 'h.0.tls', 'hg.0.tls',
  #   'hgeom.0.tls', 'hharm.0.tls'; 'num.points', 'num.points.est',
  #   'num.points.hom', 'num.points.hom.est'; 'P01', 'P05', 'P10', 'P20', 'P25',
  #   'P30', 'P40', 'P50', 'P60', 'P70', 'P75', 'P80', 'P90', 'P95', 'P99';
  #   'mean', 'max', 'min', 'sd', 'var', 'mode', 'kurtosis', 'skewness',
  #   'perc_on_mode', 'perc_on_mean', 'weibull_c', 'weibull_b'.
  # - Field variables (5): 'N'; 'G'; 'V', 'V.user' (6); 'W.user' (6); 'd', 'dg',
  #   'dgeom', 'dharm'; 'h', 'hg', 'hgeom', 'hharm'; 'd.0', 'dg.0', 'dgeom.0',
  #   'dharm.0'; 'h.0', 'hg.0', 'hgeom.0', 'hharm.0'.
  # Remarks. (1) For function case 'est' only 'N.tls' and 'G.tls' are computed.
  # (2) Only for circular fixed area and k-tree plot designs if scan approach is
  # 'simple' and 'tree.list.ds' is not NULL. (3) Only for circular fixed area
  # and k-tree plot design if scan approach is 'simple'. (4) Only for
  # angle-count plot design if scan approach is 'simple'. (5) Only for function
  # case 'sim'. (6) Only if corresponding trees variable ('v.user' for 'V.user',
  # 'w.user' for 'W.user') are provided by user in 'tree.list.field' argument.
  .var.metr <- list(tls = NULL, field = NULL)
  .var.metr$tls <- c(sapply(c("N", "G", "V"), paste,
                            c("tls", names(.ds.meth), "sh", "pam"), sep = "."),
                     t(sapply(names(c(.mean.d, .mean.h)),
                              paste, c(".tls", ".0.tls"), sep = "")),
                     paste("num.points", c("", ".est", ".hom", ".hom.est"),
                           sep = ""),
                     sprintf("P%02i", .prob), "mean", "max", "min", "sd", "var",
                     "mode", "kurtosis", "skewness", "perc_on_mode",
                     "perc_on_mean", "weibull_c", "weibull_b")
  .var.metr$field <- c("N", "G", "V", names(.var.field.user),
                       t(sapply(names(c(.mean.d, .mean.h)), paste, c("", ".0"),
                                sep = "")))
  
  # Define available calculations for trees volume. Currently available
  # calculations: 'coeff' and 'parab'
  .v.calc <- c("coeff", "parab")
  
  # Define values by default for certain columns in 'plot.parameters'
  .plot.parameters <- data.frame(radius.incr = 0.1, k.incr = 1, BAF.incr = 0.1,
                                 num.trees = 100)
  
  
  # Check arguments 'funct', 'plot.design', 'scan.approach', 'v.calc', ----
  # 'dbh.min', 'h.min', 'max.dist', 'save.result', 'dir.data', 'dir.result'
  
  # 'funct' must be a character string containing the 'function  case'
  funct <- .check.class(x = funct, x.class = "character", name = "funct",
                        n = 1, val = names(.funct))
  
  # 'plot.design' must be a character vector containing the plot designs to be
  # computed
  plot.design <- .check.class(x = plot.design, x.class = "character",
                              name = "plot.design", n = Inf,
                              val = names(.plot.design))
  
  # 'scan.approach' must be a character string indicating the scan approach used
  # for obtaining LAS files
  scan.approach <- .check.class(x = scan.approach, x.class = "character",
                                name = "scan.approach", n = 1,
                                val = names(.scan.approach))
  
  # 'v.calc' must be a character string indicating how to calculate trees volume
  v.calc <- .check.class(x = v.calc, x.class = "character", name = "v.calc",
                         n = 1, val = .v.calc)
  
  # 'dbh.min', 'h.min' and 'max.dist' must be positive numeric values
  for (.i in c("dbh.min", "h.min", "max.dist")) {
    
    # Check object class and dimension
    assign(.i, .check.class(x = get(.i), x.class = "numeric", name = .i, n = 1))
    
    # Check value is  no-NA and strictly positive
    if (is.na(get(.i)) || get(.i) <= 0)
      stop(.quot.past(.i), " argument must be strictly positive.")
    
  }
  
  # 'save.result' must be a logical
  save.result <- .check.class(x = save.result, x.class = "logical",
                              name = "save.result", n = 1)
  
  # 'dir.data' and 'dir.result' must be NULL or a character string containing
  # the absolute path to a existing directory
  # Remark. If 'dir.data' is NULL, working directory is assigned by default. If
  # 'save.result' is FALSE, 'dir.result' is forced to be NULL; otherwise, if
  # 'dir.result' is NULL, working directory is assigned by default
  for (.i in c("dir.data", "dir.result")) {
    
    if (.i %in% "dir.result" & !save.result) assign(.i, NULL)
    else {
      
      if (is.null(get(.i))) assign(.i, getwd())
      else {
        
        # Check object class and dimension
        assign(.i, .check.class(x = get(.i), x.class = "character", name = .i,
                                n = 1))
        
        
        # Check directory existence
        if (!dir.exists(get(.i))) stop(.quot.past(.i), " directory must exist.")
        
      }
      
    }
    
  }
  
  
  # Check argument 'var.metr' ----
  
  # Restrict available TLS variables and metrics according to 'function case',
  # 'tree.list.ds', scan approach and plot design(s)
  if (!"tls" %in% names(.funct[[funct]])) .var.metr$tls <- NULL
  else {
    
    if (funct %in% "est")
      .var.metr$tls <- .var.metr$tls[.var.metr$tls %in%
                                       paste(c("N", "G"), "tls", sep = ".")]
    if (is.null(tree.list.ds))
      .var.metr$tls <- .var.metr$tls[!.var.metr$tls %in%
                                       sapply(c("N", "G", "V"), paste,
                                              names(.ds.meth), sep = ".")]
    if (!scan.approach %in% "single")
      .var.metr$tls <- .var.metr$tls[!.var.metr$tls %in%
                                       sapply(c("N", "G", "V"), paste,
                                              c(names(.ds.meth), "sh", "pam"),
                                              sep = ".")]
    .var.metr$tls <- matrix(TRUE, nrow = length(plot.design),
                            ncol = length(.var.metr$tls),
                            dimnames = list(plot.design, .var.metr$tls))
    .var.metr$tls[!rownames(.var.metr$tls) %in% c("fixed.area", "k.tree"),
                  colnames(.var.metr$tls) %in% sapply(c("N", "G", "V"), paste,
                                                      c(names(.ds.meth), "sh"),
                                                      sep = ".")] <- FALSE
    .var.metr$tls[!rownames(.var.metr$tls) %in% "angle.count",
                  colnames(.var.metr$tls) %in% paste(c("N", "G", "V"), "pam",
                                                     sep = ".")] <- FALSE
    .var.metr$tls <- .var.metr$tls[, apply(.var.metr$tls, 2, any), drop = FALSE]
    
  }
  
  # Restrict field variables according to 'function case' and trees variables
  # provided by user in 'tree.list.field'
  if (!"field" %in% names(.funct[[funct]])) .var.metr$field <- NULL
  else {
    
    if (!is.data.frame(tree.list.field) ||
        any(!.var.field.user %in% colnames(tree.list.field)))
      .var.metr$field <-
        .var.metr$field[!.var.metr$field %in%
                          names(.var.field.user)[!.var.field.user %in%
                                                   colnames(tree.list.field)]]
    .var.metr$field <- matrix(TRUE, nrow = length(plot.design),
                              ncol = length(.var.metr$field),
                              dimnames = list(plot.design, .var.metr$field))
    
  }
  
  # 'var.metr' must be a list with one or two elements named 'tls' and 'field'
  var.metr <- .check.class(x = var.metr, x.class = class(.var.metr),
                           name = 'var.metr', n = Inf)
  if (length(var.metr) != length(.var.metr))
    stop("'var.metr' must have length equal to ", length(.var.metr), ".")
  .miss <- names(.var.metr)[! names(.var.metr) %in% names(var.metr)]
  if (length(.miss) > 0)
    stop("Element(s) named ", .quot.past(.miss), " is(are) missing in ",
         "'var.metr' argument.")
  
  # 'var.metr$tls' must be NULL or a character vector containing the TLS
  # variables and metrics to be computed
  # Remark. If '.var.metr$tls' is NULL, 'var.metr$tls' is forced to be NULL;
  # otherwise, if 'var.metr$tls', all possible TLS variables and metrics
  # are assigned by default
  for (.i in names(var.metr)) {
    
    if (!.i %in% names(.var.metr)) var.metr[[.i]] <- NULL
    else {
      
      if (is.null(var.metr[[.i]]))
        var.metr[[.i]] <- colnames(.var.metr[[.i]])
      else
        var.metr[[.i]] <- .check.class(x = var.metr[[.i]],
                                       x.class = "character",
                                       name = paste("var.metr", .i, sep ="$"),
                                       n = Inf, val = colnames(.var.metr[[.i]]))
      
    }
    
  }
  
  # Remove plot design(s) with no metric/variable to be computed
  .inval <- do.call(cbind, .var.metr)[, do.call(c, var.metr), drop = FALSE]
  .inval <- rownames(.inval)[apply(!.inval, 1, all)]
  if (length(.inval) > 0) {
    
    warning("Plot design(s) ", .quot.past(.inval), " in 'plot.design' ",
            "argument was(were) discarded, since no associated metric or ",
            "variable needs to be computed.", immediate. = TRUE)
    plot.design <- plot.design[!plot.design %in% .inval]
    
  }
  
  
  # Check 'plot.parameters' argument ----
  
  # Define a logical value indicating if results will be computed by stratum
  # (only for function case 'metr' when 'stratum' column is included in
  # 'plot.parameters' argument)
  .by.stratum <- funct %in% "metr" &
    (is.data.frame(plot.parameters) & "stratum" %in% colnames(plot.parameters))
  
  # Define a logical vector containing all possible mandatory columns in
  # 'plot.parameters' according to 'function case', '.by.stratum' and plot
  # design(s)
  if (funct %in% c("sim", "est")) {
    
    .col.mand <- sapply(.plot.design[plot.design], paste, c("max", "incr"),
                        sep = ".")
    if (funct %in% "sim") .col.mand <- c(.col.mand, "num.trees")
    .col.mand <- matrix(FALSE, nrow = 1, ncol = length(.col.mand),
                        dimnames = list(NULL, .col.mand))
    .col.mand[, paste(.plot.design[plot.design], "max", sep = ".")] <- TRUE
    
  } else {
    
    .col.mand <- c(.plot.design[plot.design], "num.trees")
    .col.mand <- matrix(FALSE, nrow = 1, ncol = length(.col.mand),
                        dimnames = list(NULL, .col.mand))
    .col.mand[, .plot.design[plot.design]] <- TRUE
    
  }
  if (.by.stratum) .col.mand <- cbind("stratum" = TRUE, .col.mand)
  
  # Define a character vector containing (if necessary) the class of all
  # possible mandatory columns in 'plot.parameters'
  .col.class <- rep("numeric", ncol(.col.mand))
  names(.col.class) <- colnames(.col.mand)
  .col.class[names(.col.class) %in% "stratum"] <- NA
  .col.class[names(.col.class) %in%
               c(paste("k", c("max", "incr"), sep = "."), "num.trees")] <-
    "integer"
  .col.class <- .col.class[!is.na(.col.class)]
  
  # 'plot.parameters' must be a data.frame with at least one row
  plot.parameters <- .check.class(x = plot.parameters, x.class = "data.frame",
                                  name = "plot.parameters",
                                  n = ifelse(funct %in% c("sim", "est"), 1,
                                             Inf))
  
  # Check mandatory columns existence, assign values by default (if necessary),
  # remove non necessary columns, check mandatory columns classes and check all
  # values are no-NA
  plot.parameters <- .check.col(x = plot.parameters, x.mand = .col.mand,
                                x.class = .col.class, name = "plot.parameters",
                                def = .plot.parameters)
  
  # Check there are no duplicated strata
  if (.by.stratum) {
    
    # Define stratum code
    .stratum <- plot.parameters[, "stratum", drop = FALSE]
    .stratum <- cbind(code = apply(.stratum, 1, .quot.past), .stratum)
    
    # Check duplicated
    .dupl <- unique(.stratum[duplicated(.stratum[, "code"]), "code"])
    if (length(.dupl) > 0)
      stop("'plot.parameters' argument has several different rows for the ",
           "stratum(strata) ", .quot.past(.dupl, quot = ""), ".")
    
  }
  
  # Check 'radius', 'radius.max', 'k', 'k.max', 'BAF', 'BAF.max' and/or
  # 'num.trees' are strictly positive
  .inval <- c(sapply(.plot.design, paste, c("", ".max"), sep = ""), "num.trees")
  .inval <- colnames(plot.parameters)[colnames(plot.parameters) %in% .inval]
  .inval <- .inval[!apply(plot.parameters[, .inval, drop = FALSE] > 0, 2, all)]
  if (length(.inval) > 0)
    stop("All values in column(s) ", .quot.past(.inval), " of ",
         "'plot.parameters' argument must be strictly positive.")
  
  # Check 'radius.incr', 'k.incr" and 'BAF.incr' are positive
  .inval <- paste(.plot.design, "incr", sep = ".")
  .inval <- colnames(plot.parameters)[colnames(plot.parameters) %in% .inval]
  .inval <- .inval[!apply(plot.parameters[, .inval, drop = FALSE] >= 0, 2, all)]
  if (length(.inval) > 0)
    stop("All values in column(s) ", .quot.past(.inval), " of ",
         "'plot.parameters' argument must be positive.")
  
  # Add number of decimals places to be considered
  for (.i in .plot.design[plot.design]) {
    
    .dec <- sapply(plot.parameters[, ifelse(funct %in% c("sim", "est"),
                                            paste(.i, "incr", sep = "."), .i)],
                   .decimals)
    .dec <- matrix(.dec, ncol = 1,
                   dimnames = list(NULL, paste(.i, "dec", sep = ".")))
    plot.parameters <- cbind(plot.parameters, .dec)
    
  }
  
  
  # Check 'tree.list.tls' argument ----
  
  if (!"tls" %in% names(.funct[[funct]])) tree.list.tls <- NULL
  else {
    
    # Define a logical value indicating if results will be computed by stratum
    # (only for function case 'metr' when 'stratum' column is included in
    # 'plot.parameters' argument) or if strata will be differentiated in charts
    # with different colours (only for function case 'est' when 'stratum' column
    # is included in 'tree.list.tls' argument)
    .by.stratum.2 <- .by.stratum |
      (funct %in% "est" & (is.data.frame(tree.list.tls) &
                             "stratum" %in% colnames(tree.list.tls)))
    
    # Define a logical matrix containing all possible mandatory columns in
    # 'tree.list.tls' according to 'function case', '.by.stratum.2' and
    # 'var.metr'
    .col.mand <- c("id", "file", "tree", unlist(.scan.approach), "h.dist",
                   "dbh", "h",
                   paste("num.points", c("", ".est", ".hom", ".hom.est"),
                         sep = ""),
                   "partial.occlusion")
    if (.by.stratum.2) .col.mand <- c("stratum", .col.mand)
    .col.mand <- matrix(FALSE, nrow = length(var.metr$tls),
                        ncol = length(.col.mand),
                        dimnames = list(var.metr$tls, .col.mand))
    .col.mand[, colnames(.col.mand) %in%
                c("stratum", "id", "tree", "h.dist", "dbh", "h")] <- TRUE
    .col.mand[rownames(.col.mand) %in%
                c(sprintf("P%02i", .prob), "mean", "max", "min", "sd", "var",
                  "mode", "kurtosis", "skewness", "perc_on_mode",
                  "perc_on_mean", "weibull_c", "weibull_b"),
              colnames(.col.mand) %in% "file"] <- TRUE
    .col.mand[rownames(.col.mand) %in% paste(c("N", "G", "V"), "sh", sep = "."),
              colnames(.col.mand) %in%
                c(unlist(.scan.approach), "partial.occlusion")] <- TRUE
    for (.i in paste("num.points", c("", ".est", ".hom", ".hom.est"), sep = ""))
      .col.mand[rownames(.col.mand) %in% .i,
                colnames(.col.mand) %in% .i] <- TRUE
    .col.mand <- .col.mand[, apply(.col.mand, 2, any), drop = FALSE]
    
    # Define a character vector containing (if necessary) the class of all
    # possible mandatory columns in 'tree.list.tls'
    .col.class <- rep("numeric", ncol(.col.mand))
    names(.col.class) <- colnames(.col.mand)
    .col.class[names(.col.class) %in% c("stratum", "id", "file", "tree")] <- NA
    .col.class <- .col.class[!is.na(.col.class)]
    
    # 'tree.list.tls' must be a data.frame with at least one row
    tree.list.tls <- .check.class(x = tree.list.tls, x.class = "data.frame",
                                  name = "tree.list.tls", n = Inf)
    
    # Check mandatory columns existence, remove non necessary columns, check
    # mandatory columns classes and check all values are no-NA
    tree.list.tls <- .check.col(x = tree.list.tls, x.mand = .col.mand,
                                x.class = .col.class, name = "tree.list.tls")
    
    # Define code for plot IDs
    .id <- unique(tree.list.tls[, "id", drop = FALSE])
    .id <- cbind(code = apply(.id, 1, .quot.past), .id)
    
    # Check strata
    if (.by.stratum.2) {
      
      # Check strata are included in 'plot.parameters' (only for function case
      # 'metr' when 'stratum' column is included in 'plot.parameters' argument)
      if (.by.stratum) {
        
        .miss <- apply(unique(tree.list.tls[, "stratum", drop = FALSE]), 1,
                       .quot.past)
        .miss <- .miss[!.miss %in%
                         .stratum[.stratum[, "code"] %in% .miss, "code"]]
        if (length(.miss) > 0)
          stop("Stratum(strata) ", .quot.past(.miss, quot = ""), " in ",
               "'tree.list.tls' argument is(are) missing in 'plot.parameters'",
               " argument.")
        
      }
      
      # Define code for unique pairs ('stratum', 'id')
      .stratum.id <- unique(tree.list.tls[, c("stratum", "id"), drop = FALSE])
      .stratum.id <- cbind(code = paste("(", apply(.stratum.id, 1, .quot.past),
                                        ")", sep = ""),
                           .stratum.id)
      
      # Check each plot ID has an unique associated stratum
      .dupl <- unique(.stratum.id[duplicated(.stratum.id[, "id"]), "id"])
      if (length(.dupl) > 0)
        stop("Plot ID(s) ", .quot.past(.dupl), " has(have) several different ",
             "strata associated in 'tree.list.tls' argument.")
      
      # Save code for pairs ('stratum', 'id'), and remove from 'tree.list.tls'
      .stratum.id <- .stratum.id[match(.id[, "id"], .stratum.id[, "id"]), ,
                                 drop = FALSE]
      .id <- cbind(.id, stratum.id.code = .stratum.id[, "code"],
                   stratum = .stratum.id[, "stratum"])
      tree.list.tls <- tree.list.tls[, !colnames(tree.list.tls) %in% "stratum",
                                     drop = FALSE]
      
    }
    
    # Check TXT files
    if ("file" %in% colnames(tree.list.tls)) {
      
      # Define code for unique pairs ('id', 'file')
      .id.file <- unique(tree.list.tls[, c("id", "file"), drop = FALSE])
      .id.file <- cbind(code = paste("(", apply(.id.file, 1, .quot.past), ")",
                                     sep = ""),
                        .id.file)
      
      # Check each plot ID has an unique associated TXT file
      .dupl <- unique(.id.file[duplicated(.id.file[, "id"]), "id"])
      if (length(.dupl) > 0)
        stop("Plot ID(s) ", .quot.past(.dupl), " has(have) several different ",
             "TXT files associated in 'tree.list.tls' argument.")
      
      # Check each TXT file has an unique associated plot ID
      .dupl <- unique(.id.file[duplicated(.id.file[, "file"]), "file"])
      if (length(.dupl) > 0)
        stop("TXT file(s) ", .quot.past(.dupl), " has(have) several different ",
             "plot IDs associated in 'tree.list.tls' argument.")
      
      # Check TXT files existence
      .miss <- .id.file[, "file"]
      .miss <- .miss[!file.exists(file.path(dir.data, .miss))]
      if (length(.miss) > 0)
        stop("TXT file(s) ", .quot.past(.miss), " in 'tree.list.tls' argument ",
             "does(do) not exists in the directory specified in 'dir.data' ",
             "argument.")
      
      # Save code for pairs ('id', 'file'), and remove file from 'tree.list.tls'
      .id.file <- .id.file[match(.id[, "id"], .id.file[, "id"]), , drop = FALSE]
      .id <- cbind(.id, id.file.code = .id.file[, "code"],
                   file = .id.file[, "file"])
      tree.list.tls <- tree.list.tls[, !colnames(tree.list.tls) %in% "file",
                                     drop = FALSE]
      
    }
    
    # Define code for pairs ('id', 'tree') and check they are not duplicated
    .id.tree <- tree.list.tls[, c("id", "tree"), drop = FALSE]
    .id.tree <- cbind(code = paste("(", apply(.id.tree, 1, .quot.past), ")",
                                   sep = ""),
                      .id.tree)
    .dupl <- unique(.id.tree[duplicated(.id.tree[, "code"]), "code"])
    if (length(.dupl) > 0)
      stop("'tree.list.tls' argument has several different rows for the ",
           "pair(s) of plot ID and tree ", .quot.past(.dupl, quot = ""), ".")
    
    # Check 'dbh, 'h', 'num.points', 'num.points.hom', 'num.points.est' and/or
    # 'num.points.hom.est' are strictly positive
    .inval <- c("dbh", "h", paste("num.points",
                                  c("", ".est", ".hom", ".hom.est"), sep = ""))
    .inval <- colnames(tree.list.tls)[colnames(tree.list.tls) %in% .inval]
    .inval <- .inval[!apply(tree.list.tls[, .inval, drop = FALSE] > 0, 2, all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.list.tls' argument must be strictly positive.")
    
    # Check 'num.points.hom'/'num.points.hom.est' is greater or equal to
    # 'num.points'/'num.points.est'
    .inval <- t(sapply(paste("num.points", c("", ".hom"), sep = ""), paste,
                       c("", ".est"), sep = ""))
    .inval <- .inval[, apply(apply(.inval, 1:2, "%in%",
                                   colnames(tree.list.tls)), 2, all),
                     drop = FALSE]
    if (ncol(.inval) > 0) {
      
      .inval <- .inval[, apply(.inval, 2,
                               function(dat.par, dat) {
                                 any(dat[, dat.par[1]] < dat[, dat.par[2]])
                               },
                               dat = tree.list.tls),
                       drop = FALSE]
      if (ncol(.inval) > 0)
        stop(paste("All values in column(s) ", .quot.past(.inval[1, ]), " of ",
                   "'tree.list.tls' argument must be greater or equal to the ",
                   "corresponding ones in column(s) ", .quot.past(.inval[2, ]),
                   ". ", sep = ""))
      
    }
    
    # Check 'h.dist' is positive
    .inval <- "h.dist"
    .inval <- colnames(tree.list.tls)[colnames(tree.list.tls) %in% .inval]
    .inval <- .inval[!apply(tree.list.tls[, .inval, drop = FALSE] >= 0, 2, all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.list.tls' argument must be positive.")
    
    # Check 'phi.left' and 'phi.right' belong to the interval [0, 2 * pi]
    .inval <- unlist(.scan.approach)
    .inval <- colnames(tree.list.tls)[colnames(tree.list.tls) %in% .inval]
    .inval <- .inval[!apply(tree.list.tls[, .inval, drop = FALSE] >= 0 &
                              tree.list.tls[, .inval, drop = FALSE] <= 2 * pi,
                            2, all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.list.tls' argument must be in the interval [0, 2 * pi].")
    
    # Check 'partial.occlusion' is 0 or 1
    .inval <- "partial.occlusion"
    .inval <- colnames(tree.list.tls)[colnames(tree.list.tls) %in% .inval]
    .inval <- .inval[!apply(tree.list.tls[, .inval, drop = FALSE] == 0 |
                              tree.list.tls[, .inval, drop = FALSE] == 1,
                            2, all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.list.tls' argument must be 0 or 1.")
    
  }
  
  
  # Check 'tree.list.ds' argument ----
  
  if (!"tls" %in% names(.funct[[funct]])) tree.list.ds <- NULL
  else if (!is.null(tree.list.ds)) {
    
    # Define a logical matrix containing all possible mandatory columns in
    # 'tree.list.ds' according to 'var.metr'
    .col.mand <- c("id", "tree", .ds.meth)
    .col.mand <- matrix(FALSE, nrow = length(var.metr$tls),
                        ncol = length(.col.mand),
                        dimnames = list(var.metr$tls, .col.mand))
    for (.i in names(.ds.meth))
      .col.mand[rownames(.col.mand) %in% paste(c("N", "G", "V"), .i, sep = "."),
                colnames(.col.mand) %in%
                  c("id", "tree", paste("P", .i, sep = "."))] <- TRUE
    .col.mand <- .col.mand[, apply(.col.mand, 2, any), drop = FALSE]
    
    # Define a character vector containing (if necessary) the class of all
    # possible mandatory columns in 'tree.list.ds'
    .col.class <- rep("numeric", ncol(.col.mand))
    names(.col.class) <- colnames(.col.mand)
    .col.class[names(.col.class) %in% c("id", "tree")] <- NA
    .col.class <- .col.class[!is.na(.col.class)]
    
    # If there are no mandatory columns, 'tree.list.ds' is forced to be NULL;
    # otherwise, checking process is continued
    if (ncol(.col.mand) == 0) tree.list.ds <- NULL
    else {
      
      # 'tree.list.ds' must be a list with at least an element named 'tree'
      tree.list.ds <- .check.class(x = tree.list.ds, x.class = "list",
                                   name = "tree.list.ds", n = Inf)
      .miss <- c("tree")[!c("tree") %in% names(tree.list.ds)]
      if (length(.miss) > 0)
        stop("Element(s) named ", .quot.past(.miss), " is(are) missing in ",
             "'tree.list.ds' argument.")
      tree.list.ds <- tree.list.ds$tree
      
      # 'tree.list.ds$tree' must be a data.frame with at least one row
      tree.list.ds <- .check.class(x = tree.list.ds, x.class = "data.frame",
                                   name = "tree.list.ds$tree", n = Inf)
      
      # Check mandatory columns existence, remove non necessary columns, check
      # mandatory columns classes and check all values are no-NA
      tree.list.ds <- .check.col(x = tree.list.ds, x.mand = .col.mand,
                                 x.class = .col.class,
                                 name = "tree.list.ds$tree")
      
      # Discard rows corresponding to plot IDs not included in 'tree.list.tls'
      .miss <- unique(tree.list.ds[!tree.list.ds[, "id"] %in% .id[, "id"],
                                   "id"])
      if (length(.miss) > 0) {
        
        warning("Row(s) of 'tree.list.ds$tree' argument associated to plot ",
                "ID(s) ", .quot.past(.miss, quot = ""), " was(were) discarded ",
                "because this(these) plot(s) is(are) missing in ",
                "'tree.list.tls' argument.", immediate. = TRUE)
        tree.list.ds <- tree.list.ds[!tree.list.ds[, "id"] %in% .miss, ,
                                     drop = FALSE]
        
      }
      
      # Check plot IDs included in 'tree.list.tls' are no missing
      .miss <- .id[, "id"]
      .miss <- .miss[!.miss %in% unique(tree.list.ds[, "id"])]
      if (length(.miss) > 0)
        stop("Plot ID(s) ", .quot.past(.miss), " in 'tree.list.tls' ",
             "argument is(are) missing in 'tree.list.ds$tree' argument.")
      
      # Check there are no duplicated pairs ('id', 'tree')
      .dupl <- tree.list.ds[, c("id", "tree"), drop = FALSE]
      .dupl <- paste("(", apply(.dupl, 1, .quot.past), ")", sep = "")
      .dupl <- unique(.dupl[duplicated(.dupl)])
      if (length(.dupl) > 0)
        stop("'tree.list.ds$tree' argument has several different rows for the ",
             "pair(s) of plot ID and tree ", .quot.past(.dupl, quot = ""), ".")
      
      # Check pairs ('id', 'tree') in 'tree.list.tls' are included in
      # 'tree.list.ds'
      .miss <- tree.list.ds[, c("id", "tree"), drop = FALSE]
      .miss <- paste("(", apply(.miss, 1, .quot.past), ")", sep = "")
      .miss <- .id.tree[, "code"][!.id.tree[, "code"] %in% .miss]
      if (length(.miss) > 0)
        stop("Pair(s) of plot ID and tree ", .quot.past(.miss, quot = ""),
             " in 'tree.list.tls' argument is(are) missing in ",
             "'tree.list.ds$tree' argument.")
      
      # Check 'P.hn', 'P.hr', 'P.hn.cov' and 'P.hr.cov' belong to the interval
      # [0, 1]
      .inval <- .ds.meth
      .inval <- colnames(tree.list.ds)[colnames(tree.list.ds) %in% .inval]
      .inval <- .inval[!apply(tree.list.ds[, .inval, drop = FALSE] >= 0 &
                                tree.list.ds[, .inval, drop = FALSE] <= 1,
                              2, all)]
      if (length(.inval) > 0)
        stop("All values in column(s) ", .quot.past(.inval), " of ",
             "'tree.list.ds$tree' argument must be in the interval [0, 1].")
      
    }
    
  }
  
  
  # Check 'tree.list.field' argument ----
  
  if (!"field" %in% names(.funct[[funct]])) tree.list.field <- NULL
  else {
    
    # Define a logical matrix containing all possible mandatory columns in
    # 'tree.list.field' according to 'var.metr'
    .col.mand <- c("id", "tree", "h.dist", "dbh", "h", unlist(.var.field.user))
    .col.mand <- matrix(FALSE, nrow = length(var.metr$field),
                        ncol = length(.col.mand),
                        dimnames = list(var.metr$field, .col.mand))
    .col.mand[, colnames(.col.mand) %in%
                c("id", "tree", "h.dist", "dbh", "h")] <- TRUE
    for (.i in names(.var.field.user))
      .col.mand[rownames(.col.mand) %in% .i,
                colnames(.col.mand) %in% .var.field.user[[.i]]] <- TRUE
    .col.mand <- .col.mand[, apply(.col.mand, 2, any), drop = FALSE]
    
    # Define a character vector containing (if necessary) the class of all
    # possible mandatory columns in 'tree.list.field'
    .col.class <- rep("numeric", ncol(.col.mand))
    names(.col.class) <- colnames(.col.mand)
    .col.class[names(.col.class) %in% c("id", "tree")] <- NA
    .col.class <- .col.class[!is.na(.col.class)]
    
    # 'tree.list.field' must be a data.frame with at least one row
    tree.list.field <- .check.class(x = tree.list.field, x.class = "data.frame",
                                    name = "tree.list.field", n = Inf)
    
    # Check mandatory columns existence, remove non necessary columns, check
    # mandatory columns classes and check all values are no-NA
    tree.list.field <- .check.col(x = tree.list.field, x.mand = .col.mand,
                                  x.class = .col.class,
                                  name = "tree.list.field")
    
    # Discard rows corresponding to plot IDs not included in 'tree.list.tls'
    .miss <- unique(tree.list.field[!tree.list.field[, "id"] %in% .id[, "id"],
                                    "id"])
    if (length(.miss) > 0) {
      
      warning("Row(s) of 'tree.list.field' argument associated to plot ID(s) ",
              .quot.past(.miss, quot = ""), " was(were) discarded because ",
              "this(these) plot(s) is(are) missing in 'tree.list.tls' ",
              "argument.", immediate. = TRUE)
      tree.list.field <- tree.list.field[!tree.list.field[, "id"] %in% .miss, ,
                                         drop = FALSE]
      
    }
    
    # Check plot IDs included in 'tree.list.tls' are no missing
    .miss <- .id[, "id"]
    .miss <- .miss[!.miss %in% unique(tree.list.field[, "id"])]
    if (length(.miss) > 0)
      stop("Plot ID(s) ", .quot.past(.miss), " in 'tree.list.tls' ",
           "argument is(are) missing in 'tree.list.field' argument.")
    
    # Check there are no duplicated pairs ('id', 'tree')
    .dupl <- tree.list.field[, c("id", "tree"), drop = FALSE]
    .dupl <- paste("(", apply(.dupl, 1, .quot.past), ")", sep = "")
    .dupl <- unique(.dupl[duplicated(.dupl)])
    if (length(.dupl) > 0)
      stop("'tree.list.field' argument has several different rows for the ",
           "pair(s) of plot ID and tree ", .quot.past(.dupl, quot = ""), ".")
    
    # Check 'dbh, 'h', 'v.user' and/or 'w.user' are strictly positive
    .inval <- c("dbh", "h", unlist(.var.field.user))
    .inval <- colnames(tree.list.field)[colnames(tree.list.field) %in% .inval]
    .inval <- .inval[!apply(tree.list.field[, .inval, drop = FALSE] > 0, 2,
                            all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.list.field' argument must be strictly positive.")
    
    # Check 'h.dist' is positive
    .inval <- "h.dist"
    .inval <- colnames(tree.list.field)[colnames(tree.list.field) %in% .inval]
    .inval <- .inval[!apply(tree.list.field[, .inval, drop = FALSE] >= 0, 2,
                            all)]
    if (length(.inval) > 0)
      stop("All values in column(s) ", .quot.past(.inval), " of ",
           "'tree.list.field' argument must be positive.")
    
  }
  
  
  # Create a list containing empty data.frames where results will be saved ----
  # and restrict trees' database(s) according to 'dbh.min', 'h.min' and/or
  # 'max.dist'
  
  # Create a list containing empty data.frames where results will be saved for
  # each plot design
  stand <- vector("list", length(plot.design))
  names(stand) <- plot.design
  for (.i in names(stand)) {
    
    .col.names <- c("id", .plot.design[.i])
    for (.j in names(.funct[[funct]]))
      .col.names <- c(.col.names,
                      var.metr[[.j]][.var.metr[[.j]][.i, var.metr[[.j]]]])
    if (.by.stratum.2) .col.names <- c("stratum", .col.names)
    stand[[.i]] <- data.frame(matrix(numeric(), ncol = length(.col.names),
                                     dimnames = list(NULL, .col.names)),
                              stringsAsFactors = FALSE)
    
  }
  
  # Restrict trees' database(s) according to 'dbh.min', 'h.min' and/or
  # 'max.dist'
  for (.i in names(.funct[[funct]])) {
    
    .tree <- get(paste("tree.list", .i, sep ="."))
    .inval <- which(.tree[, "dbh"] < dbh.min | .tree[, "h"] < h.min |
                      .tree[, "h.dist"] > max.dist)
    if (length(.inval) > 0) {
      
      warning(length(.inval), " tree(s) in 'tree.list.", .i, "' argument ",
              "was(were) discarded, since dbh or height were less than minima ",
              "values indicated in 'dbh.min' and 'h.min' arguments, or ",
              "horizontal distance is greater than 'max.dist' argument.",
              immediate. = TRUE)
      .tree <- .tree[- .inval, , drop = FALSE]
      assign(paste("tree.list", .i, sep ="."), .tree)
      
      # Remove plot ID(s) without trees for calculations below
      .miss <- .id[!.id[, "id"] %in% unique(.tree[, "id"]), "id"]
      if (length(.miss) > 0) {
        
        if (length(.miss) < nrow(.id)) {
          
          warning("Plot ID(s) ", .quot.past(.miss), " in 'tree.list.", .i,
                  "' argument was(were) discarded, since all its(their) trees ",
                  "have dbh or height less than minima values indicated in ",
                  "'dbh.min' and 'h.min' arguments, or horizontal distance ",
                  "greater than 'max.dist' argument.", immediate. = TRUE)
          .id <- .id[!.id[, "id"] %in% .miss, , drop = FALSE]
          
        } else
          stop("All plot ID(s) in 'tree.list.", .i, "' argument was(were) ",
               "discarded, since all trees have dbh or height less than ",
               "minima values indicated in 'dbh.min' and 'h.min' arguments, ",
               "or horizontal distance greater than 'max.dist' argument.")
        
      }
      
    }
    
  }
  
  
  # Loop for each TLS plot ----
  
  for (.i in 1:nrow(.id)) {
    
    # Define initial time, and print message with the plot ID
    t0 <- Sys.time()
    message("Computing ",
            ifelse(funct %in% c("sim", "est"), "simulations", "metrics"),
            " for plot: ", .id[.i, "code"])
    
    # Define list where results associated to the plot will be saved
    .stand <- vector("list", length(.funct[[funct]]))
    names(.stand) <- names(.funct[[funct]])
    
    # Select plot parameters according to '.by.stratum'
    .str <- ifelse(.by.stratum,
                   match(.id[.i, "stratum"], plot.parameters[, "stratum"]), 1)
    .par <- plot.parameters[.str, , drop = FALSE]
    
    
    # Loop for each trees' database (TLS and field data) ----
    
    for (.j in names(.stand)) {
      
      # Define list where results associated to the plot will be saved
      .stand[[.j]] <- vector("list", length(plot.design))
      names(.stand[[.j]]) <- plot.design
      
      
      # Create trees' database (TLS and field data) ----
      
      # Select data corresponding to the plot from the complete trees' database
      .tree <- get(paste("tree.list", .j, sep ="."))
      .tree <- .tree[.tree[, "id"] %in% .id[.i, "id"], , drop = FALSE]
      rownames(.tree) <- NULL
      
      # Convert dbh (cm) to International System of Units (m)
      .tree[, "dbh"] <- .tree[, "dbh"] / 100
      
      # Compute radius, k and BAF, and tree variables according to plot
      # design(s) and 'tree.var'. Currently available tree variables: basal area
      # (g) and volume (v)
      .col.mand <- c("g", "v")
      .col.mand <- matrix(FALSE, nrow = length(var.metr[[.j]]),
                          ncol = length(.col.mand),
                          dimnames = list(var.metr[[.j]], .col.mand))
      .col.mand[rownames(.col.mand) %in%
                  c(paste("G", c("tls", names(.ds.meth), "sh", "pam"),
                          sep = "."), "G"),
                colnames(.col.mand) %in% "g"] <- TRUE
      .col.mand[rownames(.col.mand) %in%
                  c(paste("V", c("tls", names(.ds.meth), "sh", "pam"),
                          sep = "."), "V"),
                colnames(.col.mand) %in% "v"] <- TRUE
      .col.mand <- colnames(.col.mand)[apply(.col.mand, 2, any)]
      .tree <- .tree.calc(tree = .tree, plot.design = .plot.design[plot.design],
                          tree.var = .col.mand, v.calc = v.calc)
      
      # Compute angular aperture (TLS data)
      if (.j == "tls" & all(c("phi.left", "phi.right") %in% colnames(.tree))) {
        
        # Compute angular aperture
        .wide <- .tree[, "phi.right"] - .tree[, "phi.left"]
        .wide <- ifelse(.wide < 0, (2 * pi) + .wide, .wide)
        
        # Remove angular coordinates, and add angular aperture
        .tree <- .tree[, !colnames(.tree) %in% c("phi.left", "phi.right"),
                       drop = FALSE]
        .tree <- cbind(.tree, wide = .wide)
        
      }
      
      # Add detection probability from distance sampling database (TLS data)
      if (.j == "tls" & !is.null(tree.list.ds)) {
        
        .ds <- tree.list.ds[tree.list.ds[, "id"] %in% .id[.i, "id"], ,
                            drop = FALSE]
        .tree <- cbind(.tree, .ds[match(.tree[, "tree"], .ds[, "tree"]),
                                  colnames(.ds) %in% .ds.meth, drop = FALSE])
        
      }
      
      
      # Create points database(s) (TLS data) ----
      
      .data.tls <- NULL
      if (.j == "tls" & "file" %in% colnames(.id)) {
        
        # Read points' database
        .data.tls <-
          suppressMessages(vroom::vroom(file.path(dir.data, .id[.i, "file"]),
                                        col_select = c("x", "y", "z", "rho"),
                                        progress = FALSE))
        
        # Discard points according to 'max.dist'
        .inval <- which(.data.tls[, "rho"] > max.dist)
        if (length(.inval) > 0) {
          
          warning(length(.inval), " point(s) in point cloud associated to ",
                  "plot ", .id[.i, "code"], " was(were) discarded, since ",
                  "horizontal distance is greater than 'max.dist' argument.",
                  immediate. = TRUE)
          .data.tls <- .data.tls[- .inval, , drop = FALSE]
          
        }
        
        # Order .data.tls by rho, select columns required for calculations
        # below, and convert to matrix
        .data.tls <- as.matrix(.data.tls)[order(.data.tls[, "rho"],
                                                decreasing = FALSE),
                                          c("z", "rho"), drop = FALSE]
        
      }
      
      
      # Loop for each plot design (TLS and field data) ----
      
      for (.k in plot.design) {
        
        # Order trees' database by radius/k/BAF in order to simplify
        # calculations below
        if (.k %in% c("fixed.area", "k.tree"))
          .tree <- .tree[order(.tree[, .plot.design[.k]], decreasing = FALSE), ,
                         drop = FALSE]
        else if (.k %in% "angle.count")
          .tree <- .tree[order(.tree[, .plot.design[.k]], decreasing = TRUE), ,
                         drop = FALSE]
        
        
        # Create radius/k/BAF sequence according to 'function case', ----
        # 'plot.parameters' argument, and radius/k/BAF values in trees' database
        
        # Select number of decimals
        .par.dec <- .par[, paste(.plot.design[.k], "dec", sep = ".")]
        
        # Compute range of values in trees' database
        .funct.rang <- switch(.k, fixed.area = ".customCeiling",
                              k.tree = ".customCeiling",
                              angle.count = ".customFloor")
        .db.rang <- range(get(.funct.rang)(.tree[, .plot.design[.k]],
                                           Decimals = .par.dec))
        
        # Select minimum value (ensuring at least one tree is always selected
        # for fixed area and k-tree plot designs)
        .par.min <- .db.rang[1]
        
        # Select maximum value
        if (funct %in% c("sim", "est"))
          .par.max <- .par[, paste(.plot.design[.k], "max", sep = ".")]
        else .par.max <- .par[, .plot.design[.k]]
        # Adjust maximum value according to 'max.dist' for fixed area plot
        # design, number of trees for k-tree plot design, and maximum BAF for
        # angle-count plot design
        .ref.max <- switch(.k, fixed.area = max.dist, k.tree = .db.rang[2],
                           angle.count = .db.rang[2])
        if (.par.max > .ref.max) {
          
          .par.max <- .ref.max
          warning("For ", .funct[[funct]][.j], " plot ", .id[.i, "code"], ", ",
                  ifelse(funct %in% c("sim", "est"), "maximum ", ""),
                  .quot.past(.plot.design[.k]), " was reduced to ",
                  .par.max,
                  switch(.k,
                         fixed.area = paste(", since this is the maximum",
                                            "horizontal distance indicated in",
                                            "'max.dist' argument"),
                         k.tree = paste(", since this is the number of trees",
                                        "in the plot"),
                         angle.count = paste(" to ensure that at least one",
                                             "tree is always selected")),
                  ".", immediate. = TRUE)
          
        }
        # Adjust maximum value according to minimum one
        if (.par.min > .par.max) {
          
          .par.max <- .par.min
          warning("For ", .funct[[funct]][.j], " plot ", .id[.i, "code"], ", ",
                  ifelse(funct %in% c("sim", "est"), "maximum ", ""),
                  .quot.past(.plot.design[.k]), " was increased to ", .par.max,
                  switch(.k,
                         fixed.area = paste(" to ensure that at least one",
                                            "tree is always selected"),
                         k.tree = paste(" to ensure that at least one",
                                        "tree is always selected"),
                         angle.count = paste(", since this ensures all trees",
                                             "in the plot are selected")),
                  ".", immediate. = TRUE)
          
        }
        
        # Force coincidence between minimum and maximum value, and select
        # increment according to 'function case'
        if (!funct %in% c("sim", "est")) {
          
          .par.min <- .par.max
          .par.incr <- 0
          
        } else .par.incr <- .par[, paste(.plot.design[.k], "incr", sep = ".")]
        
        # Create sequence
        .par.seq <- seq(from = .par.max, to = .par.min, by = - .par.incr)
        .par.seq <- sort(unique(round(.par.seq, .par.dec)))
        names(.par.seq) <- .format.numb(x = .par.seq, dec = .par.dec)
        
        
        # Compute stand variables and/or metrics according to plot design  ----
        # and 'var.metr'
        
        # Compute expansion factors with occlusion corrections and estimate
        # stand variables per ha, and compute mean diameters and heights
        .col.names <- var.metr[[.j]][.var.metr[[.j]][.k, var.metr[[.j]]]]
        .col.names <-
          .col.names[.col.names %in%
                       c(sapply(c("N", "G", "V"), paste,
                                c("tls", names(.ds.meth), "sh", "pam"),
                                sep = "."),
                         paste(names(c(.mean.d, .mean.h)), "tls", sep = "."),
                         paste("num.points", c("", ".est", ".hom", ".hom.est"),
                               sep = ""),
                         "N", "G", "V", names(.var.field.user),
                         names(c(.mean.d, .mean.h)))]
        .stand[[.j]][[.k]] <- lapply(.par.seq, .stand.calc,
                                     plot.design = .plot.design[.k],
                                     tree = .tree, var.metr = .col.names,
                                     ds.meth = .ds.meth,
                                     var.field.user = .var.field.user,
                                     mean.d = .mean.d, mean.h = .mean.h)
        .stand[[.j]][[.k]] <- do.call(rbind, .stand[[.j]][[.k]])
        
        # Compute mean dominant diameters and heights
        .col.names <- paste(names(c(.mean.d, .mean.h)), ".0", sep = "")
        if (.j %in% "tls") .col.names <- paste(.col.names, .j, sep = ".")
        names(.col.names) <- paste(names(c(.mean.d, .mean.h)), ".0", sep = "")
        .col.names <-
          .col.names[.col.names %in%
                       var.metr[[.j]][.var.metr[[.j]][.k, var.metr[[.j]]]]]
        if (length(.col.names) > 0) {
          
          if (.k %in% "fixed.area")
            .dh.0 <- fixed_area_cpp(radius_seq = .par.seq,
                                    hdist = .tree[, "h.dist"],
                                    d = .tree[, "dbh"], h = .tree[, "h"],
                                    num = .par[, "num.trees"])
          else if (.k %in% "k.tree")
            .dh.0 <- k_tree_cpp(k_seq = .par.seq,
                                radius_seq = .tree[, "radius.k"],
                                k = .tree[, "k"], d = .tree[, "dbh"],
                                h = .tree[, "h"], num = .par[, "num.trees"])
          else if (.k %in% "angle.count")
            .dh.0 <- angle_count_cpp(baf_seq = .par.seq, baf = .tree[, "BAF"],
                                     d = .tree[, "dbh"], h = .tree[, "h"],
                                     num = .par[, "num.trees"])
          rownames(.dh.0) <- .format.numb(x = .dh.0[, .plot.design[.k]],
                                          dec = .par.dec)
          colnames(.dh.0)[colnames(.dh.0) %in% names(.col.names)] <-
            .col.names[colnames(.dh.0)[colnames(.dh.0) %in% names(.col.names)]]
          .stand[[.j]][[.k]] <- cbind(.stand[[.j]][[.k]],
                                      .dh.0[rownames(.stand[[.j]][[.k]]),
                                            .col.names, drop = FALSE])
          
        }
        
        # Compute rho sequence for z coordinate metrics
        if (.k %in% "fixed.area") .rho.seq <- .par.seq
        else if (.k %in% "k.tree")
          .rho.seq <- sapply(.par.seq,
                             function(k, tree) {
                               max(tree[tree[, "k"] <= k, "radius.k"])
                             },
                             tree = .tree)
        else if (.k %in% "angle.count")
          .rho.seq <- sapply(.par.seq,
                             function(BAF, tree) {
                               max(tree[tree[, "BAF"] >= BAF, "h.dist"])
                             },
                             tree = .tree)
        
        # Compute percentiles of z coordinate
        .col.names <- sprintf("P%02i", .prob)
        .col.names <-
          .col.names[.col.names %in%
                       var.metr[[.j]][.var.metr[[.j]][.k, var.metr[[.j]]]]]
        if (length(.col.names) > 0) {
          
          .perc <- height_perc_cpp(rho_seq = .rho.seq, z = .data.tls[, "z"],
                                   rho = .data.tls[, "rho"])
          rownames(.perc) <- names(.rho.seq)
          .stand[[.j]][[.k]] <- cbind(.stand[[.j]][[.k]],
                                      .perc[rownames(.stand[[.j]][[.k]]),
                                            .col.names, drop = FALSE])
          
        }
        
        # Compute descriptive statistics of z coordinate, percentage of points
        # above mode and mean, and parameters of fitted Weibull distribution
        .col.names <- c("mean", "max", "min", "sd", "var", "mode", "kurtosis",
                        "skewness", "perc_on_mode", "perc_on_mean", "weibull_c",
                        "weibull_b")
        .col.names <-
          .col.names[.col.names %in%
                       var.metr[[.j]][.var.metr[[.j]][.k, var.metr[[.j]]]]]
        if (length(.col.names) > 0) {
          
          .pts.met <- .points.metrics(rho_seq = .rho.seq, data = .data.tls,
                                      metr = .col.names)
          .stand[[.j]][[.k]] <- cbind(.stand[[.j]][[.k]],
                                      .pts.met[rownames(.stand[[.j]][[.k]]),
                                               .col.names, drop = FALSE])

        }
        
        
        # Convert diameters from International System of Units (m) to cm
        .col.names <- sapply(names(.mean.d), paste,
                             c(".tls", ".0.tls", "", ".0"), sep = "")
        .col.names <-
          colnames(.stand[[.j]][[.k]])[colnames(.stand[[.j]][[.k]]) %in%
                                         .col.names]
        .stand[[.j]][[.k]][.col.names] <- .stand[[.j]][[.k]][.col.names] * 100
        
        # Add stratum
        if (.by.stratum.2)
          .stand[[.j]][[.k]] <- cbind(stratum = .id[.i, "stratum"],
                                      .stand[[.j]][[.k]])
        
      }
      
    }
    
    
    # Save fixed area/k-tree/angle-count plot results, and write csv files ----
    # containing them if 'save.result' is TRUE
    
    for (.j in plot.design) {
      
      
      # Merge fixed area/k-tree/angle-count plot results for all trees'
      # database(s)
      
      # Results for first trees' database(s)
      .stand.j <- .stand[[1]][[.j]]
      
      # Add results for the rest of trees' databases (if any)
      for (.k in names(.stand)[-1]) {
        
        # Restrict rows to those existing in trees' database to be added
        .stand.j <- .stand.j[rownames(.stand.j) %in% rownames(.stand[[.k]][[.j]]), ,
                             drop = FALSE]
        
        # Add new columns
        .stand.j <- cbind(.stand.j,
                          .stand[[.k]][[.j]][rownames(.stand.j),
                                             !colnames(.stand[[.k]][[.j]]) %in%
                                               colnames(.stand.j),
                                             drop = FALSE])
      }
      
      # Save reordered merged results
      rownames(.stand.j) <- NULL
      stand[[.j]] <- rbind(stand[[.j]],
                           .stand.j[, colnames(stand[[.j]]), drop = FALSE])
      stand[[.j]] <- stand[[.j]][order(stand[[.j]][, .plot.design[[.j]]],
                                       stand[[.j]][, "id"]), ,
                                 drop = FALSE]
      rownames(stand[[.j]]) <- NULL
      
      # Write CSV file
      if (save.result) {
        
        .file <- paste(switch(funct, sim = "simulations",
                              metr = "metrics.variables",
                              est = "estimation.plot.size"),
                       .j, "plot.csv", sep = ".")
        utils::write.csv(stand[[.j]], file = file.path(dir.result, .file),
                         row.names = FALSE)
        
      }
      
    }
    
    # Define final time, and print message
    t1 <- Sys.time()
    message(" (", format(round(difftime(t1, t0, units = "secs"), 2)), ")")
    
  }
  
  return(stand)
  
}


# Compute Pearson/Spearman correlations

.cor.pearson <- function(x, y){
  
  stats::cor.test(x = x, y = y, method = 'pearson')
  
}

.cor.spearman <- function(x, y){
  
  stats::cor.test(x = x, y = y, method = 'spearman', exact = FALSE)
  
}


##############################################################################
# Cálculo del índice de curvatura
##############################################################################


# .pol create polygons as regular squares overlapped 0.05

.pol <- function(data, dist){
  
  .xpol <- c(data$x.1-dist, data$x.2+dist, data$x.3+dist, data$x.4-dist)
  .ypol <- c(data$y.1-dist, data$y.2-dist, data$y.3+dist, data$y.4+dist)
  
  return(sp::Polygons(list(sp::Polygon(cbind(.xpol,.ypol))), ID = data$id))
  
}


# This functions apply calculated ncr index for all points.
# For that purpose, voxelize point cloud by means of regular
# grid in x any z coordinates

.ncr.remove.slice.double <- function(data){
  
  # Select necessary fields from original txt file of point cloud
  
  .data <- data[,c("point", "x", "y", "z")]
  
  
  # Create x and y coordinates for grid
  
  .x = seq(min(.data$x), max(.data$x))
  .y = seq(min(.data$y), max(.data$y))
  
  # Empty data frame where coordinates neccesaries for
  # creating grid will be saved
  
  .grid <- data.frame(id = as.character(),
                      x.1 = as.numeric(), x.2 = as.numeric(),
                      x.3 = as.numeric(), x.4 = as.numeric(),
                      y.1 = as.numeric(), y.2 = as.numeric(),
                      y.3 = as.numeric(), y.4 = as.numeric())
  
  
  # Fill the empty data frame .grid created just before
  
  for (i in 1:(length(.x)-1)) {
    for (j in 1:(length(.y)-1)) {
      
      .out <- data.frame(id = paste("x", i, "y", j, sep="."),
                         x.1 = .x[i], x.2 = .x[i+1], x.3 = .x[i+1], x.4 = .x[i],
                         y.1 = .y[j], y.2 = .y[j], y.3 = .y[j+1], y.4 = .y[j+1])
      
      .grid <- rbind(.grid, .out)
      
    }
  }
  
  .row.names <- .grid$id
  
  # Split .grid by id into a list for using lapply function and
  # applying .pol function over all elements of the list
  
  .grid <- split(.grid, .grid$id)
  
  # Apply .pol function to every element of the list and create
  # and object SpatialPolygons
  
  .grid <- sp::SpatialPolygons(lapply(.grid, .pol, dist = 0.05))
  
  
  # Generate and SpatialPoints object to extract those points
  # which are avor the polygons of .grid
  
  .pts <- .data[, c("x", "y")]
  # dimnames(.pts)[[1]] <- c(.data$point)
  .pts <- sp::SpatialPoints(.pts)
  
  .attributes <- .data[, c("point", "x", "y", "z")]
  
  .pts = sp::SpatialPointsDataFrame(.pts, .attributes)
  
  .attributes <- data.frame(row.names = .row.names)
  
  .grid = sp::SpatialPolygonsDataFrame(.grid, .attributes)
  
  .pts <- sp::over(.grid, .pts, returnList = TRUE)
  
  .dat <- lapply(.pts, as.matrix, ncol = 4)
  .dat <- .dat[names(which(lapply(.dat, length) > 4))]
  # .dat <- .dat[names(which(lapply(.dat, length) < 40000))]
  
  
  .ncr <- do.call(rbind, lapply(.dat, ncr_point_cloud_double))
  
  .ncr <- .ncr[which(.ncr$ncr > 0 & .ncr$ncr < 9999), ]
  
  .data <- merge(data, .ncr, by = "point", all = TRUE)
  
  .data <- .data[!duplicated(.data), ]
  
  return(.data)
  
}
