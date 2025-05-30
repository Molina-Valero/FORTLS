
metrics.variables <- function(tree.tls, tree.ds = NULL, tree.field = NULL,
                              plot.design = c("fixed.area", "k.tree", "angle.count"),
                              plot.parameters = data.frame(radius = 25, k = 10, BAF = 2),
                              scan.approach = "single",
                              var.metr = list(tls = NULL, field = NULL),
                              v.calc = "parab", dbh.min = 4, h.min = 1.3,
                              max.dist = Inf, dir.data = NULL,
                              save.result = TRUE, dir.result = NULL) {

  if(!is.null(tree.field) & !is.null(var.metr$tls) & !is.null(var.metr$field))
    var.metr = list(tls = var.metr$tls, field = var.metr$field)

  if(!is.null(tree.field) & !is.null(var.metr$tls) & is.null(var.metr$field))
    var.metr = list(tls = var.metr$tls, field = NULL)

  if(!is.null(tree.field) & !is.null(var.metr$field) & is.null(var.metr$tls))
    var.metr = list(tls = NULL, field = var.metr$field)

  # Call internal function
  metr <- .sim.calc(funct = "metr", tree.tls = tree.tls,
                    tree.ds = tree.ds, tree.field = tree.field,
                    plot.design = plot.design,
                    plot.parameters = plot.parameters,
                    scan.approach = scan.approach,
                    var.metr = var.metr, v.calc = v.calc,
                    dbh.min = dbh.min, h.min = h.min, max.dist = max.dist,
                    dir.data = dir.data, save.result = save.result,
                    dir.result = dir.result)

  return(metr)

}
