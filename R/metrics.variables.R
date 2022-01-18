
metrics.variables <- function(tree.list.tls, tree.list.ds = NULL,
                              plot.design = c("fixed.area", "k.tree",
                                              "angle.count"),
                              plot.parameters = data.frame(radius = 25,
                                                           k = 50, BAF = 4),
                              scan.approach = "single", var.metr = NULL,
                              v.calc = "parab", dbh.min = 7.5, h.min = 1.3,
                              max.dist = Inf, dir.data = NULL,
                              save.result = TRUE, dir.result = NULL) {
  
  # Call internal function
  metr <- .sim.calc(funct = "metr", tree.list.tls = tree.list.tls,
                    tree.list.ds = tree.list.ds, tree.list.field = NULL,
                    plot.design = plot.design,
                    plot.parameters = plot.parameters,
                    scan.approach = scan.approach,
                    var.metr = list(tls = var.metr), v.calc = v.calc,
                    dbh.min = dbh.min, h.min = h.min, max.dist = max.dist,
                    dir.data = dir.data, save.result = save.result,
                    dir.result = dir.result)
  
  return(metr)
  
}
