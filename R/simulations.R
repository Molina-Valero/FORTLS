
simulations <- function(tree.list.tls, tree.list.ds = NULL, tree.list.field,
                        plot.design = c("fixed.area", "k.tree", "angle.count"),
                        plot.parameters = data.frame(radius.max = 25,
                                                     k.max = 50, BAF.max = 4),
                        scan.approach = "single",
                        var.metr = list(tls = NULL, field = NULL),
                        v.calc = "parab", dbh.min = 7.5, h.min = 1.3,
                        max.dist = Inf, dir.data = NULL, save.result = TRUE,
                        dir.result = NULL) {
  
  # Call internal function
  sim <- .sim.calc(funct = "sim", tree.list.tls = tree.list.tls,
                   tree.list.ds = tree.list.ds,
                   tree.list.field = tree.list.field, plot.design = plot.design,
                   plot.parameters = plot.parameters,
                   scan.approach = scan.approach, var.metr = var.metr,
                   v.calc = v.calc, dbh.min = dbh.min, h.min = h.min,
                   max.dist = max.dist, dir.data = dir.data,
                   save.result = save.result, dir.result = dir.result)
  return(sim)
  
}
