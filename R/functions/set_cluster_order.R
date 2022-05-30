set_cluster_order <- function(obj, order) {
  Idents(obj) <- factor(x = Idents(obj), levels = order)
  return(obj)
}