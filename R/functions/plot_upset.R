plot_upset <- function(gene_lists,
                       order = c("0.25", "0.5", "1", "2", "12", "33", "63"),
                       order_by = "freq",
                       keep_order = TRUE) {
  upset(
    fromList(gene_lists),
    sets = order,
    order.by = order_by,
    keep.order = keep_order
  )
}
