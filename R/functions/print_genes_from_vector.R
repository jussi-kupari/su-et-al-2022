print_gene_list_from_vector <- function(genes) {
  for (i in seq_along(genes)) {
    cat(genes[i], sep = "\n")
  }
}
