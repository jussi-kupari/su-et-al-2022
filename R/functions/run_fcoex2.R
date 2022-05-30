run_fcoex2 <- function(data, n_genes) {
  
  exprs <- data.frame(GetAssayData(data, slot = "data"))
  target <- Idents(data)
  
  fc <- new_fcoex(data.frame(exprs), target)
  fc <- fc %>% 
    discretize() %>% 
    find_cbf_modules(
      n_genes_selected_in_first_step = n_genes, verbose = FALSE, is_parallel = TRUE
      ) %>% 
    get_nets()
  fc
}