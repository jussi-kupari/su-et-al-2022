run_fcoex <-
  function(celltype,
           seurat_obj,
           idents_from = "scpred_no_rejection",
           idents_to = "timepoint_days",
           n_genes = 200,
           discretize_method = "median") {
    small_data <-
      seurat_obj %>%
      set_idents(idents_from) %>%
      subset(idents = celltype) %>%
      set_idents(idents_to)

    exprs <- data.frame(GetAssayData(small_data))
    target <- Idents(small_data)

    fc <- new_fcoex(data.frame(exprs), target)
    fc <- discretize(fc, method = discretize_method)
    fc <-
      find_cbf_modules(
        fc,
        n_genes_selected_in_first_step = n_genes,
        verbose = FALSE,
        is_parallel = TRUE
      )
  }
