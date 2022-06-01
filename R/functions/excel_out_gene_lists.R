excel_out_gene_lists <- function(gene_list, out_dir, prefix) {
  
  if (!(prefix %in% c("up", "down", "regulated"))) {
    stop()
  }
  
  gene_list <- compact(gene_list)
  
  for (i in seq_along(gene_list)) {
    name <- names(gene_list[i])
    element <- gene_list[[i]]
    
    openxlsx::write.xlsx(
      element, file = here::here(out_dir, str_glue("{prefix}_{name}.xlsx"))
    )
  }
}
