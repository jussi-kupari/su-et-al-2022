get_markers2 <- function(seurat_obj, celltype, time) {
  subset <-
    seurat_obj %>%
    subset(idents = celltype) %>%
    set_idents("timepoint_days") %>%
    subset(idents = c(0, time))

  markers <-
    subset %>%
    set_idents("treatment") %>%
    FindMarkers(ident.1 = "RA", ident.2 = "Control", assay = "RNA") %>%
    mutate(
      score = abs(log2((pct.1 + 10^-6) / (pct.2 + 10^-6)) * avg_log2FC)
    ) %>%
    filter(p_val_adj <= 0.05) %>%
    mutate(celltype = celltype, timepoint = time) %>%
    rownames_to_column("gene")

  markers
}
