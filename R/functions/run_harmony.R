run_harmony <- function(seurat_object) {
  seurat_object %<>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = 100, verbose = FALSE) %>%
    RunHarmony("orig.ident", max.iter.harmony = 50)
}
