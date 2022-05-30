plot_modules_ifit1 <- function(celltype) {
  neurons %>% 
    set_idents("usoskin_id") %>% 
    subset(idents = celltype) %>% 
    set_idents("timepoint_days") %>%
    ScaleData(features = rownames(neurons)) %>%
    AverageExpression(return.seurat = TRUE) %>% 
    set_cluster_order(c(0, 0.25, 0.5, 1, 2, 12, 33, 63)) %>% 
    DoHeatmap(features = mods$Ifit1, draw.lines = FALSE, size = 3, raster = FALSE) + 
    NoLegend() + ggtitle(paste("Ifit1_mod", celltype, sep = " "))
}