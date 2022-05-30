cluster_seurat <-
  function(seurat_object, reduction = "harmony", dims, resolution) {
    seurat_object %<>%
      RunUMAP(reduction = reduction, dims = dims) %>%
      FindNeighbors(reduction = reduction, dims = dims) %>%
      FindClusters(resolution = resolution)
  }
