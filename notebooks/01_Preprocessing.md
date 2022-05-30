---
title: "01 Preprocessing"
author: "Jussi Kupari"
date: "Last edited: 2022-05-30"
output: html_notebook
---


```r
# Load libraries and custom functions
suppressPackageStartupMessages({
  library(here)
  library(fs)
  library(tidyverse)
  library(magrittr)
  library(Seurat)
})

R.utils::sourceDirectory(here("R", "functions"), modifiedOnly=FALSE)
```


```r
# Create single Seurat object from individual datasets
sample_names <-
  dir_ls(here("data", "raw")) %>%
  map_chr(str_extract, "P2.*")

drg_ra <-
  map(dir_ls(here("data", "raw")), Read10X) %>%
  map(CreateSeuratObject) %>%
  set_names(sample_names)

drg_ra %<>%
  map2(sample_names, ~
    { .x@meta.data %<>%
      mutate(orig.ident = .y)
      return(.x) })

drg_ra %<>% reduce(merge)
```

```
## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.

## Warning in CheckDuplicateCellNames(object.list = objects): Some cell names are
## duplicated across objects provided. Renaming to enforce unique cell names.
```


```r
# Add metadata: percent.mt, sample_info, timepoint_info (days since induction of RA)

drg_ra[["percent.mt"]] <- PercentageFeatureSet(drg_ra, pattern = "^mt-")

metadata <- readRDS(here("data", "metadata.rds"))

drg_ra@meta.data %<>%
  rownames_to_column("cell") %>%
  left_join(metadata, by = c("orig.ident" = "NGI_ID")) %>%
  column_to_rownames("cell")

drg_ra@meta.data %<>%
  rownames_to_column("cell") %>%
  mutate(
    treatment = ifelse(str_detect(Sample_information, "Control"), "Control", "RA"),
    timepoint_days = case_when(
      str_detect(Sample_information, "6hr_") ~ 0.25,
      str_detect(Sample_information, "12hr_") ~ 0.5,
      str_detect(Sample_information, "24hr_") ~ 1,
      str_detect(Sample_information, "48hr_") ~ 2,
      str_detect(Sample_information, "D3_") ~ 3,
      str_detect(Sample_information, "D12_") ~ 12,
      str_detect(Sample_information, "D33_") ~ 33,
      str_detect(Sample_information, "D63_") ~ 63,
      TRUE ~ 0
    )) %>%
  column_to_rownames("cell")

drg_ra %<>% 
  set_idents("timepoint_days") %>% 
  subset(idents = c("0", "0.25", "0.5", "1", "2", "12", "33", "63"))
```


```r
# Save preprocessed seurat object
saveRDS(drg_ra, file = here("data", "proc", "drg_ra_preprocessed.rds"))
```


```r
# Sessioninfo
sessionInfo()
```

```
## R version 4.1.1 (2021-08-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS:   /sw/apps/R/4.1.1/rackham/lib64/R/lib/libRblas.so
## LAPACK: /sw/apps/R/4.1.1/rackham/lib64/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] SeuratObject_4.0.4 Seurat_4.1.0       magrittr_2.0.2     forcats_0.5.1     
##  [5] stringr_1.4.0      dplyr_1.0.8        purrr_0.3.4        readr_2.1.2       
##  [9] tidyr_1.2.0        tibble_3.1.6       ggplot2_3.3.5      tidyverse_1.3.1   
## [13] fs_1.5.2           here_1.0.1        
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.15            colorspace_2.0-3      deldir_1.0-6         
##   [4] ellipsis_0.3.2        ggridges_0.5.3        rprojroot_2.0.2      
##   [7] spatstat.data_2.1-2   rstudioapi_0.13       leiden_0.3.9         
##  [10] listenv_0.8.0         ggrepel_0.9.1         fansi_1.0.2          
##  [13] lubridate_1.8.0       xml2_1.3.3            codetools_0.2-18     
##  [16] splines_4.1.1         R.methodsS3_1.8.1     knitr_1.37           
##  [19] polyclip_1.10-0       jsonlite_1.8.0        broom_0.7.12         
##  [22] ica_1.0-2             cluster_2.1.2         dbplyr_2.1.1         
##  [25] png_0.1-7             R.oo_1.24.0           uwot_0.1.11          
##  [28] spatstat.sparse_2.1-0 sctransform_0.3.3     shiny_1.7.1          
##  [31] compiler_4.1.1        httr_1.4.2            backports_1.4.1      
##  [34] lazyeval_0.2.2        assertthat_0.2.1      Matrix_1.3-4         
##  [37] fastmap_1.1.0         cli_3.2.0             later_1.3.0          
##  [40] htmltools_0.5.2       tools_4.1.1           igraph_1.2.11        
##  [43] gtable_0.3.0          glue_1.6.2            reshape2_1.4.4       
##  [46] RANN_2.6.1            Rcpp_1.0.8            scattermore_0.8      
##  [49] cellranger_1.1.0      styler_1.7.0.9001     vctrs_0.3.8          
##  [52] nlme_3.1-153          lmtest_0.9-39         spatstat.random_2.1-0
##  [55] xfun_0.30             globals_0.14.0        rvest_1.0.2          
##  [58] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.1      
##  [61] irlba_2.3.5           goftest_1.2-3         future_1.24.0        
##  [64] MASS_7.3-54           zoo_1.8-9             scales_1.1.1         
##  [67] spatstat.core_2.4-0   spatstat.utils_2.3-0  hms_1.1.1            
##  [70] promises_1.2.0.1      parallel_4.1.1        RColorBrewer_1.1-2   
##  [73] gridExtra_2.3         pbapply_1.5-0         reticulate_1.24      
##  [76] rpart_4.1-15          stringi_1.7.6         rlang_1.0.2          
##  [79] pkgconfig_2.0.3       matrixStats_0.61.0    evaluate_0.15        
##  [82] lattice_0.20-45       tensor_1.5            ROCR_1.0-11          
##  [85] htmlwidgets_1.5.4     patchwork_1.1.1       cowplot_1.1.1        
##  [88] tidyselect_1.1.2      parallelly_1.30.0     RcppAnnoy_0.0.19     
##  [91] plyr_1.8.6            R6_2.5.1              generics_0.1.2       
##  [94] DBI_1.1.2             mgcv_1.8-38           pillar_1.7.0         
##  [97] haven_2.4.3           withr_2.5.0           fitdistrplus_1.1-8   
## [100] abind_1.4-5           survival_3.2-13       future.apply_1.8.1   
## [103] modelr_0.1.8          crayon_1.5.0          KernSmooth_2.23-20   
## [106] utf8_1.2.2            spatstat.geom_2.3-2   plotly_4.10.0        
## [109] tzdb_0.2.0            grid_4.1.1            readxl_1.3.1         
## [112] data.table_1.14.2     reprex_2.0.1          digest_0.6.29        
## [115] xtable_1.8-4          R.cache_0.15.0        httpuv_1.6.5         
## [118] R.utils_2.11.0        munsell_0.5.0         viridisLite_0.4.0    
## [121] ezknitr_0.6.1
```


```r
# The End
clear_libraries()
```
