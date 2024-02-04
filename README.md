## HCC_SPH

This is the code for HCC spatially punctuated heterogeneity distribution.

#### *Directory structure*

The code for HCC spatially punctuated heterogeneity distribution (HCC_SPH) is organized into different directories and scripts.

The directory structure is as follows:

* Data :  This directory contains example data that is required for regenerating plots.
* DNA : This directory contains the analysis scripts that are based on DNA-related data. 
* RNA : This directory contains the analysis scripts that are based on RNA-related data. 
* Figure :  This directory contains the scripts that are related to reproducing the main figures.

#### *Dependencies*

Session info:

```{R}
> sessionInfo()
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22621)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggfortify_0.4.16            apTreeshape_1.5-0.1         GSVA_1.50.0                
 [4] fgsea_1.28.0                data.table_1.14.10          DESeq2_1.42.0              
 [7] SummarizedExperiment_1.32.0 MatrixGenerics_1.14.0       matrixStats_1.2.0          
[10] GenomicRanges_1.54.1        GenomeInfoDb_1.38.2         colorRamps_2.3.1           
[13] dendextend_1.17.1           scales_1.3.0                ggsci_3.0.0                
[16] msigdbr_7.5.1               org.Hs.eg.db_3.18.0         AnnotationDbi_1.64.1       
[19] IRanges_2.36.0              S4Vectors_0.40.2            Biobase_2.62.0             
[22] BiocGenerics_0.48.1         clusterProfiler_4.10.0      circlize_0.4.15            
[25] ConsensusClusterPlus_1.66.0 survminer_0.4.9             survival_3.5-7             
[28] RColorBrewer_1.1-3          ComplexHeatmap_2.18.0       ape_5.7-1                  
[31] ggpubr_0.6.0                lubridate_1.9.3             forcats_1.0.0              
[34] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.2                
[37] readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1               
[40] ggplot2_3.4.4               tidyverse_2.0.0            

loaded via a namespace (and not attached):
  [1] cubature_2.1.0              splines_4.3.2               bitops_1.0-7               
  [4] ggplotify_0.1.2             polyclip_1.10-6             graph_1.80.0               
  [7] XML_3.99-0.16.1             lifecycle_1.0.4             rstatix_0.7.2              
 [10] doParallel_1.0.17           vroom_1.6.5                 lattice_0.21-9             
 [13] MASS_7.3-60                 backports_1.4.1             magrittr_2.0.3             
 [16] pbapply_1.7-2               cowplot_1.1.2               DBI_1.1.3                  
 [19] abind_1.4-5                 zlibbioc_1.48.0             ggraph_2.1.0               
 [22] RCurl_1.98-1.13             yulab.utils_0.1.1           tweenr_2.0.2               
 [25] GenomeInfoDbData_1.2.11     KMsurv_0.1-5                enrichplot_1.22.0          
 [28] ggrepel_0.9.4               irlba_2.3.5.1               tidytree_0.4.6             
 [31] MatrixModels_0.5-3          annotate_1.80.0             DelayedMatrixStats_1.24.0  
 [34] codetools_0.2-19            DelayedArray_0.28.0         DOSE_3.28.2                
 [37] ggforce_0.4.1               tidyselect_1.2.0            shape_1.4.6                
 [40] aplot_0.2.2                 farver_2.1.1                ScaledMatrix_1.10.0        
 [43] viridis_0.6.4               jsonlite_1.8.8              GetoptLong_1.0.5           
 [46] tidygraph_1.2.3             iterators_1.0.14            foreach_1.5.2              
 [49] ggnewscale_0.4.9            tools_4.3.2                 treeio_1.26.0              
 [52] Rcpp_1.0.11                 glue_1.6.2                  gridExtra_2.3              
 [55] SparseArray_1.2.2           mgcv_1.9-0                  xfun_0.41                  
 [58] qvalue_2.34.0               HDF5Array_1.30.0            withr_2.5.2                
 [61] fastmap_1.1.1               rhdf5filters_1.14.1         fansi_1.0.6                
 [64] SparseM_1.81                rsvd_1.0.5                  digest_0.6.33              
 [67] timechange_0.2.0            R6_2.5.1                    gridGraphics_0.5-1         
 [70] colorspace_2.1-0            GO.db_3.18.0                RSQLite_2.3.4              
 [73] utf8_1.2.4                  generics_0.1.3              graphlayouts_1.0.2         
 [76] httr_1.4.7                  S4Arrays_1.2.0              scatterpie_0.2.1           
 [79] pkgconfig_2.0.3             gtable_0.3.4                blob_1.2.4                 
 [82] SingleCellExperiment_1.24.0 XVector_0.42.0              survMisc_0.5.6             
 [85] shadowtext_0.1.2            carData_3.0-5               GSEABase_1.64.0            
 [88] clue_0.3-65                 png_0.1-8                   ggfun_0.1.3                
 [91] knitr_1.45                  km.ci_0.5-6                 rstudioapi_0.15.0          
 [94] tzdb_0.4.0                  reshape2_1.4.4              rjson_0.2.21               
 [97] coda_0.19-4                 nlme_3.1-163                rhdf5_2.46.1               
[100] cachem_1.0.8                zoo_1.8-12                  GlobalOptions_0.1.2        
[103] parallel_4.3.2              HDO.db_0.99.1               pillar_1.9.0               
[106] vctrs_0.6.5                 BiocSingular_1.18.0         car_3.1-2                  
[109] beachmat_2.18.0             xtable_1.8-4                cluster_2.1.4              
[112] cli_3.6.2                   locfit_1.5-9.8              compiler_4.3.2             
[115] rlang_1.1.2                 crayon_1.5.2                ggsignif_0.6.4             
[118] labeling_0.4.3              plyr_1.8.9                  fs_1.6.3                   
[121] stringi_1.8.3               viridisLite_0.4.2           BiocParallel_1.36.0        
[124] babelgene_22.9              munsell_0.5.0               Biostrings_2.70.1          
[127] lazyeval_0.2.2              quantreg_5.97               GOSemSim_2.28.0            
[130] Matrix_1.6-4                hms_1.1.3                   patchwork_1.1.3            
[133] sparseMatrixStats_1.14.0    bit64_4.0.5                 Rhdf5lib_1.24.1            
[136] KEGGREST_1.42.0             igraph_1.6.0                broom_1.0.5                
[139] memoise_2.0.1               ggtree_3.10.0               fastmatch_1.1-4            
[142] bit_4.0.5                   gson_0.1.0 
```

#### *Data source*

The raw exome and RNA sequencing data have been deposited in the Genome Sequence Archive (GSA) hosted at National Genomics Data Center [HRA002112](), processed data can be found at source data file from the article.
