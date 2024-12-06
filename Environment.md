# Main R Environment

R version 4.0.5 (2021-03-31)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS/LAPACK: /opt/conda/lib/libopenblasp-r0.3.15.so

locale:
 [1] LC_CTYPE=en_US.UTF-8      
 [2] LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8      
 [8] LC_NAME=C                 
 [9] LC_ADDRESS=C              
[10] LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8
[12] LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices
[6] utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.2.0                                    
 [2] ggnewscale_0.5.0                                   
 [3] lmerTest_3.1-3                                     
 [4] lme4_1.1-29                                        
 [5] Matrix_1.4-1                                       
 [6] stringr_1.5.1                                      
 [7] writexl_1.5.0                                      
 [8] readxl_1.4.3                                       
 [9] arrow_13.0.0.1                                     
[10] feather_0.3.5                                      
[11] tidyr_1.3.1                                        
[12] dplyr_1.1.4                                        
[13] IlluminaHumanMethylationEPICv2anno.20a1.hg38_0.99.0
[14] IlluminaHumanMethylationEPICv2manifest_0.99.1      
[15] wateRmelon_1.34.0                                  
[16] illuminaio_0.32.0                                  
[17] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0 
[18] ROC_1.66.0                                         
[19] lumi_2.42.0                                        
[20] methylumi_2.36.0                                   
[21] FDb.InfiniumMethylation.hg19_2.2.0                 
[22] org.Hs.eg.db_3.12.0                                
[23] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2            
[24] GenomicFeatures_1.42.3                             
[25] AnnotationDbi_1.52.0                               
[26] ggplot2_3.5.1                                      
[27] reshape2_1.4.4                                     
[28] scales_1.3.0                                       
[29] limma_3.46.0                                       
[30] ENmix_1.26.10                                      
[31] doParallel_1.0.17                                  
[32] minfi_1.36.0                                       
[33] bumphunter_1.32.0                                  
[34] locfit_1.5-9.4                                     
[35] iterators_1.0.14                                   
[36] foreach_1.5.2                                      
[37] Biostrings_2.58.0                                  
[38] XVector_0.30.0                                     
[39] SummarizedExperiment_1.20.0                        
[40] Biobase_2.50.0                                     
[41] MatrixGenerics_1.2.1                               
[42] matrixStats_1.3.0                                  
[43] GenomicRanges_1.42.0                               
[44] GenomeInfoDb_1.26.7                                
[45] IRanges_2.24.1                                     
[46] S4Vectors_0.28.1                                   
[47] BiocGenerics_0.36.1                                

loaded via a namespace (and not attached):
  [1] utf8_1.2.4                   
  [2] tidyselect_1.2.1             
  [3] RSQLite_2.2.9                
  [4] grid_4.0.5                   
  [5] BiocParallel_1.24.1          
  [6] lpSolve_5.6.20               
  [7] munsell_0.5.1                
  [8] codetools_0.2-20             
  [9] preprocessCore_1.61.0        
 [10] nleqslv_3.3.5                
 [11] withr_3.0.1                  
 [12] colorspace_2.1-1             
 [13] knitr_1.48                   
 [14] GenomeInfoDbData_1.2.4       
 [15] bit64_4.0.5                  
 [16] rhdf5_2.34.0                 
 [17] vctrs_0.6.5                  
 [18] generics_0.1.3               
 [19] xfun_0.47                    
 [20] BiocFileCache_1.14.0         
 [21] R6_2.5.1                     
 [22] bitops_1.0-8                 
 [23] rhdf5filters_1.2.1           
 [24] cachem_1.1.0                 
 [25] reshape_0.8.9                
 [26] DelayedArray_0.16.3          
 [27] assertthat_0.2.1             
 [28] promises_1.3.0               
 [29] gtable_0.3.5                 
 [30] affy_1.68.0                  
 [31] rlang_1.1.4                  
 [32] genefilter_1.72.1            
 [33] splines_4.0.5                
 [34] rtracklayer_1.50.0           
 [35] impute_1.64.0                
 [36] GEOquery_2.58.0              
 [37] BiocManager_1.30.24          
 [38] yaml_2.3.10                  
 [39] httpuv_1.6.6                 
 [40] tools_4.0.5                  
 [41] nor1mix_1.3-3                
 [42] affyio_1.60.0                
 [43] gplots_3.1.3.1               
 [44] RColorBrewer_1.1-3           
 [45] siggenes_1.64.0              
 [46] dynamicTreeCut_1.63-1        
 [47] Rcpp_1.0.9                   
 [48] plyr_1.8.9                   
 [49] sparseMatrixStats_1.2.1      
 [50] progress_1.2.3               
 [51] zlibbioc_1.36.0              
 [52] purrr_1.0.2                  
 [53] RCurl_1.98-1.16              
 [54] prettyunits_1.2.0            
 [55] openssl_2.1.1                
 [56] cluster_2.1.6                
 [57] magrittr_2.0.3               
 [58] data.table_1.15.4            
 [59] hms_1.1.3                    
 [60] mime_0.12                    
 [61] xtable_1.8-4                 
 [62] RPMM_1.25                    
 [63] XML_3.99-0.17                
 [64] mclust_6.1.1                 
 [65] irr_0.84.1                   
 [66] compiler_4.0.5               
 [67] biomaRt_2.57.1               
 [68] tibble_3.2.1                 
 [69] KernSmooth_2.23-24           
 [70] crayon_1.5.3                 
 [71] minqa_1.2.8                  
 [72] htmltools_0.5.8.1            
 [73] mgcv_1.9-1                   
 [74] later_1.3.2                  
 [75] tzdb_0.4.0                   
 [76] geneplotter_1.68.0           
 [77] DBI_1.2.3                    
 [78] ExperimentHub_1.16.1         
 [79] dbplyr_2.5.0                 
 [80] MASS_7.3-60.0.1              
 [81] rappdirs_0.3.3               
 [82] boot_1.3-30                  
 [83] readr_2.1.5                  
 [84] cli_3.6.3                    
 [85] quadprog_1.5-8               
 [86] pkgconfig_2.0.3              
 [87] GenomicAlignments_1.26.0     
 [88] numDeriv_2016.8-1.1          
 [89] xml2_1.3.6                   
 [90] annotate_1.68.0              
 [91] rngtools_1.5.2               
 [92] multtest_2.46.0              
 [93] beanplot_1.3.1               
 [94] doRNG_1.8.6                  
 [95] scrime_1.3.5                 
 [96] digest_0.6.37                
 [97] base64_2.0.1                 
 [98] cellranger_1.1.0             
 [99] DelayedMatrixStats_1.12.3    
[100] curl_5.0.0                   
[101] shiny_1.9.1                  
[102] Rsamtools_2.6.0              
[103] gtools_3.9.5                 
[104] nloptr_2.1.1                 
[105] lifecycle_1.0.4              
[106] nlme_3.1-166                 
[107] jsonlite_1.8.8               
[108] Rhdf5lib_1.12.1              
[109] askpass_1.2.0                
[110] fansi_1.0.6                  
[111] pillar_1.9.0                 
[112] lattice_0.22-6               
[113] fastmap_1.2.0                
[114] httr_1.4.7                   
[115] survival_3.7-0               
[116] interactiveDisplayBase_1.28.0
[117] glue_1.7.0                   
[118] BiocVersion_3.12.0           
[119] bit_4.0.5                    
[120] stringi_1.7.6                
[121] HDF5Array_1.18.1             
[122] blob_1.2.4                   
[123] AnnotationHub_2.22.1         
[124] caTools_1.18.2               
[125] memoise_2.0.1

# R Environment for the MuMIn package

R version 4.2.2 (2022-10-31)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS/LAPACK: /opt/conda/envs/EpiClock_R/lib/libopenblasp-r0.3.21.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] MuMIn_1.48.4

loaded via a namespace (and not attached):
[1] compiler_4.2.2  Matrix_1.5-3    nlme_3.1-161    grid_4.2.2     
[5] stats4_4.2.2    lattice_0.20-45

# Python environment

Python Session Information
Python Version: 3.10.12 (main, Nov  6 2024, 20:22:13) [GCC 11.4.0]
Platform: Linux-5.15.0-76-generic-x86_64-with-glibc2.35
OS: Linux
Release: 5.15.0-76-generic
Machine: x86_64
Processor: 

Installed Packages:
appdirs==1.4.4
asttokens==2.4.1
biolearn==0.3.1
certifi==2023.11.17
charset-normalizer==3.3.2
comm==0.2.1
debugpy==1.8.0
decorator==5.1.1
exceptiongroup==1.2.0
executing==2.0.1
feather-format==0.4.1
idna==3.6
ipykernel==6.29.0
ipython==8.21.0
jedi==0.19.1
joblib==1.3.2
jupyter_client==8.6.0
jupyter_core==5.7.1
matplotlib-inline==0.1.6
mysql-connector==2.2.9
nest-asyncio==1.6.0
numpy==1.26.3
packaging==23.2
pandas==2.2.0
parso==0.8.3
pexpect==4.9.0
pip==23.3.2
platformdirs==4.2.0
prompt-toolkit==3.0.43
psutil==5.9.8
ptyprocess==0.7.0
pure-eval==0.2.2
pyarrow==15.0.0
Pygments==2.17.2
python-dateutil==2.8.2
pytz==2023.4
PyYAML==6.0.1
pyzmq==25.1.2
requests==2.31.0
scikit-learn==1.4.0
scipy==1.12.0
six==1.16.0
stack-data==0.6.3
threadpoolctl==3.2.0
tornado==6.4
traitlets==5.14.1
tzdata==2023.4
urllib3==2.2.0
wcwidth==0.2.13

