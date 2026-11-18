# Tested Software Environment

Recorded for the local audit on 18 November 2026. This is not a reconstruction
of every detail of the historical Bunya environment.

Hardware: Apple M2 Pro, Mac14,9, 32 GiB RAM (recorded with the completed mock fit).
The sandboxed inventory could not repeat the hardware query; that failure is
retained in output/check/environment.rds rather than substituted with invented metadata.
Operating system: ProductName:		macOS; ProductVersion:		27.2.1; BuildVersion:		26B131
R: R version 4.5.2 (2025-10-31)
CmdStan: 2.40.0

## Required R Packages

| Package | Tested version |
|---|---|
| tidyverse | 2.0.0 |
| survival | 3.8-3 |
| nlme | 3.1-168 |
| MASS | 7.3-65 |
| cmdstanr | 0.9.0 |
| posterior | 1.6.1 |
| loo | 2.9.0 |
| patchwork | 1.3.2 |
| JMbayes2 | 0.6-0 |

The inventory also loads the components used through tidyverse and JMbayes2:

| Package | Tested version |
|---|---|
| dplyr | 1.2.0 |
| ggplot2 | 4.0.2 |
| purrr | 1.2.1 |
| readr | 2.2.0 |
| stringr | 1.6.0 |
| tibble | 3.3.1 |
| tidyr | 1.3.2 |
| forcats | 1.0.1 |
| lubridate | 1.9.5 |
| coda | 0.19.4.1 |

## Compiler And Included Libraries

```text
Apple clang version 21.0.0 (clang-2100.4.1.1)
Target: arm64-apple-darwin27.2.0
Thread model: posix
InstalledDir: /Library/Developer/CommandLineTools/usr/bin
Apple clang version 21.0.0 (clang-2100.4.1.1)
Target: arm64-apple-darwin27.2.0
Thread model: posix
InstalledDir: /Library/Developer/CommandLineTools/usr/bin
Apple clang version 21.0.0 (clang-2100.4.1.1)
Target: arm64-apple-darwin27.2.0
Thread model: posix
InstalledDir: /Library/Developer/CommandLineTools/usr/bin
GNU Make 3.81
Copyright (C) 2006  Free Software Foundation, Inc.
This is free software; see the source for copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

This program built for i386-apple-darwin11.3.0
boost_1.87.0
eigen_5.0.1
sundials_6.1.1
tbb
tbb_2020.3
upgrade-boost.sh
upgrade-sundials.sh
upgrade-tbb.sh
```

Per-fit build.rds, compile.log and linked-libraries.txt record the compiler flags,
source/executable checksums and linked libraries actually used. The no-PCH option
is an explicit build choice. No packages or CmdStan installations are rebuilt
automatically by the audit entry script.

## Loaded-Package Session Information

The following was captured after the inventory loaded all required package
namespaces. Package provenance, including available remote commit hashes, is
also recorded in output/check/environment.rds.

```text
R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS 27.2.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: Australia/Brisbane
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] tensorA_0.36.2.1     generics_0.1.4       tidyr_1.3.2         
 [4] stringi_1.8.7        lattice_0.22-7       survC1_1.0-3        
 [7] hms_1.1.4            magrittr_2.0.4       evaluate_1.0.5      
[10] grid_4.5.2           timechange_0.4.0     RColorBrewer_1.1-3  
[13] Matrix_1.7-4         GLMMadaptive_0.9-7   tidyverse_2.0.0     
[16] backports_1.5.0      survival_3.8-3       gridExtra_2.3       
[19] ps_1.9.1             purrr_1.2.1          scales_1.4.0        
[22] abind_1.4-8          cli_3.6.5            rlang_1.1.7         
[25] cmdstanr_0.9.0       parallelly_1.46.1    splines_4.5.2       
[28] otel_0.2.0           parallel_4.5.2       tools_4.5.2         
[31] tzdb_0.5.0           coda_0.19-4.1        checkmate_2.3.4     
[34] dplyr_1.2.0          ggplot2_4.0.2        forcats_1.0.1       
[37] JMbayes2_0.6-0       vctrs_0.7.1          posterior_1.6.1     
[40] R6_2.6.1             matrixStats_1.5.0    lifecycle_1.0.5     
[43] lubridate_1.9.5      stringr_1.6.0        MASS_7.3-65         
[46] pkgconfig_2.0.3      pillar_1.11.1        gtable_0.3.6        
[49] Rcpp_1.1.1           loo_2.9.0            glue_1.8.0          
[52] xfun_0.56            tibble_3.3.1         tidyselect_1.2.1    
[55] knitr_1.51           farver_2.1.2         patchwork_1.3.2     
[58] nlme_3.1-168         readr_2.2.0          compiler_4.5.2      
[61] S7_0.2.1             distributional_0.6.0
```

## Historical Production Environment

The manuscript simulations used UQ Bunya EPYC3/EPYC4 nodes and CmdStan 2.37.0.
Complete historical OS/compiler/package versions were not retained and must not
be inferred from the local record above. Exact cross-platform MCMC equality
is not claimed. The source, initialization and data-generation equivalence
checks are separate from comparisons of posterior Monte Carlo results.
