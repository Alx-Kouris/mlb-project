# Machine Learning in Computational Biology
# Final Project from Alex Kouris

Due to the file sizes and github limitations, some of the files need to be downloaded separately and couldn't be included in the repository.

## Reference data download

The reference dataset must be downloaded from the scMatch github repository:
```bat
https://github.com/asrhou/scMatch/tree/master/refDB/FANTOM5
```
The files to download are:
- 9606_map.csv
- 9606_symbol.csv.zip

After download the files, extract the symbol file and copy over the map and the symbol files into: `data/fantom5`
Another file called `map.csv` already exists in that folder, it is needed too.

## Query data download
For the 68k PBMC dataset that has been used as query dataset, it must be downloaded from the 10X Genomics platform.
```bat
https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz
```
After finishing downloading the file, extract it and copy over the `filtered_matrices_mex` folder into:
```bat
data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices
```
The annotation data for this dataset has been already been included in the above folder.

## Data preprocess
In the notebooks folder, a `data_preprocess.ipynb` notebook exists, that read the reference dataset and creates the gene expressions csvs, that are needed as input for the R scripts.
After running the notebook, under the data directory, a new directory called `preprocessed` should have been created, with two files inside:
- `gene_expressions.csv`
- `gene_expressions_augmented.csv`

## Running scPred
1. To install the scPred tool follow the instructions found here:
   ```bat
   https://powellgenomicslab.github.io/scPred/index.html
   ```
2. Under the `src` folder, there are R scripts that run the different study cases:
   - `scpred_example_run.R`: Runs the example datasets from the scPred introduction, in order to familirize with the tool.
   - `scpred_fantom5.R`: Runs the case where FANTOM5 is used as reference dataset and the 68k PBMC dataset as query dataset.
   - `scpred_augmented_fantom5.R`: Runs the case where the augmented FANTOM5 is used as reference dataset and the 68k PBMC dataset as query dataset.
   - `scpred_pbmc.R`: Runs the case where the 68k PBMC dataset is used as reference dataset, it's split into validate, train, test sets and then it's used to classify the FANTOM5 dataset.
   - `utils.R`: A file with helper functions that are used throughout the other scpred scripts. It's sourced in all the files.
  Each file can be separately, there are no dependecies between them.

## Session info for R
```bat
> devtools::session_info()
─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.1 (2024-06-14)
 os       Ubuntu 22.04.4 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Athens
 date     2024-09-05
 rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
 pandoc   NA

─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 package          * version    date (UTC) lib source
 abind              1.4-5      2016-07-21 [2] CRAN (R 4.4.1)
 adabag             5.0        2023-05-31 [1] CRAN (R 4.4.1)
 base64enc          0.1-3      2015-07-28 [2] CRAN (R 4.4.1)
 beeswarm           0.4.0      2021-06-01 [2] CRAN (R 4.4.1)
 cachem             1.1.0      2024-05-16 [2] CRAN (R 4.4.1)
 caret            * 6.0-94     2023-03-21 [2] CRAN (R 4.4.1)
 class              7.3-22     2023-05-03 [4] CRAN (R 4.3.1)
 cli                3.6.3      2024-06-21 [2] CRAN (R 4.4.1)
 cluster            2.1.6      2023-12-01 [4] CRAN (R 4.3.2)
 codetools          0.2-20     2024-03-31 [4] CRAN (R 4.4.1)
 colorspace         2.1-1      2024-07-26 [2] CRAN (R 4.4.1)
 ConsRank           2.1.4      2024-01-24 [1] CRAN (R 4.4.1)
 cowplot            1.1.3      2024-01-22 [2] CRAN (R 4.4.1)
 crayon             1.5.3      2024-06-20 [2] CRAN (R 4.4.1)
 data.table         1.16.0     2024-08-27 [2] CRAN (R 4.4.1)
 dbscan             1.2-0      2024-06-28 [2] CRAN (R 4.4.1)
 deldir             2.0-4      2024-02-28 [2] CRAN (R 4.4.1)
 devtools           2.4.5      2022-10-11 [2] CRAN (R 4.4.1)
 digest             0.6.37     2024-08-19 [2] CRAN (R 4.4.1)
 doParallel         1.0.17     2022-02-07 [1] CRAN (R 4.4.1)
 dotCall64          1.1-1      2023-11-28 [2] CRAN (R 4.4.1)
 dplyr            * 1.1.4      2023-11-17 [2] CRAN (R 4.4.1)
 ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.4.1)
 fansi              1.0.6      2023-12-08 [2] CRAN (R 4.4.1)
 farver             2.1.2      2024-05-13 [2] CRAN (R 4.4.1)
 fastDummies        1.7.4      2024-08-16 [2] CRAN (R 4.4.1)
 fastmap            1.2.0      2024-05-15 [2] CRAN (R 4.4.1)
 fitdistrplus       1.2-1      2024-07-12 [2] CRAN (R 4.4.1)
 FNN                1.1.4      2024-01-12 [2] CRAN (R 4.4.1)
 foreach            1.5.2      2022-02-02 [2] CRAN (R 4.4.1)
 fs                 1.6.4      2024-04-25 [2] CRAN (R 4.4.1)
 future             1.34.0     2024-07-29 [2] CRAN (R 4.4.1)
 future.apply       1.11.2     2024-03-28 [2] CRAN (R 4.4.1)
 generics           0.1.3      2022-07-05 [2] CRAN (R 4.4.1)
 ggbeeswarm         0.7.2      2023-04-29 [2] CRAN (R 4.4.1)
 ggplot2          * 3.5.1      2024-04-23 [2] CRAN (R 4.4.1)
 ggrepel            0.9.5      2024-01-10 [2] CRAN (R 4.4.1)
 ggridges           0.5.6      2024-01-23 [2] CRAN (R 4.4.1)
 globals            0.16.3     2024-03-08 [2] CRAN (R 4.4.1)
 glue               1.7.0      2024-01-09 [2] CRAN (R 4.4.1)
 goftest            1.2-3      2021-10-07 [2] CRAN (R 4.4.1)
 gower              1.0.1      2022-12-22 [2] CRAN (R 4.4.1)
 gridExtra          2.3        2017-09-09 [2] CRAN (R 4.4.1)
 gtable             0.3.5      2024-04-22 [2] CRAN (R 4.4.1)
 gtools             3.9.5      2023-11-20 [2] CRAN (R 4.4.1)
 hardhat            1.4.0      2024-06-02 [2] CRAN (R 4.4.1)
 harmony            1.2.0      2024-08-24 [2] Github (immunogenomics/harmony@f054b03)
 htmltools          0.5.8.1    2024-04-04 [2] CRAN (R 4.4.1)
 htmlwidgets        1.6.4      2023-12-06 [2] CRAN (R 4.4.1)
 httpuv             1.6.15     2024-03-26 [2] CRAN (R 4.4.1)
 httr               1.4.7      2023-08-15 [2] CRAN (R 4.4.1)
 ica                1.0-3      2022-07-08 [2] CRAN (R 4.4.1)
 igraph             2.0.3      2024-03-13 [2] CRAN (R 4.4.1)
 ipred              0.9-15     2024-07-18 [2] CRAN (R 4.4.1)
 irlba              2.3.5.1    2022-10-03 [2] CRAN (R 4.4.1)
 iterators          1.0.14     2022-02-05 [2] CRAN (R 4.4.1)
 jsonlite           1.8.8      2023-12-04 [2] CRAN (R 4.4.1)
 kernlab            0.9-33     2024-08-13 [2] CRAN (R 4.4.1)
 KernSmooth         2.23-24    2024-05-17 [4] CRAN (R 4.4.0)
 knitr              1.48       2024-07-07 [2] CRAN (R 4.4.1)
 labeling           0.4.3      2023-08-29 [2] CRAN (R 4.4.1)
 later              1.3.2      2023-12-06 [2] CRAN (R 4.4.1)
 lattice          * 0.22-6     2024-03-20 [4] CRAN (R 4.4.1)
 lava               1.8.0      2024-03-05 [2] CRAN (R 4.4.1)
 lazyeval           0.2.2      2019-03-15 [2] CRAN (R 4.4.1)
 leiden             0.4.3.1    2023-11-17 [2] CRAN (R 4.4.1)
 lifecycle          1.0.4      2023-11-07 [2] CRAN (R 4.4.1)
 listenv            0.9.1      2024-01-29 [2] CRAN (R 4.4.1)
 lmtest             0.9-40     2022-03-21 [2] CRAN (R 4.4.1)
 lubridate          1.9.3      2023-09-27 [2] CRAN (R 4.4.1)
 magrittr         * 2.0.3      2022-03-30 [2] CRAN (R 4.4.1)
 MASS               7.3-61     2024-06-13 [4] CRAN (R 4.4.1)
 Matrix             1.6-5      2024-01-11 [4] CRAN (R 4.3.3)
 matrixStats        1.3.0      2024-04-11 [2] CRAN (R 4.4.1)
 memoise            2.0.1      2021-11-26 [2] CRAN (R 4.4.1)
 mime               0.12       2021-09-28 [2] CRAN (R 4.4.1)
 miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.4.1)
 MLmetrics        * 1.1.3      2024-04-13 [1] CRAN (R 4.4.1)
 ModelMetrics       1.2.2.2    2020-03-17 [2] CRAN (R 4.4.1)
 munsell            0.5.1      2024-04-01 [2] CRAN (R 4.4.1)
 nlme               3.1-166    2024-08-14 [4] CRAN (R 4.4.1)
 nnet               7.3-19     2023-05-03 [4] CRAN (R 4.3.1)
 parallelly         1.38.0     2024-07-27 [2] CRAN (R 4.4.1)
 patchwork          1.2.0      2024-01-08 [2] CRAN (R 4.4.1)
 pbapply            1.7-2      2023-06-27 [2] CRAN (R 4.4.1)
 pillar             1.9.0      2023-03-22 [2] CRAN (R 4.4.1)
 pkgbuild           1.4.4      2024-03-17 [2] CRAN (R 4.4.1)
 pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.4.1)
 pkgload            1.4.0      2024-06-28 [2] CRAN (R 4.4.1)
 plotly             4.10.4     2024-01-13 [2] CRAN (R 4.4.1)
 plyr               1.8.9      2023-10-02 [2] CRAN (R 4.4.1)
 png                0.1-8      2022-11-29 [2] CRAN (R 4.4.1)
 polyclip           1.10-7     2024-07-23 [2] CRAN (R 4.4.1)
 pROC             * 1.18.5     2023-11-01 [2] CRAN (R 4.4.1)
 prodlim            2024.06.25 2024-06-24 [2] CRAN (R 4.4.1)
 profvis            0.3.8      2023-05-02 [2] CRAN (R 4.4.1)
 progressr          0.14.0     2023-08-10 [2] CRAN (R 4.4.1)
 promises           1.3.0      2024-04-05 [2] CRAN (R 4.4.1)
 proxy              0.4-27     2022-06-09 [2] CRAN (R 4.4.1)
 purrr              1.0.2      2023-08-10 [2] CRAN (R 4.4.1)
 R6                 2.5.1      2021-08-19 [2] CRAN (R 4.4.1)
 RANN               2.6.1      2019-01-08 [2] CRAN (R 4.4.1)
 RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.4.1)
 Rcpp               1.0.13     2024-07-17 [2] CRAN (R 4.4.1)
 RcppAnnoy          0.0.22     2024-01-23 [2] CRAN (R 4.4.1)
 RcppHNSW           0.6.0      2024-02-04 [2] CRAN (R 4.4.1)
 recipes            1.1.0      2024-07-04 [2] CRAN (R 4.4.1)
 remotes            2.5.0      2024-03-17 [2] CRAN (R 4.4.1)
 reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.4.1)
 reticulate         1.38.0     2024-06-19 [2] CRAN (R 4.4.1)
 rgl                1.3.1      2024-03-05 [1] CRAN (R 4.4.1)
 RhpcBLASctl        0.23-42    2023-02-11 [2] CRAN (R 4.4.1)
 rlang              1.1.4      2024-06-04 [2] CRAN (R 4.4.1)
 rlist              0.4.6.2    2021-09-03 [1] CRAN (R 4.4.1)
 ROCR               1.0-11     2020-05-02 [2] CRAN (R 4.4.1)
 rpart              4.1.23     2023-12-05 [4] CRAN (R 4.3.2)
 RSpectra           0.16-2     2024-07-18 [2] CRAN (R 4.4.1)
 rstudioapi         0.16.0     2024-03-24 [2] CRAN (R 4.4.1)
 Rtsne              0.17       2023-12-07 [2] CRAN (R 4.4.1)
 scales             1.3.0      2023-11-28 [2] CRAN (R 4.4.1)
 scattermore        1.2        2023-06-12 [2] CRAN (R 4.4.1)
 scPred           * 1.9.2      2024-08-27 [2] Github (powellgenomicslab/scPred@9f407b7)
 sctransform        0.4.1      2023-10-19 [2] CRAN (R 4.4.1)
 sessioninfo        1.2.2      2021-12-06 [2] CRAN (R 4.4.1)
 Seurat           * 5.1.0      2024-05-10 [2] CRAN (R 4.4.1)
 SeuratObject     * 5.0.2      2024-05-08 [2] CRAN (R 4.4.1)
 shiny              1.9.1      2024-08-01 [2] CRAN (R 4.4.1)
 smotefamily      * 1.4.0      2024-03-14 [2] CRAN (R 4.4.1)
 sp               * 2.1-4      2024-04-30 [2] CRAN (R 4.4.1)
 spam               2.10-0     2023-10-23 [2] CRAN (R 4.4.1)
 spatstat.data      3.1-2      2024-06-21 [2] CRAN (R 4.4.1)
 spatstat.explore   3.3-2      2024-08-21 [2] CRAN (R 4.4.1)
 spatstat.geom      3.3-2      2024-07-15 [2] CRAN (R 4.4.1)
 spatstat.random    3.3-1      2024-07-15 [2] CRAN (R 4.4.1)
 spatstat.sparse    3.1-0      2024-06-21 [2] CRAN (R 4.4.1)
 spatstat.univar    3.0-0      2024-06-28 [2] CRAN (R 4.4.1)
 spatstat.utils     3.1-0      2024-08-17 [2] CRAN (R 4.4.1)
 stringi            1.8.4      2024-05-06 [2] CRAN (R 4.4.1)
 stringr            1.5.1      2023-11-14 [2] CRAN (R 4.4.1)
 survival           3.7-0      2024-06-05 [4] CRAN (R 4.4.0)
 tensor             1.5        2012-05-05 [2] CRAN (R 4.4.1)
 tibble             3.2.1      2023-03-20 [2] CRAN (R 4.4.1)
 tidyr              1.3.1      2024-01-24 [2] CRAN (R 4.4.1)
 tidyselect         1.2.1      2024-03-11 [2] CRAN (R 4.4.1)
 timechange         0.3.0      2024-01-18 [2] CRAN (R 4.4.1)
 timeDate           4032.109   2023-12-14 [2] CRAN (R 4.4.1)
 urlchecker         1.0.1      2021-11-30 [2] CRAN (R 4.4.1)
 usethis            3.0.0      2024-07-29 [2] CRAN (R 4.4.1)
 utf8               1.2.4      2023-10-22 [2] CRAN (R 4.4.1)
 uwot               0.2.2      2024-04-21 [2] CRAN (R 4.4.1)
 vctrs              0.6.5      2023-12-01 [2] CRAN (R 4.4.1)
 vipor              0.4.7      2023-12-18 [2] CRAN (R 4.4.1)
 viridisLite        0.4.2      2023-05-02 [2] CRAN (R 4.4.1)
 withr              3.0.1      2024-07-31 [2] CRAN (R 4.4.1)
 xfun               0.47       2024-08-17 [2] CRAN (R 4.4.1)
 xtable             1.8-4      2019-04-21 [2] CRAN (R 4.4.1)
 zoo                1.8-12     2023-04-13 [2] CRAN (R 4.4.1)
```
  
