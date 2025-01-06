# -----------------------------------------------------------------------------------------------------------
# For reproducible research, please install the following R packages 
# and make sure the R and BiocManager versions are correct
# ─ Session info ───────────────────────────────────────────────
# setting  value
# version  R version 4.4.2 (2024-10-31)
# os       Ubuntu 22.04.5 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-11-20
# rstudio  2024.04.2+764.pro1 Chocolate Cosmos (server)
# pandoc   3.1.11 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# The following packages should be load before running the analysis
# -----------------------------------------------------------------------------------------------------------
list.of.packages <- c(
  "caret",
  "clusterProfiler",
  "corpcor",
  "devtools",
  "dplyr",
  "doParallel",
  "genefilter",
  "GEOquery",
  "glmnet",
  "ggplot2",
  "ggpubr",
  "InterSIM",
  "janitor",
  "MASS",
  "matrixStats",
  "mixOmics",
  "msigdbr", 
  "plyr",
  "purrr",
  "pheatmap",
  "PMA",
  "randomForestSRC",
  "readr",
  "RGCCA",
  "rstatix",
  "survival",
  "survminer",
  "TCGAbiolinks",
  "tibble",
  "umap"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  for(new in new.packages){
    if(new %in% available.packages()[,1]){
      install.packages(new)
    } else BiocManager::install(new)
  }
} 

for (pkg in list.of.packages) {
  library(pkg,character.only = TRUE)
}