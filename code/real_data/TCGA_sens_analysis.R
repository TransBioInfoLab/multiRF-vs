# Load related function
dir.base <- "."
script <- list.files( 
  path = file.path(dir.base,"function"),
  pattern = "[.]R$", 
  full.names = T,
  recursive = T
)
for (f in script) source(f)

dir.data <- file.path(dir.base, "$data path$")

# Data folder
cohort <- "BRCA"
dir.data.raw <- file.path(dir.data, cohort, "raw/")
dir.data.processed <- file.path(dir.data, cohort, "processed/")

# Results folder
dir.results <- file.path(dir.base, "results/real_data", cohort)

# Load require library
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(survival)
library(survminer)
library(ggvenn)
library(clusterProfiler)
library(msigdbr)
library(mixOmics)

load(file.path(dir.data.processed, "BRCA_three_omics_top_2000.rda"))


size <- run_different_seed(tcga_brca,
                           status = tcga_brca_clinical$OS,
                           time = tcga_brca_clinical$OS.time,
                           rep = 30,
                           scale = F,
                           connect_list = list( c("gene", "methy"), c("mirna", "methy")))

write_csv(
  size,
  file.path(dir.results, "BRCA_sens_analysis.csv")
)

save(
  size,
  file = file.path(dir.results, "BRCA_sens_analysis.rda")
)

# Data folder
cohort <- "COAD"
dir.data.raw <- file.path(dir.data, cohort, "raw/")
dir.data.processed <- file.path(dir.data, cohort, "processed/")

# Results folder
dir.results <- file.path(dir.base, "results/real_data", cohort)

load(file.path(dir.data.processed, "COAD_three_omics_top_2000.rda"))

size <- run_different_seed(tcga_coad,
                           status = tcga_coad_clinical$OS,
                           time = tcga_coad_clinical$OS.time,
                           rep = 30,
                           connect_list = list( c("mirna", "methy"), c("mirna", "gene")),
                           scale = F)

write_csv(
  size,
  file.path(dir.results, "COAD_sens_analysis.csv")
)

save(
  size,
  file = file.path(dir.results, "COAD_sens_analysis.rda")
)
