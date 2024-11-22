dir.base <- "~/TBL Dropbox/Wei Zhang/mrf/"
dir.data <- file.path(dir.base, "data_processed/PAN/")
dir.data.processed <- file.path(dir.data, "processed")
dir.results <- file.path(dir.base, "data_results/PAN")
# --------------------------------------------------------------------------------------------------------------------------
library(multiRF)
library(SummarizedExperiment)
library(tidyverse)
library(IntNMF)
# --------------------------------------------------------------------------------------------------------------------------
# Add annotation
meta <- data.table::fread(
  file.path(dir.data, "raw/TCGA_ATAC_peak.all.probeMap")
) 
meta2 <- data.table::fread(
  file.path(dir.data, "raw/TCGA-ATAC_PanCancer_PeakSet.txt")
) %>% dplyr::mutate(id = name, 
                    chrom = seqnames, 
                    .keep = "unused")
atac_anno <- left_join(
  meta,
  meta2
)
# --------------------------------------------------------------------------------------------------------------------------
# Aux function
get_most_variable_gene <- function(exp, zero_thres = 1, top = 2000, transpose = T){
  
  exp_zero <- rowSums(exp == 0)/ncol(exp)
  exp_filtered <- exp[exp_zero < zero_thres,]
  
  if(nrow(exp_filtered) > top){
    iqr <- order(matrixStats::rowVars(as.matrix(exp_filtered)), decreasing = T)[1:top]
    exp_filtered <- exp_filtered[iqr,]
  }
  
  if(transpose){
    exp_filtered <- exp_filtered %>% t() %>% as.data.frame()
  }
  
  return(exp_filtered)
}
# --------------------------------------------------------------------------------------------------------------------------
# Load data
load(
  file.path(dir.data.processed, "PAN_tpm_atac_clinical.rda")
)
# distal
atac_anno_dis <- atac_anno %>% 
  filter(!annotation %in% "Promoter")
pan_atac_dis <- PAN_atac[atac_anno_dis$id,]
# promoter
pan_atac_pro <- PAN_atac[!rownames(PAN_atac) %in% atac_anno_dis$id,]
PAN_tpm <- log2(PAN_tpm + 1)
# --------------------------------------------------------------------------------------------------------------------------
message("50000 + 5000")
## 50000 + 5000
# pan_tpm_filtered <- get_most_variable_gene(PAN_tpm, top = 5000)
pan_atac_filtered <- get_most_variable_gene(PAN_atac, top = 50000)
pan_filtered <- list(
  atac = pan_atac_filtered,
  rna = pan_tpm_filtered
)
mods <- mrf3_init(pan_filtered, connect_list = list(c("rna", "atac")), ntree = 300, scale = T, yprob = .5)
save(
  mods,
  file = file.path(dir.results, "PAN_all_atac_50000_rna_5000_mod.rda")
)
gc()
# Block
pan_atac_filtered <- get_most_variable_gene(PAN_atac, top = 50000)
pan_atac_list <- plyr::llply(
  seq(1,ncol(pan_atac_filtered), by = 5000),
  .fun = function(i) {
    pan_atac_filtered[,i:(i + 4999)]
  }
)
names(pan_atac_list) <- paste0("block", 1:10, "atac")
pan_filtered <- c(
  pan_atac_list,
  list(rna = pan_tpm_filtered)
)
mods <- mrf3_init(pan_filtered, connect_list = lapply(names(pan_atac_list), function(i) c("rna",i)), ntree = 300, scale = T, yprob = .5)
save(
  mods,
  file = file.path(dir.results, "PAN_all_block_atac_50000_rna_5000_mod.rda")
)
gc()
