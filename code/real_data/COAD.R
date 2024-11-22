dir.base <- "."
dir.data <- file.path(dir.base, "data_processed/")

# Data folder
cohort <- "COAD"
dir.data.raw <- file.path(dir.data, cohort, "raw/")
dir.data.processed <- file.path(dir.data, cohort, "processed/")

# Results folder
dir.results <- file.path(dir.base, "data_results/", cohort)

# Load require library
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(survival)
library(survminer)
library(ggvenn)
library(multiRF)
library(IntNMF)

load(file.path(dir.data.processed, "COAD_three_omics_top_2000.rda"))

# function
source("paper1/code/utility.R")

df <- NULL
mod_list <- list()

for (i in c(1:10)) {
  mods <- mrf3_init(
    tcga_coad,
    connect_list = list( c("gene", "mirna"), c("methy", "mirna"), c("mirna", "gene")),
    scale = F,
    ntree = 300
  )
  
  mods_ls <- plyr::llply(
    c("filter", "mixture", "test"),
    .fun = function(j) {
      c1 <- "normal"
      mods_vs <-  mrf3_vs(mods, dat.list = tcga_coad,
                          method = j,
                          c1 = c1,
                          scale = F)
      mods_vs
    }
  )
  
  names(mods_ls) <- c("filter", "mixture", "test")
  
  num_var <- plyr::ldply(
    mods_ls,
    .fun = function(m) {
      num_var <- m$weights %>% 
        purrr::map(., ~length(.[. > 0])) %>% 
        unlist()
      num_var <- data.frame(num_var_selected = num_var,
                            dataset = names(num_var))
      num_var
    }, .id = "model"
  )
  
  num_var$rep <- i
  
  # NMF clustering
  
  coad_nmf <- plyr::llply(
    mods_ls, 
    .fun = function(mod) {
      #k <- nmf.opt.k(df, k.range = 2:9, n.runs = 10, n.fold = 3, make.plot = FALSE)
      fit <- nmf.mnnals(dat=mod$dat.list %>% purrr::map(~as.matrix(.)), k=2, ini.nndsvd=TRUE, seed=TRUE)
      fit
    }
  )
  names(coad_nmf) <- c("filter", "mixture", "test")
  
  # Logrank test
  
  logrank <- plyr::ldply(
    coad_nmf,
    .fun = function(fit) {
      num_cl <- length(unique(fit$clusters))
      df <- data.frame(scores = factor(fit$clusters), time = tcga_coad_clinical$OS.time, death = tcga_coad_clinical$OS)
      fo <- as.formula(paste0("Surv(time, death) ~ scores"))
      fit <- surv_fit(fo, data = df)
      data.frame(cluster = num_cl, pval = survdiff(fo, data = df)$p)
    }, .id = "model"
  )
  
  num_var <- left_join(num_var, logrank)
  
  df <- rbind(df, num_var)
    
  mod_list <- c(mod_list, list(list(mrf = mods_ls, nmf = coad_nmf)))
}

save(
  mod_list,
  df,
  file = file.path(dir.results, "COAD_10_rep_results.rda")
)

