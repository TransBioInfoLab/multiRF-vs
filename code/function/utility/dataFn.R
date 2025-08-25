
run_different_seed <- function(dat, time, status, base_seed = 123, rep = 10, scale = F, ...) {
  
  # doParallel::registerDoParallel(10)
  plyr::ldply(
    1:rep,
    .fun = function(i) {
      mods <- mrf3_init(
        dat,
        scale = scale,
        ntree = 300,
        seed = base_seed + i,
        ...
      )
      
      results <- plyr::ldply(
        c("filter", "mixture", "test"),
        .fun = function(j) {
          mods_vs <-  mrf3_vs(mods, dat.list = dat,
                              method = j,
                              scale = scale,
                              re_fit = F, se = 1.5)
          eva <- eval_prognosis(mods_vs$dat.list, time = time, status = status)
          Size <- mods_vs$dat.list %>% purrr::map(., ~ncol(.)) %>% unlist()
          df <- cbind(Size %>% t() %>% data.frame(), eva)
          df$var_selected <- list(mods_vs$dat.list %>% purrr::map(., ~colnames(.)))
          df
        }
      )

      results$method <- c("filter", "mixture", "test")
      results$seed <- base_seed + i
      results
    }, .parallel = F
  )
  
 
}

# -----------------------------
# Prognostic evaluation on selected features
# -----------------------------
eval_prognosis <- function(data_list, time, status) {
  # df <- data_list %>% purrr::map(
  #   ., ~{
  #     if(scale) . <- scale(.)
  #     if(!all(. > 0)) {
  #       m <- abs(min(.))
  #       . <- pmax(. + m, 0)
  #       as.matrix(./max(.))
  #     }
  #   }
  # )
  # NMF rank=2
  nmf_mod  <- IntNMF::nmf.mnnals(data_list %>% purrr::map(., ~as.matrix(.)), 2)
  clusters <- nmf_mod$clusters
  
  surv_obj <- Surv(time, status)
  cox_mod  <- coxph(surv_obj ~ factor(clusters))
  c_index  <- summary(cox_mod)$concordance[1]
  lr       <- survdiff(surv_obj ~ factor(clusters))
  p_lr     <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  df <- data.frame(
    C_index   = c_index,
    logrank_p = p_lr,
    cluster = I(list(clusters))
  )
}


jaccard <- function(a, b) {
  if (length(a)==0 && length(b)==0) return(1)
  intersect_len <- length(intersect(a, b))
  union_len     <- length(union(a, b))
  intersect_len / union_len
}

overlap <- function(a, b) {
  if (length(a)==0 && length(b)==0) return(1)
  intersect_len <- length(intersect(a, b))
  min_len     <- min(length(a), length(b))
  intersect_len / min_len
}

