library(xgboost)
library(randomForestSRC)
library(mixOmics)
library(PMA)
library(RGCCA)
library(PRROC)
library(doParallel)
library(tidyverse)
library(plyr)
library(gbm3)
library(pROC)
#### Functions for variable selection evaluation
get_all_imp <- function(dat, keep.list, ...) {
  
  ## MRF
  imd <- sim.fn.mrf3.m(dat = dat, ...)

  spls_imp <- plyr::llply(
    c("PMDCCA", "SPLS", "RGCCA"),
    .fun = function(method) {
      sim.fn.other(dat, cca.model = method, keep.list = keep.list)
    }
  )
  names(spls_imp) <- c("PMDCCA", "SPLS", "RGCCA")
  
  c(imd, spls_imp)
}

evaluate_varsel <- function(imp_vec, true_idx) {
  imp_vec <- abs(imp_vec)
  # 1) Precision‐Recall AUC

  pr <- pr.curve(
    scores.class0 = imp_vec[ true_idx == 1],  # positives
    scores.class1 = imp_vec[ true_idx == 0],  # negatives
    curve = FALSE
  )
  prauc_val <- pr$auc.integral
  
  rocobj  <- roc(true_idx, imp_vec, quiet = TRUE)
  auc_val <- auc(rocobj)
  
  # 2) Recall and precision
  k <- sum(true_idx == 1)
  selected <- which(imp_vec != 0)
  tp   <- sum(selected %in% which(true_idx == 1)) 
  prec <- if (length(selected)>0) tp/length(selected) else NA
  rec  <- tp/length(true_idx[true_idx == 1])
  f1        <- if ((prec + rec) == 0) 0 else 
    2 * prec * rec / (prec + rec)
  
  c(PRAUC       = as.numeric(prauc_val),
    AUC = auc_val,
    Recall    = rec,
    Precision = prec,
    F1 = f1,
    Size = length(selected))
    
}

#### Functions for importance evaluation
get_importances <- function(dat, type = "binary", connect_list = list(c("X2", "X1"))) {
  
  if(type == "binary") {
    dat$X$X2 <- lapply(dat$X$X2, as.factor) %>% as.data.frame()
    objective <- "binary:logistic"
    eval_metric <- "auc"
    distribution <- "bernoulli"
  } else {
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
    distribution <- "gaussian"
  }
  X <- as.matrix(dat$X$X1); Y <- dat$X$X2
  p <- ncol(X); m <- ncol(Y)
  
  # — MRF-IMD and MD —
  mrf_mod  <- mrf3_init(dat$X, connect_list = connect_list, yprob = 1, calc = "X", scale = F)
  wl <- mrf_mod$weights_ls[[1]] %>% purrr::map("X")
  wl_m <- plyr::llply(seq(0.1,0.9,0.1),
                      .fun = function(a) {
                        sapply(wl, function(w) {
                          new_w <- w^a
                        }) %>% rowMeans()
                        
                      })
  names(wl_m) <- paste0("IMD_a0", 1:9)
  mrf_imd <- mrf_mod$weights$X1
  vs      <- var.select(mrf_mod$mod[[1]], conservative = "medium")
  mrf_md <- 1/vs$md.obj$order[,1]
  
  # — Gradient Boosting (xgboost) —
  gb_imp_mat <- matrix(0, p, m, dimnames = list(colnames(X), colnames(Y)))
  params <- list(objective = objective,
                 eval_metric = eval_metric)
  
  for (j in seq_len(m)) {
    if(type == "binary") {
      y <- as.numeric(Y[, j]) - 1
    } else y <- Y[,j]
    
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y)
    model  <- xgboost::xgb.train(params, dtrain, nrounds = 10, verbose = 0)
    imp_tbl <- xgboost::xgb.importance(model = model)[, c("Feature","Gain")]
    gb_imp_mat[imp_tbl$Feature, j] <- imp_tbl$Gain
  }
  gb_imp <- rowMeans(gb_imp_mat)
  
  gbm2_imp_mat <- matrix(0, p, m, dimnames = list(colnames(X), colnames(Y)))
  for (j in seq_len(m)) {
    gbm2_mod <- gbm3::gbm(
      formula   = Y ~ .,
      data      = data.frame(Y = Y[, j], X),
      distribution    = distribution
    )
    relinf <- summary(gbm2_mod, plot = FALSE)
    # summary gives data.frame(var, rel.inf)
    gbm2_imp_mat[relinf$var, j] <- relinf$rel_inf
  }
  gbm2_imp <- rowMeans(gbm2_imp_mat)
  
  # — Univariate RF —
  rf_uni_mat <- matrix(0, p, m, dimnames = list(colnames(X), colnames(Y)))
  for (j in seq_len(m)) {
    rf_uni <- randomForest::randomForest(x = X, y = Y[, j],
                                         ntree = 300,
                                         importance = TRUE)
    imp_j   <- randomForest::importance(rf_uni, type = 2)
    rf_uni_mat[, j] <- imp_j
  }
  rfuni_imp <- rowMeans(rf_uni_mat)
  
  if(type != "binary") {
    # — SPLS —
    spls_mod <- mixOmics::spls(X, as.matrix(Y))
    spls_imp <- apply(abs(spls_mod$loadings$X), 1, max)
    
    # — PMDCCA —
    pmdcca_mod <- PMA::CCA(x = X, z = Y,
                           typez = "standard",
                           typex = "standard")
    pmdcca_imp <- apply(abs(pmdcca_mod$u), 1, max)
    
    # — RGCCA —
    blocks <- list(X = X, Y = Y)
    rgcca_mod <- RGCCA::rgcca(blocks,
                              connection = matrix(c(0,1,1,0), 2, 2),
                              tau = c(1,1),
                              scheme = "centroid",
                              scale = TRUE)
    rgcca_imp <- apply(abs(rgcca_mod$a[[1]]), 1, max)
    ls <- list(
      IMD    = mrf_imd,
      MD     = mrf_md,
      XGBoost    = gb_imp,
      GBM = gbm2_imp,
      RFuni  = rfuni_imp,
      SPLS   = spls_imp,
      PMDCCA = pmdcca_imp,
      RGCCA  = rgcca_imp
    )
    
  } else {
    ls <- list(
      IMD    = mrf_imd,
      MD     = mrf_md,
      XGBoost    = gb_imp,
      GBM = gbm2_imp,
      RFuni  = rfuni_imp
    )
  }
  
  
  c(wl_m, ls)
}

evaluate_global <- function(imp_vec, true_idx) {
  # Binary truth vector for “signal” predictors
  truth_vec <- seq_along(imp_vec) %in% true_idx
  
  # 1) Precision‐Recall AUC
  pr <- pr.curve(
    scores.class0 = imp_vec[ truth_vec],  # positives
    scores.class1 = imp_vec[!truth_vec],  # negatives
    curve = FALSE
  )
  prauc_val <- pr$auc.integral
  
  rocobj  <- roc(truth_vec, imp_vec, quiet = TRUE)
  auc_val <- auc(rocobj)
  
  # 2) Top‐k metrics
  k           <- length(true_idx)
  ranked      <- order(imp_vec, decreasing = TRUE)
  topk        <- ranked[1:k]
  tp          <- sum(truth_vec[topk])
  recall_k    <- tp / k
  precision_k <- tp / k
  f1_k        <- if ((precision_k + recall_k) == 0) 0 else 
    2 * precision_k * recall_k / (precision_k + recall_k)
  
  c(PRAUC       = as.numeric(prauc_val),
    AUC = auc_val,
    Recall    = recall_k,
    Precision = precision_k,
    F1        = f1_k)
  
}

