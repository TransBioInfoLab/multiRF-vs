#' MRF variable selection
#' 
#' @export
# This function performs variable selection for multi-omics data using MRF 

mrf3_vs <- function(mod, 
                    dat.list,
                    method = "filter",       # Selection method: "filter", "test", "mixture"
                    se = NULL,                  # Standard error multiplier for thresholding. Use if method = "filter". If NULL, select by best OOB error
                    c1 = "normal",           # Parameter for mixture model. Use if method = "mixture"
                    c2 = "normal",           # Parameter for mixture model. Use if method = "mixture"
                    level = 0.05,            # Significance level for t-test. Use if method = "test"
                    tscore = F,              # Whether to use t scores
                    re_weights = F,          # Whether to recalculate weights
                    re_fit = T,              # Whether to refit the model
                    ntree = 300,             # Number of trees for random forests
                    scale = F,               # Scale the data
                    k = 3,                   # Number of iterations for OOB error calculation
                    tol = 0.01,              # Tolerance for stopping OOB error calculation. Use if method = "filter"
                    iter = 1000,             # Max iterations for mixture model. Use if method = "mixture"
                    eps = 1e-05,             # Convergence threshold for mixture model. Use if method = "mixture"
                    normalized = T,          # Normalize weights after selection
                    select = "ALL",          # Subset of variables to select
                    ...){                    # Additional arguments in rfsrc package
   
  
  weights <- mod$weights
  new_dat <- dat.list
  dat_names <- names(new_dat)
  connect_list <- mod$connection

  message("Variable selection..")
  
  # T-score method
  if(method == 'test'|tscore) {

    thres <- test_fn( 
      wl = mod$weights_ls,
      connection = connect_list,
      dat_names = dat_names,
      sig.thres = level)
    
  }
  
  # Mixture model method
  if(method == "mixture") {
    if(tscore) ts <- thres$ts
    
    thres <- plyr::llply(
      names(weights),
      .fun = function(t) {
        if(tscore) {
          ww0 <- ts[[t]]
        } else {
          ww0 <- weights[[t]]
        }
        
        param <- get_param(ww0)
        e <- em(ww0, p = param$p, mu = param$mu, sigma = param$sigma, iter = iter, eps = eps, c1 = c1, c2 = c2)
        cbind(e$par$post[,1]/rowSums(e$par$post),
              e$par$post[,2]/rowSums(e$par$post),
              e$par$post[,3]/rowSums(e$par$post))
        
      }
    )
    
    names(thres) <- names(weights)

  }
  
  # Filter method (default)
  if(method == "filter") {
    weights_norm <- purrr::map(weights, ~./sqrt(sum(.^2)))
    thres <- choose_thres2(
      weights,
      connection = connect_list,
      new_dat = new_dat, 
      yprob = mod$yprob, 
      ntree = mod$ntree, 
      type = mod$type,
      oob_init = mod$oob_err, 
      k = k, 
      select = select,
      tol = tol,
      se = se,
      ...)
    
  }

  # Update weights based on thresholding results
  weights_new <- plyr::llply(
    dat_names,
    .fun = function(i) {
      w <- weights[[i]]
      if(!method %in% c("mixture", "test") ) {
        if(select != 'ALL') {
          
          if(i %in% select) {
            w[w < thres[i]] <- 0
            w[w > thres[i]] <- w[w > thres[i]] - thres[i]
          } 
        } else {
          w[w < thres[i]] <- 0
          w[w > thres[i]] <- w[w > thres[i]] - thres[i]
        }
      } else {
        if(method == "mixture") {
          t <- ifelse(thres[[i]][,1] < level,1,0)
        } else {
          t <- thres$keep_idx[[i]]
          
        }

        if(select != 'ALL') {
          
          if(i %in% select) {
            
            w <- w * t
          } 
        } else {

          w <- w * t
        }
        
      }
      
      
      if(normalized) {
        w <- w/sqrt(sum(w^2))
      }
      w
    }
  )

  names(weights_new) <- dat_names
  mod$weights <- weights_new
  message("Refit model..")
    
  # Update datasets based on selected variables
  new_dat <- plyr::llply(
    dat_names,
    .fun = function(i){
      if(select != 'ALL') {
        if(i %in% select) {
          new_dat[[i]] <- new_dat[[i]][,weights_new[[i]] != 0]
          new_dat[[i]]
        } else {
          new_dat[[i]]
        }
      } else {
        new_dat[[i]] <- new_dat[[i]][,weights_new[[i]] != 0]
        new_dat[[i]]
      }
    }
  )
  names(new_dat) <- dat_names
  m <- purrr::map(weights_new, ~.[.>0])
  
  # Scale data if specified
  if(scale) {
    new_dat2 <- purrr::map(new_dat, ~scale(.))
  } else {
    new_dat2 <- new_dat
  }
    
  # Optionally refit the model
  if(!re_weights) {
    m <- NULL
  } 
  if(re_fit) {
    refit <- fit_multi_rfsrc(new_dat2, connect_list = connect_list, ntree = ntree, 
                             type = mod$type, var.wt = m, yprob = mod$yprob, ...)
    
    oob_err <- purrr::map(refit, ~get_r_sq(.))
    oob_err <- Reduce("+", oob_err)
    
    mod$oob_err <- oob_err
    mod$mod <- refit
  }
  
  # Finalize model
  mod$dat.list <- new_dat2
  mod$thres <- thres
  mod$ntree <- ntree
  
  class(mod) <- c(class(mod), "vs")
  
  mod
  
}

# Function to compute node distributions
get_d <- function(node) {
  t <- table(node)
  ind <- order(as.numeric(names(t)), decreasing = T)
  t <- t[ind]

  ld <- t
  Ld <- c(1, cumsum(t)[-length(cumsum(t))])
  
  list(ind = as.numeric(names(t)), ld = ld, Ld = Ld, tb = t)
}

# Function to calculate probability for given parameters
calc_prob <- function(prob, LD, ld) {
  
  ((1-prob)^LD)*(1-(1-prob)^ld)
} 
 
# Function to calculate mean and standard deviation of matrix columns
calc_mean_sd <- function(mat) {
  
  p <- nrow(mat)
  # p <- sum(rowMeans(mat != 0) < mean(mat != 0))
  mean_mat <- apply(
    mat, 
    2,
    function(m) {
      d <- get_d(m)
      # p <- nrow(mat) - d$t[length(d$t)]
      ex <- sapply(1:(length(d$ind) - 1), function(i) {
        pr <- calc_prob(1/p, d$Ld[i], d$ld[i])
        d$ind[i] * pr
      })
      
      mm <- sum(ex)
      
      exx <- sapply(1:(length(d$ind) - 1), function(i) {
        pr <- calc_prob(1/p, d$Ld[i], d$ld[i])
        (d$ind[i] - mm)^2 * pr
      })
      
      s <- sqrt(sum(exx))
      
      list(mean = mm, sd = s)
    }
  )

  mm <- unlist(purrr::map(mean_mat, "mean"))
  s <- unlist(purrr::map(mean_mat, "sd"))
  
  
  list(mean = mean(mm),
       se = s/sqrt(length(s)))
  
}

# Function to choose optimal thresholds by minimizing OOB error
choose_thres2 <- function(weights, connection, new_dat, yprob, ntree, type, oob_init, k = 3, tol = 0.01, select = "ALL", se = se,  ...) {
  
  if(is.null(se)) {
    num <- seq(.8,3.1,by = 0.1)
    oob <- c()
    i <- 1
    while(i <= length(num)) {
      t <- num[i]
      
      weights_new <- plyr::llply(
        weights,
        .fun = function(w) {
          
          thres <- t * sd(w)
          if(select != 'ALL') {
            
            if(i %in% select) {
              w[w < thres] <- 0
              w[w > thres] <- w[w > thres] - thres
            } 
          } else {
            w[w < thres] <- 0
            w[w > thres] <- w[w > thres] - thres
          }
          
          w
          
        }
      )
      
      new_dat2 <- plyr::llply(
        names(new_dat),
        .fun = function(i){
          if(select != 'ALL') {
            if(i %in% select) {
              new_dat[[i]] <- new_dat[[i]][,weights_new[[i]] != 0]
              new_dat[[i]]
            } else {
              new_dat[[i]]
            }
          } else {
            new_dat[[i]] <- new_dat[[i]][,weights_new[[i]] != 0]
            new_dat[[i]]
          }
        }
      )
      names(new_dat2) <- names(new_dat)
      
      oob_err0 <- plyr::laply(1:k,
                              .fun = function(k) {
                                refit <- fit_multi_rfsrc(new_dat2, connect_list = connection, 
                                                         ntree = ntree, type = type, yprob = yprob, 
                                                         var.wt = NULL,
                                                         forest.wt = "oob")
                                
                                oob_err <- purrr::map(refit, ~get_r_sq(.))
                                Reduce("+", oob_err)/length(new_dat2)
                              })
      
      oob_err <- mean(unlist(oob_err0))
      oob <- c(oob, oob_err)
      oob_init <- oob_err
      i <- i + 1
    }
    diff <- abs(diff(oob[-1]))
    t <- sort(num[-1][which(diff < tol)])[1]
  } else { t <- se}
  

 
  
  message("Choose ", t, " times sd")
  thres <- plyr::llply(
    weights,
    .fun = function(w) {
      
      t * sd(w)
    }
  )
  
  unlist(thres)
  
}

# Function for t-test on weights
test_fn <- function(wl, connection, dat_names, sig.thres = 0.05) {

  res <- plyr::llply(1:length(wl),
              .fun = function(i) {
                w <- wl[[i]]
                conn <- rev(connection[[i]])
                pls <- plyr::llply(
                  names(w[[1]]),
                  .fun = function(j) {
                    mat <- Reduce(cbind, purrr::map(w, j) )

                    mu <- mean(mat)

                    se <- apply(mat, 1, function(k){
                      if(length(k[k!=0]) > 1) {
                        sd(k)/sqrt(length(k) - 1)
                      } else {
                        -1
                      }

                    } )

                    z <- ifelse(se != -1, (rowMeans(mat) - mu)/se, 0)

                    p <- pt(z,df = nrow(mat) - 1 ,lower.tail = F)

                    keep <- ifelse(p < sig.thres, 1, 0)

                    list(keep_idx = keep,
                         pval = p,
                         ts = z)
                  }
                )

                p <- purrr::map(pls, "pval")
                keep_idx <- purrr::map(pls, "keep_idx")
                ts <- purrr::map(pls, "ts")
                names(p) <- names(keep_idx) <- names(ts) <- conn


                list(keep_idx = keep_idx,
                     pval = p, 
                     ts = ts)
              })

  keep_idx <- purrr::map(res, "keep_idx")
  p <- purrr::map(res, "pval")
  trans_score <- ts <- purrr::map(res, "ts")
  trans_score <- plyr::llply(dat_names,
                             .fun = function(d) {
                               ts0 <- purrr::map(ts, d) %>% purrr::compact()
                               Reduce("+", ts0)/length(ts0)
                             })
  names(trans_score) <- dat_names

  keep_res <- plyr::llply(dat_names,
              .fun = function(d) {
                keep <- Reduce(cbind, purrr::map(keep_idx, d) )
                ts <- Reduce(cbind, purrr::map(ts, d))
                if(is.null(ncol(ts))) ts <- ts else ts <- colMeans(ts)


                if(!is.null(ncol(keep))) {
                  keep <- rowMeans(keep) >= 0.5
                } else {
                  keep <- keep == 1
                }
                keep_idx <- ifelse(keep, 1, 0)
                p <- Reduce(cbind, purrr::map(p, d))
                list(keep_idx = keep_idx,
                     pval = p,
                     ts = ts)

              })

  names(keep_res) <- dat_names
  keep_idx <- purrr::map(keep_res, "keep_idx")
  p <- purrr::map(keep_res, "pval")
  ts <- purrr::map(keep_res, "ts")

  list(keep_idx = keep_idx,
       pval = p,
       ts = ts,
       trans_score = trans_score)

}
