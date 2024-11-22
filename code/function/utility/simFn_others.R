##### Simulation Function for non-MRF based model #####
### Load packages
library(PMA)
library(mixOmics)
library(RGCCA)
library(caret)

spls.with.tune <- function(X, Y, keep.list,...){
  
  list.keepX <- keep.list[[1]]
  list.keepY <- keep.list[[2]]
  t <- tune.spls(X = X, Y = Y, ncomp = 1,
                        test.keepX = list.keepX,
                        test.keepY = list.keepY,
                        nrepeat = 1, folds = 3, 
                        measure = 'cor',... ) 
  
  sim.spls = spls(X = X, Y = Y, ncomp = 1, 
                  keepX = t$choice.keepX, keepY = t$choice.keepY)
  
  return(sim.spls)
  
}

o2pls.with.tune <- function(X, Y, keep.list,...){
  
  list.keepX <- keep.list[[1]]
  list.keepY <- keep.list[[2]]
  c <- crossval_sparsity(
    X = X, Y = Y, 
    n = 1, nx = 1, ny = 1,
    nr_folds = 5, 
    keepx_seq = list.keepX, 
    keepy_seq = list.keepY
  )
  
  sim.o2pls <- o2m(
    X = X, Y = Y, 
    n = 1, nx = 1, ny = 1,
    sparse = T,
    keepx = c$Best[1],
    keepy = c$Best[2],
    ...
  )
  
  return(sim.o2pls)
}

cca.with.tune <- function(X,Y,...){

  # perm.out <- CCA.permute(
  #   x = as.matrix(X),z = as.matrix(Y),trace = F,
  #   ...
  # )
  
  out <- CCA(
    x=as.matrix(X),z=as.matrix(Y),
    typex = "standard",
    typez = "standard",
    # penaltyx=perm.out$bestpenaltyx,
    # v=perm.out$v.init, 
    # penaltyz=perm.out$bestpenaltyz,
    xnames=colnames(X),
    znames=colnames(Y),
    trace = F,...
  )
  
  return(out)
}

sgcca.with.tune <- function(X,Y,...){
  
  tune <- rgcca_permutation(list(X = X, Y = Y),n_perms = 10,par_type = "sparsity",...)
  
  s <- rgcca(list(X = X, Y = Y), sparsity = tune$best_params)
  
  return(s)
}
  
  
sim.fn.other <- function(B = 200, sim.mod = "SimpleLinear", cca.model = "spls", sim.para, keep.list = NULL, parallel = F, ...){
  
  if(parallel) doParallel::registerDoParallel(10)
  
  mod <- plyr::llply(
    1:B,
    .fun = function(i){
      if(sim.mod == "SimpleLinear"){
        dat <- sim.linear3(
          n = sim.para$n,
          p = sim.para$p,
          rho = sim.para$rho,
          psel = sim.para$psel
        )
        X1 <- dat$X$X1; X2 <- dat$X$X2
        
      }else if(sim.mod == "Covariance"){
        dat <- sim.cov(
          n = sim.para$n,
          p = sim.para$p,
          ccor = sim.para$ccor,
          rho = sim.para$rho,
          psel = sim.para$psel
        )
        X1 <- data.frame(dat$X); X2 <- data.frame(dat$Y)
      }else if(sim.mod == "SimpleNonLinear"){
        dat <- sim.nonlinear2(
          n = sim.para$n,
          p = sim.para$p,
          rho = sim.para$rho,
          psel = sim.para$psel,
          j = sim.para$j
        )
        X <- dat$X %>% purrr::map(.,~scale(.))
        keep.list <-  rep(list(sim.para$psel),sim.para$j)
        names(keep.list) <- names(X)
      }else if(sim.mod == "NonLinearReg"){
        dat <- sim.nonlinear3(
          n = sim.para$n,
          p = sim.para$p,
          p.group = sim.para$p.group,
          q.group = sim.para$q.group,
          p.b = sim.para$p.b
        )
        X <- dat$X %>% purrr::map(.,~scale(.))
        keep.list <- list(X1 = 2 * sim.para$p.group, X2 = sim.para$q.group)
        names(keep.list) <- names(X)
      }else{
        dat <- sim.interaction1(
          n = sim.para$n,
          p = sim.para$p
        )
        X1 <- dat$X$X1; X2 <- dat$X$X2
      }
      
      
      tryCatch({
        if(cca.model == "spls"){
          

          mod <- block.spls(
            X = X, indY = length(X),  keepX = keep.list, ncomp = 1,
            ...
          )
          
          feature.sel.scores <- mod$loadings %>% purrr::map(.,~.[,1])
          for (i in 1:length(X)) {
            names(feature.sel.scores[[i]]) <- colnames(X[[i]])
          }
          feature.sel.labels <- feature.sel.scores %>% purrr::map(~ifelse(. != 0, 1, 0))
        } 
        if(cca.model == "o2pls"){
          mod <- o2pls.with.tune(
            X = X1, Y = X2, keep.list = keep.list,
            ...
          )
          
          feature.sel.scores <- list(X = mod$W., Y = mod$C.)
          feature.sel.labels <- feature.sel.scores %>% purrr::map(~ifelse(. != 0, 1, 0))
        }
        if(cca.model == "PMDCCA"){

          mod <- MultiCCA(
            X,
            ...
          )

          feature.sel.scores <- mod$ws %>% purrr::map(.,~.[,1]) 
          names(feature.sel.scores) <- names(X)
          for (i in 1:length(X)) {
            names(feature.sel.scores[[i]]) <- colnames(X[[i]])
          }
          feature.sel.labels <- feature.sel.scores %>% purrr::map(~ifelse(. != 0, 1, 0))
        } 
        if(cca.model == "sgcca"){
          mod <- rgcca(X, sparsity = 0.2)
          
          feature.sel.scores <- mod$a %>% purrr::map(.,~.[,1])
          for (i in 1:length(X)) {
            names(feature.sel.scores[[i]]) <- colnames(X[[i]])
          }
          feature.sel.labels <- feature.sel.scores %>% purrr::map(~ifelse(. != 0, 1, 0))
        }

        return(list(scores = feature.sel.scores, freq = feature.sel.labels))
      },error = function(e) {print(e);return(NULL)})
      
    },.parallel = parallel,.progress = "time"
  )
  mod <- mod %>% discard(is.null)
  dat_names <- names(mod[[1]]$scores)
  
  v <- plyr::llply(
    dat_names,
    .fun = function(i){
      mod %>% purrr::map("scores") %>% purrr::map(i) %>% Reduce(rbind,.)
    }
  )
  names(v) <- dat_names
  
  fr <- plyr::llply(
    dat_names,
    .fun = function(i){
      mod %>% purrr::map("freq") %>% purrr::map(i) %>% Reduce(rbind,.)
    }
  )
  names(fr) <- dat_names
  
  return(
    list(scores = v, freq = fr)
  )
}



