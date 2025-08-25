### Load packages
library(PMA)
library(mixOmics)
library(RGCCA)
library(caret)
library(truncnorm)

# MRF
sim.fn.mrf3.m <- function(dat, ...){
  
  mrfinit  <- mrf3_init(dat$X, ...)
  
  imp_init <- plyr::llply(
    c("filter", "mixture", "test"),
    .fun = function(m) {
      w <- mrf3_vs(mrfinit, dat.list = dat$X, method = m)
      Reduce(c, w$weights)
    }
  )

  names(imp_init) <- c("filter", "mixture", "test")
  
  imp_a <- plyr::llply(
    c(0.1, 0.5),
    .fun = function(a){
      wla <- plyr::llply(
        mrfinit$weights_ls,
        .fun = function(wl) {
          Xa <- sapply(wl %>% purrr::map("X"), function(w) {
            new_w <- w^a
          }) %>% rowMeans()
          Ya <- sapply(wl %>% purrr::map("Y"), function(w) {
            new_w <- w^a
          }) %>% rowMeans()
          list(X = Xa, Y = Ya)
        }
      )
      wls <- plyr::llply(
        1:length(mrfinit$connection),
        .fun = function(i) {
          connection <- mrfinit$connection[[i]]
          names(wla[[i]]) <- rev(connection)
          wla[[i]]
        }
      )
      wl <- plyr::llply(
        names(dat$X),
        .fun = function(x) {
          k <- wls %>% purrr::map(x)
          k <- purrr::compact(k)
          Reduce("+", k)/length(k)
          
        }
      )
      names(wl) <- names(dat$X)
      mrfinit$weights <- wl

      mod_a <- plyr::llply(
        c("filter", "mixture", "test"),
        .fun = function(m) {
          w <- mrf3_vs(mrfinit, dat.list = dat$X, method = m)
          Reduce(c, w$weights)
        }
      )
      names(mod_a) <- c("filter", "mixture", "test")
      mod_a
    }
    
  )

  names(imp_a) <- paste0("IMD_a0", c(1,5))
  
  M <- c(list(IMD = imp_init), imp_a)
  
  
}

# Other benchmark methods
## Tuning functions
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

## Wrapper functions

sim.fn.other <- function(dat, cca.model = "spls", keep.list = NULL, ...){

  X <- dat$X %>% purrr::map(.,~scale(.))
  if(cca.model == "SPLS"){
    
    mod <- block.spls(
      X = X, indY = length(X),  keepX = keep.list, ncomp = 1,
      ...
    )
    
    feature.sel.scores <- mod$loadings %>% purrr::map(.,~.[,1])
    for (i in 1:length(X)) {
      names(feature.sel.scores[[i]]) <- colnames(X[[i]])
    }
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
  } 
  if(cca.model == "RGCCA"){
    mod <- rgcca(X, sparsity = 0.2)
    
    feature.sel.scores <- mod$a %>% purrr::map(.,~.[,1])
    for (i in 1:length(X)) {
      names(feature.sel.scores[[i]]) <- colnames(X[[i]])
    }
  }
        
  Reduce(c,feature.sel.scores)

}



