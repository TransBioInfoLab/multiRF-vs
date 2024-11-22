library(truncnorm)
sim.fn.mrf3.m <- function(B = 200, sim.mod = "SimpleLinear", sim.para, parallel = T, j = 2, ...){
  
  if(parallel) doParallel::registerDoParallel(3)
  
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
        
      }else if(sim.mod == "Covariance"){
        dat <- sim.cov(
          n = sim.para$n,
          p = sim.para$p,
          ccor = sim.para$ccor,
          rho = sim.para$rho,
          psel = sim.para$psel
        )
      }else if(sim.mod == "SimpleNonLinear"){

        dat <- sim.nonlinear2(
          n = sim.para$n,
          p = sim.para$p,
          rho = sim.para$rho,
          psel = sim.para$psel,
          j = sim.para$j
        )
      }else if(sim.mod %in% c("NonLinearReg", "high_dimentional")){
        dat <- sim.nonlinear3(
          n = sim.para$n,
          p = sim.para$p,
          p.group = sim.para$p.group,
          q.group = sim.para$q.group,
          p.b = sim.para$p.b
        )
      }else{
        dat <- sim.interaction1(
          n = sim.para$n,
          p = sim.para$p
        )

      }
      
      tryCatch({

        modinit <- mrf3_init(
          dat.list = dat$X, 
          parallel = parallel,
          ...
        )
        
        mod <- list(
          # test =  mrf3_vs(modinit, dat.list = dat$X, method = "test")
          # thres1 = mrf3_vs(modinit, method = "thres", se = 0),
          # thres1_se = mrf3_vs(modinit, method = "thres", se = 2),
          # thres2 = mrf3_vs(modinit, method = "thres", use_distribution = F),
          # mixture = mrf3_vs(modinit, method = "mixture"),
          # mixture2 = mrf3_vs(modinit, dat.list = dat$X, method = "mixture", c1 = "normal")
          filter = mrf3_vs(modinit, dat.list = dat$X, method = "filter")
        ) 

        feature.sel.scores <- mod %>% purrr::map(., ~.$weights)
        feature.sel.labels <- plyr::llply(
          feature.sel.scores,
          .fun = function(mm) {
            mm %>% purrr::map(~ifelse(. > 0, 1, 0))
          }
        )
          
        return(list(scores = feature.sel.scores, freq = feature.sel.labels))
      },error = function(e) {print(e);return(NULL)})
      
    },.parallel = F
  )
  mod <- mod %>% purrr::discard(is.null)

  dat_names <- names(mod[[1]]$scores[[1]])

  v <- plyr::llply(
    names(mod[[1]]$scores),
    .fun = function(j) {
      mod2 <- mod %>% purrr::map("scores") %>% purrr::map(j)
      res <- plyr::llply(
        dat_names,
        .fun = function(i){
          mod2 %>% purrr::map(i) %>% Reduce(rbind,.)
        }
      )
      names(res) <- dat_names
      res
    }
  )
    
  names(v) <- names(mod[[1]]$scores)
  
  fr <- plyr::llply(
    names(mod[[1]]$freq),
    .fun = function(j) {
      mod2 <- mod %>% purrr::map("freq") %>% purrr::map(j)
      res <- plyr::llply(
        dat_names,
        .fun = function(i){
          mod2 %>% purrr::map(i) %>% Reduce(rbind,.)
        }
      )
      names(res) <- dat_names
      res
    }
  )
  names(fr) <- names(mod[[1]]$freq)
  
  return(
    list(scores = v, freq = fr)
  )
}
