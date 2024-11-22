dir.base <- "."
dir.sim <- file.path(dir.base, "results/simulation/")
model <- "MRF-IMD-t300"
dir.sim.mrf <- file.path(dir.sim, "/model/", model, "two_omics")
dir.sim.mrf.ns <- file.path(dir.sim.mrf, "/Simple_Latent_NonLinear/")
dir.sim.mrf.nlreg <- file.path(dir.sim.mrf, "/NonLinearReg/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)
# Load related function
script <- list.files( 
  path = file.path(dir.base,"/code/function"),
  pattern = "[.]R$", 
  full.names = T,
  recursive = T
)
for (f in script) source(f)
library(multiMRF)


n <- c(100, 200, 200)
p <- c(200, 500, 1000)
psel = c(20, 30, 50)

parameter <- list(
  s1 = list(rho = 0)
)

mod = "SimpleNonLinear"
for(i in c(1:3)){
  sim.para <- list(
    n = n[i],
    p = p[i],
    rho = 0,
    psel = psel[i],
    j = 2
  )

  var.weights <- sim.fn.mrf3.m(
    B = 50,
    sim.mod = mod,
    sim.para = sim.para,
    parallel = T,
    connect_list = list(c("X1", "X2")),
    ntree = 300
  )

  true_label <- rep(c(rep(1,psel[i]), rep(0,p[i] - psel[i])), times = 2)

  plyr::l_ply(
    names(var.weights[[1]]),
    .fun = function(v) {

      d <- data.frame(
        dataset = c(rep("X", p[i]), rep("Y", p[i])),
        true = true_label
      )
      b <- var.weights$scores[[v]] %>%
        purrr::map(.,~t(.)) %>%
        Reduce(rbind,.) %>%
        data.frame() %>%
        rownames_to_column("variable") %>%
        dplyr::mutate(dataset = d$dataset, .before = 1)

      f <- var.weights$freq[[v]] %>%
        purrr::map(.,~t(.)) %>%
        Reduce(rbind,.) %>%
        data.frame() %>%
        rownames_to_column("variable") %>%
        dplyr::mutate(dataset = d$dataset, .before = 1)

      if(v == "thres1") {
        mm <- "-thres"
      }
      if(v == "thres1_se") {
        mm <- "-thres-se"
      }
      if(v == "thres2") {
        mm <- "-thres2"
      }
      if(v == "filter") {
        mm <- "-filter"
      }
      if(v == "mixture") {
        mm <- "-mixture"
      }
      if(v == "test") {
        mm <- "-test"
      }
      if(v == "mixture2") {
        mm <- "-mixture-norm"
      }
      
      m <- paste0(model, mm, "_",mod, "_n", n[i], "_p", p[i], "_psel", psel[i])

      writexl::write_xlsx(
        list(scores = b, freq = f, true = d),
        file.path(dir.sim.mrf.ns, paste0(m,".xlsx"))
      )

      message(m)
    }
  )


}


n <- c(100, 200, 200)
p <- c(200, 500, 1000)

parameter <- list(
  # s1 = list(p.group = 1, q.group = 5, p.b = 0)
  s2 = list(p.group = 2, q.group = 5, p.b = 0),
  s3 = list(p.group = 2, q.group = 10, p.b = 0)
)
mod = "NonLinearReg"
for(i in c(1:3)){
  for(para in parameter){
    sim.para <- list(
      n = n[i],
      p = p[i],
      p.group = para$p.group,
      q.group = para$q.group,
      p.b = para$p.b
    )
    
    var.weights <- sim.fn.mrf3.m(
      B = 50,
      sim.mod = mod,
      sim.para = sim.para,
      parallel = T,
      connect_list = list(c("X1", "X2")),
      ntree = 300
    )
    
    true_label <- c(
      c(rep(1, (2 * sim.para$p.group * sim.para$q.group)), rep(0, sim.para$p - (2 * sim.para$p.group * sim.para$q.group))),
      c(rep(1, sim.para$q.group), rep(0, sim.para$p - sim.para$q.group))
    )
    plyr::l_ply(
      names(var.weights[[1]]),
      .fun = function(v) {
        
        d <- data.frame(
          dataset = c(rep("X", p[i]), rep("Y", p[i])),
          true = true_label
        )
        b <- var.weights$scores[[v]] %>%
          purrr::map(.,~t(.)) %>%
          Reduce(rbind,.) %>%
          data.frame() %>%
          rownames_to_column("variable") %>%
          dplyr::mutate(dataset = d$dataset, .before = 1)
        
        f <- var.weights$freq[[v]] %>%
          purrr::map(.,~t(.)) %>%
          Reduce(rbind,.) %>%
          data.frame() %>%
          rownames_to_column("variable") %>%
          dplyr::mutate(dataset = d$dataset, .before = 1)
        
        if(v == "thres1") {
          mm <- "-thres"
        }
        if(v == "thres1_se") {
          mm <- "-thres-se"
        }
        if(v == "thres2") {
          mm <- "-thres2"
        }
        if(v == "filter") {
          mm <- "-filter"
        }
        if(v == "mixture") {
          mm <- "-mixture"
        }
        if(v == "test") {
          mm <- "-test"
        }
        if(v == "mixture2") {
          mm <- "-mixture-norm"
        }
        
        m <- paste0(model, mm,"_",mod, "_n", n[i], "_p", p[i], "_pgroup", para$p.group, "_qgroup", para$q.group, "_pb", para$p.b)
        
        writexl::write_xlsx(
          list(scores = b, freq = f, true = d),
          file.path(dir.sim.mrf.nlreg, paste0(m,".xlsx"))
        )
        
        message(m)
      }
    )
    
    
  }
  
}


dir.sim.mrf <- file.path(dir.sim, "/model/", model, "three_omics")
dir.sim.mrf.ns <- file.path(dir.sim.mrf, "/Simple_Latent_NonLinear/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

n <- c(100, 200, 200)
p <- c(200, 500, 1000)
psel = c(20, 30, 50)

mod = "SimpleNonLinear"
for(i in c(1:3)){
  sim.para <- list(
    n = n[i],
    p = p[i],
    rho = 0,
    psel = psel[i],
    j = 3
  )

  var.weights <- sim.fn.mrf3.m(
    B = 50,
    sim.mod = mod,
    sim.para = sim.para,
    parallel = T
  )

  true_label <- rep(c(rep(1,psel[i]), rep(0,p[i] - psel[i])), times = 3)
  d <- data.frame(
    dataset = c(rep("X", p[i]), rep("Y", p[i]), rep("Z", p[i])),
    true = true_label
  )

  plyr::l_ply(
    names(var.weights[[1]]),
    .fun = function(v) {

      b <- var.weights$scores[[v]] %>%
        purrr::map(.,~t(.)) %>%
        Reduce(rbind,.) %>%
        data.frame() %>%
        rownames_to_column("variable") %>%
        dplyr::mutate(dataset = d$dataset, .before = 1)

      f <- var.weights$freq[[v]] %>%
        purrr::map(.,~t(.)) %>%
        Reduce(rbind,.) %>%
        data.frame() %>%
        rownames_to_column("variable") %>%
        dplyr::mutate(dataset = d$dataset, .before = 1)

      if(v == "thres1") {
        mm <- "-thres"
      }
      if(v == "thres1_se") {
        mm <- "-thres-se"
      }
      if(v == "thres2") {
        mm <- "-thres2"
      }
      if(v == "filter") {
        mm <- "-filter"
      }
      if(v == "mixture") {
        mm <- "-mixture"
      }
      if(v == "test") {
        mm <- "-test"
      }
      if(v == "mixture2") {
        mm <- "-mixture-norm"
      }

      m <- paste0(model, mm,"_",mod, "_n", n[i], "_p", p[i], "_psel", psel[i])

      writexl::write_xlsx(
        list(scores = b, freq = f, true = d),
        file.path(dir.sim.mrf.ns, paste0(m,".xlsx"))
      )

      message(m)
    }
  )



}

