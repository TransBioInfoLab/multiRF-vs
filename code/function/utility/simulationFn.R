##### Simulation Model #####
library(MASS)
library(corpcor)


####### Simple linear latent function #######



sim.linear1 = function(n, p, sd.w1 = 1, sd.w2 = 1, sigma = 0.3, psel = NULL){
  w1.0 = rnorm(n, mean = 0, sd = sd.w1)
  w1 = w1.0/sum(w1.0)
  w2.0 = rnorm(n, mean = 0, sd = sd.w2)
  w2 = w2.0/sum(w2.0)
  u1 = c(rep(1, 20), rep(-1, 20), rep(0, 60))
  u2 = c(rep(-1, 10), rep(1, 10), rep(-1, 10), rep(1, 10), rep(0, 60))
  v1 = c(rep(0, 60), rep(-1, 20), rep(1, 20))
  v2 = c(rep(0, 60), rep(1, 10), rep(-1, 10), rep(1, 10), rep(-1, 10))
  e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma^2, nrow = p))
  e2 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma^2, nrow = p))
  X1 = w1 %*% t(u1) + w2 %*% t(u2) + e1
  X2 = w1 %*% t(v1) + w2 %*% t(v2) + e2
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  return(list(X = list(X1 = X1, X2 = X2)))
}

#l1 = sim.linear1(50, sd.w1 = 1, sd.w2 = 1, sigma = 0.3)

sim.linear2 = function(n, p, j = 2, mu.sd = 2, sigma = 0.3, psel = 20){
  X = list()
  w = list()
  mu = rnorm(n, mean = 0, sd = mu.sd) 
  A = c("X", "Y", "Z", "V", "W")
  for(i in 1:j){
    w0 = runif(psel, -1, 1) 
    w[[i]] = c(w0/sum(w0), rep(0, p - psel))
    e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
    X[[i]] = data.frame(scale(mu %*% t(w[[i]]) + e1))
    colnames(X[[i]]) = sapply(1:p, function(m) paste0(A[i],m))
  }
  names(X) = sapply(1:j, function(i) paste0("X",i))
  return(list(X = X, weight = w))
}

#l2 = sim.linear2(n = 50, p = 100, j = 2, mu.sd = 2, sigma = 0.3)

### Gaussian latent variable and orthogonal weights

# sim.linear3 = function(n, p, j = 2, mu.sd = 2, latent = 5, rho = 0.9, sigma = 0.3, psel = 20){
#   X = list()
#   w = list()
#   Corr = diag(1, nrow = latent)
#   Corr[upper.tri(Corr)] = rho
#   Corr[lower.tri(Corr)] = rho
#   std = rep(mu.sd, latent)
#   Sigma = std %*% t(std) * Corr
#   mu = mvrnorm(n, mu = rep(0, latent), Sigma)
#   A = c("X", "Y", "Z", "V", "W")
#   for(i in 1:j){
#     w0 = runif(psel, -1, 1) 
#     w[[i]] = c(w0/sqrt(sum(w0^2)), rep(0, p - psel))
#     e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma^2, nrow = p))
#     X0 = lapply(1:latent, function(k) mu[,k] %*% t(w[[i]]))
#     X[[i]] = data.frame(scale( Reduce("+",X0)  + e1))
#     colnames(X[[i]]) = sapply(1:p, function(m) paste0(A[i],m))
#   }
#   names(X) = sapply(1:j, function(i) paste0("X",i))
#   return(list(X = X, weight = w))
# }

sim.linear3 = function(n, p, j = 2, mu.sd = 2, rho = 0, sigma = 0.3, psel = 20){
  X = list()
  w = list()
  Corr = diag(1, nrow = p)
  Corr[upper.tri(Corr)] = rho
  Corr[lower.tri(Corr)] = rho
  std = rep(sigma, p)
  Sigma = std %*% t(std) * Corr
  mu = rnorm(n, mean = 0, sd = mu.sd) 
  A = c("X", "Y", "Z", "V", "W")
  for(i in 1:j){
    w0 = runif(psel, -1, 1) 
    w[[i]] = c(w0/sqrt(sum(w0^2)), rep(0, p - psel))
    e1 = mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    X[[i]] = data.frame(scale(mu %*% t(w[[i]]) + e1))
    colnames(X[[i]]) = sapply(1:p, function(m) paste0(A[i],m))
  }
  names(X) = sapply(1:j, function(i) paste0("X",i))
  return(list(X = X, weight = w))
}


covarianceToeplitz <- function(p,rho = 0.9){
  Sigma <- matrix(0, nrow = p, ncol =p)
  for (i in 1:p){
    for (j in 1:p){
      power <- abs(i-j)
      Sigma[i,j] = rho^power
    }
  }
  return(Sigma)
}

sim.cov <- function(n, p, ccor = 0.9, rho = 0.9, psel = 20){
  
  SigmaX <- covarianceToeplitz(p, rho = rho)
  SigmaY <- covarianceToeplitz(p, rho = rho)

  w1 <- c(runif(psel),rep(0, p-psel))
  w2 <- c(runif(psel),rep(0, p-psel))
  SigmaXY <- ccor*(SigmaX%*%w1%*%t(w2)%*%SigmaY)
  SigmaYX <- t(SigmaXY)
  ## Generate the data-sets
  Sigma1 <- cbind(SigmaX, SigmaXY)
  Sigma2 <- cbind(SigmaYX, SigmaY)
  Sigma <- rbind(Sigma1, Sigma2)
  Sigma <- make.positive.definite(Sigma)
  XY <- mvrnorm(n = n, mu = rep(0,2*p), Sigma = Sigma)
  X <- XY[,1:p]
  Y <- XY[,-(1:p)]
  colnames(X) = sapply(1:p, function(m) paste0("X",m))
  colnames(Y) = sapply(1:p, function(m) paste0("Y",m))
  return(list("X" = X, "Y" = Y))
}

sim.linear4 = function(n, p, p.group = 3, q.group = 3, p.b = 0, psel = NULL){
  y.latent = q.group*3 
  latent = y.latent + p.b
  x = sapply(1:latent, function(i) {
    #rnorm(n = n, mean = 0, sd = 2)})
    runif(n = n, min = 0, max = 1)})
  colnames(x) = paste("x", 1:latent, sep = ".")
  
  ## simulate groups of correlated variables based on functional variables
  
  g = lapply(1:y.latent, function(i) {
    v = sapply(1:p.group, function(j) {
      if(p.group > 1){
        x[, i] + 0.01 + 0.5 * (j - 1) / (p.group - 1)  *
          rnorm(n = n, mean = 0, sd = 0.3)
      }else{
        x[, i] + 0.01 +
          rnorm(n = n, mean = 0, sd = 0.3)
      }
    })
    colnames(v) = paste("X", i, 1:p.group, sep = ".")
    return(v)
  })
  g = do.call(cbind, g)
  
  if(p.b > 0){
    g.0 = lapply((y.latent+1):latent, function(i) {
      v = sapply(1:p.group, function(j) {
        if(p.group > 1){
          x[, i] + 0.01 + 0.5 * (j - 1) / (p.group - 1)  *
            rnorm(n = n, mean = 0, sd = 0.3)
        }else{
          x[, i] + 0.01 + 
            rnorm(n = n, mean = 0, sd = 0.3)
        }
      })
      colnames(v) = paste("X0", i, 1:p.group, sep = ".")
      return(v)
    })
    g.0 = do.call(cbind, g.0)
  }else g.0 <- NULL
  
  
  if (!is.null(g.0)){
    if(ncol(g) + ncol(g.0) >= p) {
      warning("No additional independent variables are simulated!")
      ind = NULL
    } else {
      ind = mvrnorm(n = n, mu = rep(0, p - (ncol(g) + ncol(g.0))), 
                    Sigma = diag(.3^2,  p - (ncol(g) + ncol(g.0))))
      colnames(ind) = paste("ind", 1:ncol(ind), sep = ".")
    }
  } else {
    if(ncol(g) >= p) {
      warning("No additional independent variables are simulated!")
      ind = NULL
    } else {
      ind = mvrnorm(n = n, mu = rep(0, p - ncol(g)), 
                    Sigma = diag(.3^2,  p - ncol(g)))
      colnames(ind) = paste("ind", 1:ncol(ind), sep = ".")
    }
  }
  
  Y = sapply(1:q.group, function(i) {
    v = sapply(1:3, function(j) {
      x[,3*(i - 1) + j]
    })
    y = get_kennel_fn2(v[,1], v[,2], v[,3])
    return(y)
  })
  colnames(Y) = sapply(1:q.group, function(i) paste("Y", i, sep = "."))
  
  Y.ind = mvrnorm(n = n, mu = rep(0, p - q.group), 
                  Sigma = diag(.3^2,  p - q.group))
  colnames(Y.ind) = paste("Y.ind", 1:ncol(Y.ind), sep = ".")
  if (!is.null(g.0)){X1 = data.frame(scale(cbind(g, g.0, ind)))} else X1 = data.frame(scale(cbind(g, ind)))
  X2 = data.frame(scale(cbind(Y, Y.ind)))
  
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  
  return(list(X = list(X1 = X1, X2 = X2)))
}


######## nonlinear latent function #########


get_kennel_fn = function(x1, x2, x3){
  y = 0.25 * exp( 4 * x1) + 4/(1+exp(-20*(x2 - 0.5))) + 
    3 * x3 + rnorm(n = length(x1), mean = 0, sd = 0.2)
  y
}
get_kennel_fn2 = function(x1, x2, x3){
  y = 0.25 * x1 + x2 + 
    3 * x3 + rnorm(n = length(x1), mean = 0, sd = 0.2)
  y
}
get_kennel_fn3 = function(x1, x2, x3, x4){
  y = 0.25 * exp( 4 * x1) + 4/(1+exp(-20*(x2 - 0.5))) + 
    x3 * x4 + rnorm(n = length(x1), mean = 0, sd = 0.2)
  y
}
get_kennel_fn4 = function(x1, x2){
  y = 0.25 * exp( 4 * x1) + 4/(1+exp(-20*(x2 - 0.5))) + rnorm(n = length(x1), mean = 0, sd = 0.2)
  y
}

sim.nonlinear1 = function(n,p, mu.sd = 2, sigma = 0.2, p.group = 30, q.group = 4, psel = NULL){
  X1.0 = sapply(1:p.group, function(i) rnorm(n, 0, mu.sd))
  X2.0 = sapply(1:(p.group - q.group), function(i) rnorm(n, 0, mu.sd))
  Y = sapply(1:q.group, function(i) {
    v = sapply(1:3, function(j) {
      X1.0[,3*(i - 1) + j]
    })
    y = get_kennel_fn2(v[,1], v[,2], v[,3])
    return(y)
  })
  e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  e2 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  X1 = cbind(X1.0, matrix(0, ncol = p-p.group, nrow = n)) + e1
  X2 = cbind(Y, X2.0, matrix(0, ncol = p-p.group, nrow = n)) + e2
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  return(list(X = list(X1 = X1, X2 = X2)))
}

sim.nonlinear4 = function(n,p, mu.sd = 2, sigma = 0.2, p.group = 30, q.group = 4, psel = NULL){
  X1.0 = sapply(1:p.group, function(i) rnorm(n, 0, mu.sd))
  X2.0 = sapply(1:(p.group - q.group), function(i) rnorm(n, 0, mu.sd))
  #Y = sapply(1:q.group, function(i) {
  #  v = sapply(1:3, function(j) {
  #    X1.0[,3*(i - 1) + j]
  #  })
  #  y = get_kennel_fn(v[,1], v[,2], v[,3])
  #  return(y)
  #})
  X2.1 = 0.5*exp(X1.0[,1]) + X1.0[,2]^2 - X1.0[,7]
  X2.2 = 3*sin(pi*X1.0[,3]/2) - 2*X1.0[,4]^4 
  X2.3 = abs(X1.0[,5]) + 0.25*X1.0[,6]^3
  X2.4 = 4/(1 + exp(X1.0[,8]-20)) 
  e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  e2 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  X1 = cbind(X1.0, matrix(0, ncol = p-p.group, nrow = n)) + e1
  X2 = cbind(X2.1,X2.2,X2.3,X2.4, X2.0, matrix(0, ncol = p-p.group, nrow = n)) + e2
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  return(list(X = list(X1 = X1, X2 = X2)))
}

sim.nonlinear5 = function(n,p, mu.sd = 2, sigma = 0.2, p.group = 30, q.group = 4, psel = NULL){
  X1.0 = sapply(1:p.group, function(i) rnorm(n, 0, mu.sd))
  X2.0 = sapply(1:(p.group - q.group), function(i) rnorm(n, 0, mu.sd))
  #Y = sapply(1:q.group, function(i) {
  #  v = sapply(1:3, function(j) {
  #    X1.0[,3*(i - 1) + j]
  #  })
  #  y = get_kennel_fn(v[,1], v[,2], v[,3])
  #  return(y)
  #})
  X2.1 = X1.0[,1]^2 - X1.0[,2]
  X2.2 = exp(X1.0[,3]) - 2*X1.0[,4]^2 
  X2.3 = 0.25*X1.0[,5]^3
  X2.4 = 4/(1 + exp(X1.0[,6]-20)) 
  e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  e2 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  X1 = cbind(X1.0, matrix(0, ncol = p-p.group, nrow = n)) + e1
  X2 = cbind(X2.1,X2.2,X2.3,X2.4, X2.0, matrix(0, ncol = p-p.group, nrow = n)) + e2
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  return(list(X = list(X1 = X1, X2 = X2)))
}


sim.nonlinear2 = function(n, p, j = 2, mu.sd = 2, rho = 0, sigma = 0.3, psel = 20){
  
  Corr = diag(1, nrow = p)
  Corr[upper.tri(Corr)] = rho
  Corr[lower.tri(Corr)] = rho
  std = rep(sigma, p)
  Sigma = std %*% t(std) * Corr
  mu0 <- rnorm(n, mean = 0, sd = mu.sd)
  mu = list(mu0^2,
            4/exp(mu0),
            mu0)
  w <- list()
  X <- list()
  
  A <- c("X", "Y", "Z")       
  for(i in 1:j){
    w0 = runif(psel, -1, 1) 
    w[[i]] = c(w0/sqrt(sum(w0^2)), rep(0, p - psel))
    e1 = mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    X[[i]] = data.frame(scale(mu[[i]] %*% t(w[[i]]) + e1))
    colnames(X[[i]]) = sapply(1:p, function(m) paste0(A[i],m))
  }
  
  names(X) = sapply(1:j, function(i) paste0("X",i))
  return(list(X = X, weight = w))
}

sim.nonlinear3 = function(n, p, p.group = 3, q.group = 3, p.b = 0, psel = NULL){
  y.latent = q.group*3 
  latent = y.latent + p.b
  x = sapply(1:latent, function(i) {
       #rnorm(n = n, mean = 0, sd = 2)})
    runif(n = n, min = 0, max = 1)})
  colnames(x) = paste("x", 1:latent, sep = ".")
  
  ## simulate groups of correlated variables based on functional variables
  
  g = lapply(1:y.latent, function(i) {
    v = sapply(1:p.group, function(j) {
      if(p.group > 1){
        x[, i] + 0.01 + 0.5 * (j - 1) / (p.group - 1)  *
          rnorm(n = n, mean = 0, sd = 0.3)
      }else{
        x[, i] + 0.01 +
          rnorm(n = n, mean = 0, sd = 0.3)
      }
    })
    colnames(v) = paste("X", i, 1:p.group, sep = ".")
    return(v)
  })
  g = do.call(cbind, g)
  
  if(p.b > 0){
    g.0 = lapply((y.latent+1):latent, function(i) {
      v = sapply(1:p.group, function(j) {
        if(p.group > 1){
          x[, i] + 0.01 + 0.5 * (j - 1) / (p.group - 1)  *
            rnorm(n = n, mean = 0, sd = 0.3)
        }else{
          x[, i] + 0.01 + 
            rnorm(n = n, mean = 0, sd = 0.3)
        }
      })
      colnames(v) = paste("X0", i, 1:p.group, sep = ".")
      return(v)
    })
    g.0 = do.call(cbind, g.0)
  }else g.0 <- NULL
  
  
  if (!is.null(g.0)){
    if(ncol(g) + ncol(g.0) >= p) {
      warning("No additional independent variables are simulated!")
      ind = NULL
    } else {
      ind = mvrnorm(n = n, mu = rep(0, p - (ncol(g) + ncol(g.0))), 
                    Sigma = diag(.3^2,  p - (ncol(g) + ncol(g.0))))
      colnames(ind) = paste("ind", 1:ncol(ind), sep = ".")
    }
  } else {
    if(ncol(g) >= p) {
      warning("No additional independent variables are simulated!")
      ind = NULL
    } else {
      ind = mvrnorm(n = n, mu = rep(0, p - ncol(g)), 
                    Sigma = diag(.3^2,  p - ncol(g)))
      colnames(ind) = paste("ind", 1:ncol(ind), sep = ".")
    }
  }
  
  Y = sapply(1:q.group, function(i) {
    v = sapply(1:2, function(j) {
      x[,2*(i - 1) + j]
    })
    y = get_kennel_fn4(v[,1], v[,2])
    return(y)
  })
  colnames(Y) = sapply(1:q.group, function(i) paste("Y", i, sep = "."))
  
  Y.ind = mvrnorm(n = n, mu = rep(0, p - q.group), 
                  Sigma = diag(.3^2,  p - q.group))
  colnames(Y.ind) = paste("Y.ind", 1:ncol(Y.ind), sep = ".")
  if (!is.null(g.0)){X1 = data.frame(scale(cbind(g, g.0, ind)))} else X1 = data.frame(scale(cbind(g, ind)))
  X2 = data.frame(scale(cbind(Y, Y.ind)))
  
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  
  return(list(X = list(X1 = X1, X2 = X2)))
}


sim.nonlinear6 = function(n,p, mu.sd = 2, rho = 0, sigma = 0.3, psel = 20){
  w1.0 = runif(psel, -1, 1) 
  w1 = c(w1.0/sqrt(sum(w1.0^2)), rep(0, p - psel))
  w2.0 = runif(psel, -1, 1)
  w2 = c(w2.0/sqrt(sum(w2.0^2)), rep(0, p - psel))
  Corr = diag(1, nrow = p)
  Corr[upper.tri(Corr)] = rho
  Corr[lower.tri(Corr)] = rho
  std = rep(sigma, p)
  Sigma = std %*% t(std) * Corr
  mu = rnorm(n, mean = 0, sd = mu.sd) 
  # mu1 = mu^2
  # mu2 = 4/exp(mu)
  e1 = mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  e2 = mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X1 = data.frame(scale((mu %*% t(w1))^2 + e1))
  X2 = data.frame(scale(4/exp(mu %*% t(w2)) + e2))
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  return(list(X = list(X1 = X1, X2 = X2), weights = cbind(w1,w2)))
}


######## interaction function #########


sim.interaction1 = function(n, p, distribution = "normal", para = c(1,2), sigma = 0.2,
                            psel = NULL){
  if(distribution == 'gamma'){
    X1.0 = sapply(1:8, function(k) rgamma(n, para[1], para[2]))
  }
  if(distribution == 'beta'){
    X1.0 = sapply(1:8,  function(k) rbeta(n, para[1], para[2]))
  }
  if(distribution == 'uniform'){
    X1.0 = sapply(1:8, function(k) runif(n, para[1], para[2]))
  }
  if(distribution %in% c('gaussian','normal')){
    X1.0 = sapply(1:8, function(k) rnorm(n, para[1], para[2]^2))
  }
  X2.1 = X1.0[,1] * X1.0[,2]
  X2.2 = X1.0[,3] * X1.0[,4]
  X2.3 = X1.0[,5] * X1.0[,6]
  X2.4 = X1.0[,7] * X1.0[,8]
  e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  e2 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  X1 = scale(cbind(X1.0, matrix(0, ncol = p-8, nrow = n)) + e1)
  X2 = scale(cbind(X2.1,X2.2, X2.3, X2.4, matrix(0, ncol = p-4, nrow = n)) + e2)
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  return(list(X = list(X1 = X1, X2 = X2)))
}



sim.interaction2 = function(n, p, distribution = "normal", para = c(1,2), sigma = 0.2,
                            psel = NULL){
  if(distribution == 'gamma'){
    X1.0 = sapply(1:16, function(k) rgamma(n, para[1], para[2]))
  }
  if(distribution == 'beta'){
    X1.0 = sapply(1:16,  function(k) rbeta(n, para[1], para[2]))
  }
  if(distribution == 'uniform'){
    X1.0 = sapply(1:16, function(k) runif(n, para[1], para[2]))
  }
  if(distribution %in% c('gaussian','normal')){
    X1.0 = sapply(1:16, function(k) rnorm(n, para[1], para[2]^2))
  }
  X2.1 = X1.0[,1] * X1.0[,2] + X1.0[,4] * X1.0[,5] + exp(X1.0[,16])
  X2.2 = X1.0[,7] * X1.0[,6] * X1.0[,3] +  X1.0[,8] * X1.0[,9] - X1.0[,15]
  X2.3 = X1.0[,10] * X1.0[,11] * X1.0[,12] * X1.0[,13] + X1.0[,14]^2
  X2.4 = X1.0[,4] * X1.0[,5] * X1.0[,6] + X1.0[,7] * X1.0[,8]
  e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  e2 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  X1 = scale(cbind(X1.0, matrix(0, ncol = p-16, nrow = n)) + e1)
  X2 = scale(cbind(X2.1,X2.2, X2.3,  matrix(0, ncol = p-3, nrow = n)) + e2)
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  return(list(X = list(X1 = X1, X2 = X2)))
}


sim.interaction3 = function(n, p, p.group = 3, q.group = 3, psel = NULL){
  y.latent = q.group*4 
  latent = y.latent + 5
  x = sapply(1:latent, function(i) {
    #rnorm(n = n, mean = 0, sd = 2)})
    runif(n = n, min = 0, max = 1)})
  colnames(x) = paste("x", 1:latent, sep = ".")
  
  ## simulate groups of correlated variables based on functional variables
  
  g = lapply(1:y.latent, function(i) {
    v = sapply(1:p.group, function(j) {
      if(p.group > 1){
        x[, i] + 0.01 + 0.5 * (j - 1) / (p.group - 1)  *
          rnorm(n = n, mean = 0, sd = 0.3)
      }else{
        x[, i] + 0.01 +
          rnorm(n = n, mean = 0, sd = 0.3)
      }
    })
    colnames(v) = paste("X", i, 1:p.group, sep = ".")
    return(v)
  })
  g = do.call(cbind, g)
  
  g.0 = lapply((y.latent+1):latent, function(i) {
    v = sapply(1:p.group, function(j) {
      if(p.group > 1){
        x[, i] + 0.01 + 0.5 * (j - 1) / (p.group - 1)  *
          rnorm(n = n, mean = 0, sd = 0.3)
      }else{
        x[, i] + 0.01 + 
          rnorm(n = n, mean = 0, sd = 0.3)
      }
    })
    colnames(v) = paste("X0", i, 1:p.group, sep = ".")
    return(v)
  })
  g.0 = do.call(cbind, g.0)
  
  if (ncol(g) + ncol(g.0) >= p) {
    warning("No additional independent variables are simulated!")
    ind = NULL
  } else {
    ind = mvrnorm(n = n, mu = rep(0, p - (ncol(g) + ncol(g.0))), 
                  Sigma = diag(.3^2,  p - (ncol(g) + ncol(g.0))))
    colnames(ind) = paste("ind", 1:ncol(ind), sep = ".")
  }
  
  Y = sapply(1:q.group, function(i) {
    v = sapply(1:4, function(j) {
      x[,4*(i - 1) + j]
    })
    y = get_kennel_fn3(v[,1], v[,2], v[,3], v[,4])
    return(y)
  })
  colnames(Y) = sapply(1:q.group, function(i) paste("Y", i, sep = "."))
  
  Y.ind = mvrnorm(n = n, mu = rep(0, p - q.group), 
                  Sigma = diag(.3^2,  p - q.group))
  colnames(Y.ind) = paste("Y.ind", 1:ncol(Y.ind), sep = ".")
  X1 = data.frame(g, g.0, ind)
  X2 = data.frame(Y, Y.ind)
  
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  
  return(list(X = list(X1 = X1, X2 = X2),
              X0 = list(X1 = data.frame(g, g.0, ind), X2 = data.frame(Y, Y.ind))))
}

sim.interaction4 = function(n,p, mu.sd = 2, sigma = 0.2, p.group = 30, q.group = 4, psel = NULL){
  X1.0 = sapply(1:p.group, function(i) rnorm(n, 0, mu.sd))
  X2.0 = sapply(1:(p.group - q.group), function(i) rnorm(n, 0, mu.sd))
  #Y = sapply(1:q.group, function(i) {
  #  v = sapply(1:3, function(j) {
  #    X1.0[,3*(i - 1) + j]
  #  })
  #  y = get_kennel_fn(v[,1], v[,2], v[,3])
  #  return(y)
  #})
  X2.1 = 0.5*exp(X1.0[,1]) + X1.0[,2]^2 - X1.0[,7] * X1.0[,9]
  X2.2 = 3*sin(pi*X1.0[,3]/2) - 2*X1.0[,4]^4 
  X2.3 = abs(X1.0[,5]) + 0.25*X1.0[,6]^3 
  X2.4 = 4/(1 + exp(X1.0[,8]-20)) + X1.0[,10] * X1.0[,11]
  e1 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  e2 = mvrnorm(n, mu = rep(0, p), Sigma = diag(sigma, nrow = p))
  X1 = scale(cbind(X1.0, matrix(0, ncol = p-p.group, nrow = n)) + e1)
  X2 = scale(cbind(X2.1,X2.2,X2.3,X2.4, X2.0, matrix(0, ncol = p-p.group, nrow = n)) + e2)
  colnames(X1) = sapply(1:p, function(i) paste0("X",i))
  colnames(X2) = sapply(1:p, function(i) paste0("Y",i))
  return(list(X = list(X1 = X1, X2 = X2)))
}

sim.illustrate <- function(n, p1, p2){
  
  X1 <- sapply(1:p1, function(i) runif(n, -1, 1))
  X1 <- data.frame(scale(X1))
  colnames(X1) <- paste0("X", 1:p1)
  X2 <- matrix(0, nrow = n, ncol = p2)
  error <- sapply(1:2, function(i) rnorm(n, 0, 1))
  X2[,1] <- X1[,1]^2 + exp(X1[,2])  + error[,1]
  X2[,2] <- X1[,3] - exp(X1[,4]^2) + error[,2]
  X2[,3:p2] <- sapply(1:(p2 - 2), function(i) rnorm(n, 0, 1))
  X2 <- data.frame(scale(X2))
  colnames(X2) <- paste0("Y", 1:p2)
  
  return(list(X = list(X1 = X1, X2 = X2)))
  
}
