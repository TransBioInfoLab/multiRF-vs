##### Simulation Model #####
library(MASS)
library(corpcor)

######## nonlinear function #########

get_kennel_fn4 = function(x1, x2){
  y = 0.25 * exp( 4 * x1) + 4/(1+exp(-20*(x2 - 0.5))) + rnorm(n = length(x1), mean = 0, sd = 0.2)
  y
}

sim.nonlinear2 = function(n, p, j = 2, mu.sd = 2, rho = 0, sigma = 0.3, psel = 20, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  
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

sim.nonlinear3 = function(n, p, p.group = 3, q.group = 3, p.b = 0, psel = NULL, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
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



######## interaction function #########

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

simulate_multi <- function(n               = 200,
                           p               = 300,
                           beta            = 2,
                           num_signal_responses  = 5,
                           num_noise_responses   = 20,
                           noise_sd        = 0.2,
                           seed            = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Generate feature matrix X (n × p)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  
  # 2. Prepare response matrix
  total_responses <- num_signal_responses + num_noise_responses
  Y <- matrix(0, nrow = n, ncol = total_responses)
  colnames(Y) <- paste0("Y", seq_len(total_responses))
  
  for (j in seq_len(num_signal_responses)) {
    # indices for this triplet
    idx <- ((j - 1) * 3 + 1):((j - 1) * 3 + 3)
    x1 <- X[, idx[1]]
    x2 <- X[, idx[2]]
    x3 <- X[, idx[3]]
    
    # nonlinear interaction
    eta_j <- (x1^2 - 1) * (x2 + 0.5) * (x3 - .5)
    
    # continuous response with noise
    Y[, j] <- beta * eta_j + rnorm(n, mean = 0, sd = noise_sd)
  }
  # 4. Create the pure-noise Y's
  if (num_noise_responses > 0) {
    noise_cols <- (num_signal_responses + 1):total_responses
    for (j in noise_cols) {
      Y[, j] <- rnorm(n, mean = 0, sd = noise_sd)
    }
  }
  
  
  # 4. Return as data.frame: responses first, then predictors
  list(X = list(X1 = scale(X), X2 = scale(Y)))
}

simulate_multivar <- function(n         = 200,
                              p         = 300,
                              flip_prop = 0.10,
                              seed      = NULL,
                              binary = T,
                              noise_sd = 0.1) {
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Generate predictor matrix X (n × p) ~ N(0,1)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  
  # 2. Compute latent nonlinear signals for Y1 and Y2
  #    Shared triple: X1, X2, X3
  x1 <- X[, "X1"]
  x2 <- X[, "X2"]
  x3 <- X[, "X3"]
  
  #    Unique features:
  #      Y1 has X4, X5
  #      Y2 has X6, X7
  x4 <- X[, "X4"]
  x5 <- X[, "X5"]
  x6 <- X[, "X6"]
  x7 <- X[, "X7"]
  
  #    Define η1 and η2
  eta1 <- 1.5 * (x1^2 - 1) * (x2 + 0.5) * (x3 - .5) + x4 + x5
  eta2 <- 2.0 * (x1^2 - 0.5) * (x2 + 0.5) * (x3 - .5) + x6 - x7
  
  if(binary) {
    # 3. Generate binary outcomes via logistic model
    p1 <- 1 / (1 + exp(-eta1))
    p2 <- 1 / (1 + exp(-eta2))
    Y1 <- rbinom(n, size = 1, prob = p1)
    Y2 <- rbinom(n, size = 1, prob = p2)
    
    # 4. Optional label noise: flip a fraction of each outcome
    if (flip_prop > 0) {
      n_flip1 <- floor(flip_prop * n)
      n_flip2 <- floor(flip_prop * n)
      flip_idx1 <- sample.int(n, size = n_flip1)
      flip_idx2 <- sample.int(n, size = n_flip2)
      Y1[flip_idx1] <- 1 - Y1[flip_idx1]
      Y2[flip_idx2] <- 1 - Y2[flip_idx2]
    }
    
  } else {
    Y1 <- eta1 + rnorm(n, mean = 0, sd = noise_sd)
    Y2 <- eta2 + rnorm(n, mean = 0, sd = noise_sd)
  }
 
  # 5. Return as a list
  list(
    X  = list(X1 = as.data.frame(X),
    X2  = data.frame(Y1 = Y1, Y2 = Y2))
  )
}
