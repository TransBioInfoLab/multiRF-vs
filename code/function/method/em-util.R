library(truncnorm)

get_param <- function(x){
  x0 <- x[x == 0]
  x1 <- x[x > 0]
  cut <- mean(x1)
  p0 <- length(x0)/length(x)

  m1 <- mean(x1[x1 < cut])
  sd1 <- sqrt(var(x1[x1 < cut]))
  m2 <- mean(x1[x1 > cut])
  sd2 <- sqrt(var(x1[x1 > cut]))
  list(
    mu = c(m1, m2),
    sigma = c(sd1, sd2),
    p = length(x1[x1 < cut])/length(x1)
  )
}

em <- function(x1, p, mu, sigma, iter, eps, verbose = F, c1 = "truncn", c2 = "gamma"){
  i <- 0
  diff <- 1

  x0 <- x1[x1 == 0]
  x <- x1[x1 > 0]
  p0 <- length(x0)/length(x1)
  post <- cbind(rep(1-p0, length(x1)),
                rep(0, length(x1)),
                rep(p0, length(x1)))

  if(c1 == "normal") {
    Q <- sum.finite(log( p*(1-p0) )+log(dnorm(x, mu[1], sigma[1]))) + sum.finite(log( p*(1-p0) )+log(dnorm(x, mu[2], sigma[2])))
  }
  while(i < iter & diff > eps){

    # E step
    if(c1 == "truncn") {
      comp1 <- p*(1-p0) * dtruncnorm(x, a = 0, b = 1,
                                     mean = mu[1], sd = sigma[1])
    } else {
      comp1 <- p*(1-p0) * dnorm(x, mean = mu[1], sd = sigma[1])
    }


    if(c2 == "gamma") {
      shape <- (mu[2]/sigma[2])^2
      scale <- sigma[2]^2/mu[2]
      comp2 <- (1-p)*(1-p0) * dgamma(x, shape = shape, scale = scale)
    }
    if(c2 == "normal") {
      comp2 <- (1-p)*(1-p0) * dnorm(x, mean = mu[2], sd = sigma[2])
    }

    prop <- comp1/(comp1 + comp2)

    # M step
    p1 <- mean(prop, na.rm = T)

    if(c1  == "truncn") {
      delta_m <- sigma[1] * dnorm(mu[1]/sigma[1])/pnorm(mu[1]/sigma[1])
      mu1 <- c(mean(x*prop*(1-p0))/(p1*(1-p0))-delta_m,
               mean(x*(1-prop)*(1-p0))/((1-p1)*(1-p0)))
      delta_v <- sigma[1]^2-sigma[1]^2*
        (-mu1[1]/sigma[1]*dnorm(mu1[1]/sigma[1])+
           pnorm(mu1[1]/sigma[1]))/pnorm(mu1[1]/sigma[1])

      S1 <- (x - mu1[1])^2
      S2 <- (x - mu1[2])^2
      sigma1 <- c(mean(S1 * prop*(1-p0))/(p1*(1-p0)),mean(S2 * (1-prop)*(1-p0))/((1-p1)*(1-p0)))
      sigma1[1] <- sigma1[1] + delta_v
      sigma1 <- sqrt(sigma1)

      diff <- abs(p1*(1-p0)-p*(1-p0))/5 + abs(mu1[1] - mu[1])/5 +
        abs(sigma1[1] - sigma[1])/5 + abs(mu1[2]-mu[2])/5 +
        abs(sigma1[2]-sigma[2])/5
      i <- i + 1
      p <- p1
      mu <- mu1
      sigma <- sigma1

    } else {

      mu[1] <- sum.finite( prop*(1-p0) * x)/sum.finite( prop*(1-p0))
      mu[2] <- sum.finite((1-prop)*(1-p0) * x)/sum.finite((1-prop)*(1-p0))

      sigma[1] <- sqrt(sum.finite(prop*(1-p0) * (x-mu[1])^2)/sum.finite(prop*(1-p0)))
      sigma[2] <- sqrt(sum.finite((1-prop)*(1-p0) * (x-mu[2])^2)/sum.finite((1-prop)*(1-p0)))

      i <- i + 1
      p <- p1

      diff <- abs(sum(log(comp1 + comp2)) - Q)
      Q <- sum(log(comp1 + comp2))
    }


    if(verbose) {
      cat(i,diff,p,mu,sigma,length(x0)*log(p0)+sum(log(comp1+comp2)),"\n")
    }

  }
  post[x1 > 0,1:2] <- cbind(comp1, comp2)
  list(par=list(p=p*(1-p0),mu=mu,sigma=sigma,post=post),
       conv=list(iter=i,diff=diff))
}

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

find_thres <- function(x1, fdr = 0.05) {

  param <- get_param(x1)

  x <- x1[x1 > 0]

  d <- density(x)
  get_comp0 <- with(d, approxfun(x, y, rule=1))
  comp0 <- get_comp0(x)
  mode0 <- estimate_mode(x)
  idx <- which(x < mode0)
  ss <- c()
  pp <- seq(param$p, 0.99, by = 0.01)

  for(i in pp) {
    mu <- mean(sort(x)[1:ceiling(length(x)*i)])
    sigma <- sd(sort(x)[1:ceiling(length(x)*i)])
    shape <- (mu/sigma)^2
    scale <- sigma^2/mu
    comp2 <- i *  dtruncnorm(x, a = 0, b = 1,
                             mean = mu, sd = sigma)

    ss <- c(ss, sum((comp2[idx] - comp0[idx])^2))

  }

  p <- pp[which.min(ss)]
  mu <- mean(sort(x)[1:ceiling(length(x)*p)])
  sigma <- sd(sort(x)[1:ceiling(length(x)*p)])
  comp2 <- p * dtruncnorm(x, a = 0, b = 1,
                          mean = mu, sd = sigma)
  f <- comp2/comp0
  min(x[f < fdr])
}

#' @keywords internal
sum.finite <- function(x) {
  sum(x[is.finite(x)])
}

em_norm <- function(x, p, mu, sigma, iter, eps){
  i <- 0

  Q <- sum.finite(log(p)+log(dnorm(x, mu[1], sigma[1]))) + sum.finite(log(1-p)+log(dnorm(x, mu[2], sigma[2])))
  diff <- abs(Q)

  while(i < iter & diff > eps){
    # E step
    comp1 <- p * dnorm(x, mean = mu[1], sd = sigma[1])
    comp2 <- (1-p) * dnorm(x, mean = mu[2], sd = sigma[2])
    prop1 <- comp1/(comp1 + comp2)
    prop2 <- comp2/(comp1 + comp2)

    # M step
    p1 <- sum.finite(prop1)/length(x)

    mu[1] <- sum.finite(prop1 * x)/sum.finite(prop1)
    mu[2] <- sum.finite(prop2 * x)/sum.finite(prop2)

    sigma[1] <- sqrt(sum.finite(prop1 * (x-mu[1])^2)/sum.finite(prop1))
    sigma[2] <- sqrt(sum.finite(prop2 * (x-mu[2])^2)/sum.finite(prop2))

    i <- i + 1
    p <- p1

    diff <- abs(sum(log(comp1 + comp2)) - Q)
    Q <- sum(log(comp1 + comp2))
    cat(i,diff,p,mu,sigma,Q,"\n")
  }

  post <- cbind(comp1, comp2)

  list(par=list(p=p,mu=mu,sigma=sigma,post=post),
       conv=list(iter=i,diff=diff))
}
