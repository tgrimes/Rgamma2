#library("moments")
#library("MASS")

ht.saddle.gamma2 <- function(obs, sigma.sq.0, significance.level = 0.05, 
                             tail = "left", nuisance = NULL) {
  n <- length(obs)
  
  alpha.hat <- nuisance
  if(is.null(nuisance)) {
    #Approximate maximum likelihood estimator.
    log.xbar <- log(mean(obs))
    logx.bar <- mean(log(obs))
    s <- log.xbar - logx.bar 
    alpha.hat <- (3 - s + sqrt((s-3)^2 + 24*s))/(12*s)
  } 
  
  obs <- obs/alpha.hat
  sigma.sq.0 <- sigma.sq.0/alpha.hat
  
  alpha.hat <- 1
  
  #Exponential Family setup.
  t <- function(x) { -x }
  t.bar <- mean(t(obs))
  k <- function(theta) { -alpha.hat*log(theta) }
  kp <- function(theta) { -alpha.hat/theta }
  kpp <- function(theta) { alpha.hat/theta^2 }
  theta.hat <- -alpha.hat/t.bar
  theta.0 <- sqrt(alpha.hat/sigma.sq.0)
  
  temp <- 2*n*((theta.hat - theta.0)*t.bar - k(theta.hat) + k(theta.0))
  if(temp < 0) {
    print("Temp is negative. ht.gamma2")
    return(NaN)
  }
  
  ts.R <- sign(theta.hat - theta.0)*sqrt(temp)
  ts.U <- sqrt(n)*(theta.hat - theta.0)*sqrt(kpp(theta.hat))
  ts <- ts.R + 1/ts.R*log(ts.U/ts.R)
  
  
  if(is.nan(ts)) {
    #Test procedure failed.
    print("ts is NaN. ht.gamma2")
    return(NaN)
  }
  
  if(tail == "left") {
    if(ts >= qt(1 - significance.level/2, n)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else if(tail == "right") {
    if(ts < qt(significance.level/2, n)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else {
    if(abs(ts) >= qt(1 - significance.level/4, n)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}


ht.saddle.gamma <- function(obs, sigma.sq.0, significance.level = 0.05,
                            tail = "left", nuisance = NULL) {
  n <- length(obs)
  
  alpha.hat <- nuisance
  if(is.null(nuisance)) {
    #Approximate maximum likelihood estimator.
    log.xbar <- log(mean(obs))
    logx.bar <- mean(log(obs))
    s <- log.xbar - logx.bar #Previously used (1/2)(1/(log.xbar - logx.bar))?
    alpha.hat <- (3 - s + sqrt((s-3)^2 + 24*s))/(12*s)
  } 
  
  #Exponential Family setup.
  t <- function(x) { -x }
  t.bar <- mean(t(obs))
  k <- function(theta) { -alpha.hat*log(theta) }
  kp <- function(theta) { -alpha.hat/theta }
  kpp <- function(theta) { alpha.hat/theta^2 }
  theta.hat <- -alpha.hat/t.bar
  theta.0 <- sqrt(alpha.hat/sigma.sq.0)
  
  temp <- 2*n*((theta.hat - theta.0)*t.bar - k(theta.hat) + k(theta.0))
  if(temp < 0) {
    print("Temp is negative. ht.gamma")
    return(FALSE)
  }
  
  ts.R <- sign(theta.hat - theta.0)*sqrt(temp)
  ts.U <- sqrt(n)*(theta.hat - theta.0)*sqrt(kpp(theta.hat))
  ts <- ts.R + 1/ts.R*log(ts.U/ts.R)
  
  if(is.nan(ts)) {
    #Test procedure failed.
    print("ts is nan. ht.gamma")
    return(NaN)
  }
  
  if(tail == "left") {
    if(ts >= qnorm(1 - significance.level)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else if(tail == "right") {
    return(NULL)
  }
  else {
    return(NULL)
  }
}



verify.gamma <- function(n = 100, rep = 1000, sig.level = 0.05, tail = "left") {
  
  delta.set <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
  
  ht.list <- c("ht.saddle.gamma")
  alpha.set <- c(0.400, 0.600, 0.667, 0.750, 0.857, 
                 1.000,  1.200,  1.500,  6.000)
  beta.set <- c(1)
  
  rdist <- rgamma
  variance.fn <- function(alpha, beta) { alpha/beta^2 }
  
  distribution <- "gamma"
  coverage.data <- NULL
  for(alpha in alpha.set) {
    for(beta in beta.set) {
      for(delta in delta.set) {
        rdist.wpar <- function(n) { rdist(n, alpha, beta) }
        variance <- variance.fn(alpha, beta)
        coverage <- sim.dist(rdist.wpar, variance, delta, n, rep, ht.list, sig.level, tail, alpha)
        coverage.data <- rbind(coverage.data, 
                               data.frame(distribution = distribution, 
                                          test = ht.list, 
                                          delta = delta, 
                                          coverage = coverage,
                                          alpha = alpha,
                                          beta = beta))
      }
    }
  }
  return(coverage.data)
}
