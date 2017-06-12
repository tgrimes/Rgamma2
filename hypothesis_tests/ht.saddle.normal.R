ht.saddle.normal <- function(obs, sigma.sq.0, significance.level = 0.05,
                             tail = "left", nuisance = NULL) {
  n <- length(obs)
  
  if(var(obs) > sigma.sq.0) {
    return(FALSE)
  }
  
  mu <- nuisance
  if(is.null(nuisance)) {
    mu <- mean(obs)
  }
  
  s.sq <- var(obs)
  
  #Exponential Family setup:
  t <- function(x) { -0.5*(x - mu)^2 + 0.5*mu^2 }
  t.bar <- 0.5*(mu^2 - s.sq)
  k <- function(theta) { 0.5*mu^2*theta - 0.5*log(theta) }
  kp <- function(theta) { 0.5*mu^2 - 0.5/theta }
  kpp <- function(theta) { 0.5/theta^2 }
  theta.hat <- 1/s.sq
  theta.0 <- 1/sigma.sq.0
  
  temp <- 2*n*((theta.hat - theta.0)*t.bar - k(theta.hat) + k(theta.0))
  if(temp < 0) {
    print("Temp is negative. ht.normal")
    return(NaN)
  }
  
  ts.R <- sign(theta.hat - theta.0)*sqrt(temp)
  #ts.U <- sqrt(n)*(theta.hat - theta.0)*sqrt(kpp(theta.hat))
  ts.U <- sqrt(n/2)*(1-theta.0/theta.hat)
  ts <- ts.R + 1/ts.R*log(ts.U/ts.R)
  
  if(is.nan(ts)) {
    #Test procedure failed.
    print("ts is NaN. ht.normal")
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
    if(ts < qnorm(significance.level)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else {
    if(abs(ts) >= qnorm(1 - significance.level/2)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}


verify.normal <- function(n = 100, rep = 1000, sig.level = 0.05, tail = "left") {
  
  delta.set <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
  
  ht.list <- c("ht.saddle.normal")
  mu.set <- c(0.5, 3, 8)
  sigma.set <- c(0.5, 3, 8)
  
  rdist <- function(n, mu, sigma) { rnorm(n, mu + 10*sigma, sigma) }
  variance.fn <- function(mu, sigma) { sigma^2 }
  
  distribution <- "normal"
  coverage.data <- NULL
  for(mu in mu.set) {
    for(sigma in sigma.set) {
      for(delta in delta.set) {
        rdist.wpar <- function(n) { rdist(n, mu, sigma) }
        variance <- variance.fn(mu, sigma)
        coverage <- sim.dist(rdist.wpar, variance, delta, n, rep, ht.list, sig.level, tail, mu)
        coverage.data <- rbind(coverage.data, 
                               data.frame(distribution = distribution, 
                                          test = ht.list, 
                                          delta = delta, 
                                          coverage = coverage,
                                          mu = mu,
                                          sigma = sigma))
      }
    }
  }
  return(coverage.data)
}
