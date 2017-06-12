ht.saddle.exp <- function(obs, sigma.sq.0, significance.level = 0.05,
                                tail = "left", nuisance = NULL) {
  n <- length(obs)
  
  if(var(obs) > sigma.sq.0) {
    return(FALSE)
  }
  
  #Exponential Family setup.
  t <- function(x) { -x }
  t.bar <- mean(t(obs))
  k <- function(theta) { -log(theta) }
  kp <- function(theta) { -1/theta }
  kpp <- function(theta) { 1/theta^2 }
  theta.hat <- -1/t.bar
  theta.0 <- sqrt(1/sigma.sq.0)
  
  temp <- 2*n*((theta.hat - theta.0)*t.bar - k(theta.hat) + k(theta.0))
  if(temp < 0) {
    print("Temp is negative. ht.exp")
    return(FALSE)
  }
  
  ts.R <- sign(theta.hat - theta.0)*sqrt(temp)
  ts.U <- sqrt(n)*(theta.hat - theta.0)*sqrt(kpp(theta.hat))
  ts <- ts.R + 1/ts.R*log(ts.U/ts.R)
  
  
  if(is.nan(ts)) {
    #Test procedure failed.
    print("ts is NaN. ht.exp")
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


verify.exp <- function(n = 100, rep = 1000, sig.level = 0.05, tail = "left") {

  delta.set <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)

  ht.list <- "ht.saddle.exp"
  lambda.set <- c(0.2, 0.05, 1, 1.5, 2, 5, 8, 20, 50)

  rdist <- rexp
  variance.fn <- function(lambda) { 1/lambda^2 }
  
  distribution <- "exponential"
  coverage.data <- NULL
  for(lambda in lambda.set) {
    for(delta in delta.set) {
      rdist.wpar <- function(n) { rdist(n, lambda) }
      variance <- variance.fn(lambda)
      coverage <- sim.dist(rdist.wpar, variance, delta, n, rep, ht.list, sig.level, tail)
      coverage.data <- rbind(coverage.data, 
                             data.frame(distribution = distribution, 
                                        test = ht.list, 
                                        delta = delta, 
                                        coverage = coverage,
                                        lambda = lambda))
    }
  }
  return(coverage.data)
}