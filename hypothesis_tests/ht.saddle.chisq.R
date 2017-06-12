#library("distr")

ht.saddle.chisq <- function(obs, sigma.sq.0, significance.level = 0.05,
                            tail = "left", nuisance = NULL) {
  n <- length(obs)
  
  if(var(obs) > sigma.sq.0) {
    return(FALSE)
  }
  
  #Exponential Family setup.
  t <- function(x) { 0.5*log(x) }
  t.bar <- mean(t(obs))
  k <- function(theta) { theta*log(2)/2 + lgamma(theta/2)}
  kp <- function(theta) { log(2)/2 + digamma(theta/2)/2 }
  kpp <- function(theta) { trigamma(theta/2)/4 }
  theta.hat <- 2*igamma(2*t.bar - log(2))
  theta.0 <- sigma.sq.0/2
  
  temp <- 2*n*((theta.hat - theta.0)*t.bar - k(theta.hat) + k(theta.0))
  if(temp < 0) {
    #print("Temp is negative. ht.chisq")
    return(FALSE)
  }
  
  ts.R <- sign(theta.hat - theta.0)*sqrt(temp)  
  ts.U <- sqrt(n)*(theta.hat - theta.0)*sqrt(kpp(theta.hat))
  ts <- ts.R + 1/ts.R*log(ts.U/ts.R)
  
  
  if(tail == "left") {
    if(ts <= qnorm(significance.level)) {
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



verify.chisq <- function(n = 100, rep = 1000, sig.level = 0.05, tail = "left") {

  delta.set <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
  
  ht.list <- "ht.saddle.chisq"
  df.set <- c(1, 2, 3, 4, 5, 10)
  
  rdist <- rchisq
  variance.fn <- function(df) { 2*df }
  
  distribution <- "chisq"
  coverage.data <- NULL
  for(df in df.set) {
    for(delta in delta.set) {
      rdist.wpar <- function(n) { rdist(n, df) }
      variance <- variance.fn(df)
      coverage <- sim.dist(rdist.wpar, variance, delta, n, rep, ht.list, sig.level, tail)
      coverage.data <- rbind(coverage.data, 
                             data.frame(distribution = distribution, 
                                        test = ht.list, 
                                        delta = delta, 
                                        coverage = coverage,
                                        df = df))
    }
  }
  return(coverage.data)
}