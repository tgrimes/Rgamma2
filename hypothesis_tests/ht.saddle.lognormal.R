#require("MASS")

ht.saddle.lognormal <- function(obs, sigma.sq.0, significance.level = 0.05,
                              tail = "left", nuisance = NULL) {
  n <- length(obs)
  
  if(var(obs) > sigma.sq.0) {
    return(FALSE)
  }
  
  lsigma.sq <- nuisance
  if(is.null(nuisance)) {
    lsigma.sq <- tryCatch({
      lsigma.sq <- ((fitdistr(obs, "lognormal"))$estimate[[2]])^2
    }, error = function(e) { 
      print("fitdist failed, returning NULL.")
      return(NULL)
    })
    if(is.null(lsigma.sq)) {
      print("fitdist failed. returning NaN. ht.lognormal.")
      return(NaN)
    }
  }
  
  #Exponential Family setup.
  t <- function(x) { log(x)/lsigma.sq }
  t.bar <- mean(t(obs))
  k <- function(theta) { theta^2/(2 * lsigma.sq) }
  kp <- function(theta) { theta/lsigma.sq }
  kpp <- function(theta) { 1/lsigma.sq }
  theta.hat <- lsigma.sq * t.bar
  theta.0 <- (log(sigma.sq.0/(exp(2 * lsigma.sq) - exp(lsigma.sq))))/2
  
  temp <- 2*n*((theta.hat - theta.0)*t.bar - k(theta.hat) + k(theta.0))
  if(temp < 0) {
    print("Temp is negative. ht.lognormal")
    return(FALSE)
  }
  
  ts.R <- sign(theta.hat - theta.0)*sqrt(temp)
  ts.U <- sqrt(n)*(theta.hat - theta.0)*sqrt(kpp(theta.hat))
  ts <- ts.R + 1/ts.R*log(ts.U/ts.R)
  
  if(is.nan(ts)) {
    #Test procedure failed.
    print("ts is NaN. ht.lognormal")
    return(NaN)
  }
  
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



verify.lognormal <- function(n = 100, rep = 1000, sig.level = 0.05, tail = "left") {
  
  delta.set <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
  
  ht.list <- "ht.saddle.lognormal"
  lmu.set <- c(0)
  lsigma.set <- c(0.669, 0.595, 0.576, 0.554, 0.530, 
                  0.503,  0.471,  0.433,  0.240)

  rdist <- rlnorm
  variance.fn <- function(lmu, lsigma) { (exp(lsigma^2) - 1) * exp(2 * lmu + lsigma^2) }
  
  distribution <- "lognormal"
  coverage.data <- NULL
  for(lmu in lmu.set) {
    for(lsigma in lsigma.set) {
      for(delta in delta.set) {
        rdist.wpar <- function(n) { rdist(n, lmu, lsigma) }
        variance <- variance.fn(lmu, lsigma)
        coverage <- sim.dist(rdist.wpar, variance, delta, n, rep, ht.list, sig.level, tail, lsigma^2)
        coverage.data <- rbind(coverage.data, 
                               data.frame(distribution = distribution, 
                                          test = ht.list, 
                                          delta = delta, 
                                          coverage = coverage,
                                          lmu = lmu,
                                          lsigma = lsigma))
      }
    }
  }
  return(coverage.data)
}
