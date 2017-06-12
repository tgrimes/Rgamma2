#require("MASS")

ht.saddle.weibull <- function(obs, sigma.sq.0, significance.level = 0.05,
                      tail = "left", nuisance = NULL) {
  n <- length(obs)
  
  if(var(obs) > sigma.sq.0) {
    return(FALSE)
  }
  
  shape <- nuisance
  if(is.null(nuisance)) {
    shape <- tryCatch({
        shape <- (fitdistr(obs, "weibull"))$estimate[[1]]
    }, error = function(e) {
        #print("fitdist failed, returning NULL.")
        return(NULL)
    })
    if(is.null(shape)) {
      #print("fitdist failed. returning NaN. ht.weibull.")
      return(NaN)
    }
  }
  
  #Exponential Family setup.
  t <- function(x) { -(x^shape) }
  t.bar <- mean(t(obs))
  k <- function(theta) { -log(theta) }
  kp <- function(theta) { -1/theta }
  kpp <- function(theta) { 1/theta^2 }
  theta.hat <- -1/t.bar
  theta.0 <- ((gamma(1 + 2/shape) - (gamma(1 + 1/shape))^2)/sigma.sq.0)^(shape/2)

  temp <- 2*n*((theta.hat - theta.0)*t.bar - k(theta.hat) + k(theta.0))
  if(is.nan(theta.0) || theta.0 < 0) {
    print("theta.0 is nan or < 0. ht.weibull")
    return(FALSE)
  }
  if(is.nan(theta.hat) || theta.hat < 0) {
    print("theta.hat is nan or < 0. ht.weibull")
    return(FALSE)
  }
  if(is.nan(temp)) {
    #print("Temp is nan. ht.weibull")
    return(FALSE)
  }
  if(!is.numeric(temp)) {
    print("Temp is not numeric. ht.weibull")
    return(FALSE)
  }
  if(temp < 0) {
    #print("Temp is negative. ht.weibull")
    return(FALSE)
  }
  
  ts.R <- sign(theta.hat - theta.0)*sqrt(temp)
  ts.U <- sqrt(n)*(theta.hat - theta.0)*sqrt(kpp(theta.hat))
  ts <- ts.R + 1/ts.R*log(ts.U/ts.R)
  
  
  if(is.nan(ts)) {
    #Test procedure failed.
    #print("ts is NaN. ht.weibull")
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




verify.weibull <- function(n = 100, rep = 1000, sig.level = 0.05, tail = "left") {

  delta.set <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
  
  ht.list <- c("ht.saddle.weibull")
  shape.set <- c(0.764, 0.859, 0.886, 0.917, 0.955, 
                 1.000,  1.056,  1.128,  1.615)
  scale.set <- c(1)
  
  rdist <- rweibull
  variance.fn <- function(shape, scale) { scale^2 * (gamma(1 + 2/shape) - (gamma(1 + 1/shape)^2)) }

  distribution <- "weibull"
  coverage.data <- NULL
  for(shape in shape.set) {
    for(scale in scale.set) {
      for(delta in delta.set) {
        rdist.wpar <- function(n) { rdist(n, shape, scale) }
        variance <- variance.fn(shape, scale)
        coverage <- sim.dist(rdist.wpar, variance, delta, n, rep, ht.list, sig.level, tail, shape)
        coverage.data <- rbind(coverage.data, 
                               data.frame(distribution = distribution, 
                                          test = ht.list, 
                                          delta = delta, 
                                          coverage = coverage,
                                          shape = shape,
                                          scale = scale))
      }
    }
  }
  return(coverage.data)
}