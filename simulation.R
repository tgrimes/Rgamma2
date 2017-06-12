run.simulation <- function(n = 20, m = 1000, sig.level = 0.05, ht.list = NULL,
                           delta.set = c(1, 2, 3, 4), tail = "left", DEFAULT.DISTS = TRUE, DIST = NULL,
                           rdist = NULL, param = NULL, variance.fn = NULL, header = NULL) {
  
  if(is.null(ht.list)) {
    ht.list <- c("ht.chisq", "ht.robust", "ht.saddle.normal", "ht.z6", "ht.saddle3", "ht.saddle.gamma3")
  }
  
  results <- NULL

  if(!is.null(rdist)) {
    coverage <- simulate(rdist, param, variance.fn, delta.set, n, m, ht.list, sig.level,
                         tail, header)
    results <- c(results, list(coverage))
  }
  
  # NOTE: If any distributions are added, `table.coverage.R` will need to be updated.
  
  if(DEFAULT.DISTS == TRUE || any(strcmpi("normal", DIST))) {
    #Normal distribution.
    print("Normal:")
    param <- list(mu.set = c(0, 10, 100),
                  sigma.set = c(1, 10, 100))    
    header <- list("Normal (mu, sigma)")
    rdist <- function(n, mu, sigma) { rnorm(n, mu, sigma) + 1000 } # Offset for positive values.
    variance.fn <- function(mu, sigma) { sigma^2 }
    coverage <- simulate(rdist, param, variance.fn, delta.set, n, m, ht.list, sig.level,
                         tail, header)
    results <- c(results, list(coverage))
  }
  
  if(DEFAULT.DISTS == TRUE || any(strcmpi("t", DIST))) {
    #t distribution.
    print("t:")
    param <- list(df.set = c(4.4, 4.6, 4.667, 4.750, 4.857, 5, 5.2, 5.5, 10))
    header <- list("t (df)")
    rdist <- function(n, df) { rt(n, df) + 1000 } # Offset for positive values.
    variance.fn <- function(df) { df/(df - 2) }
    coverage <- simulate(rdist, param, variance.fn, delta.set, n, m, ht.list, sig.level,
                         tail, header)
    results <- c(results, list(coverage))
  }
  
  if(DEFAULT.DISTS == TRUE || any(strcmpi("chisq", DIST))) {
    #Chi-sq distribution.
    print("Chi-sq:")
    param <- list(df.set = c(0.8, 1.2, 1.333, 1.5, 1.714, 2, 2.4, 3, 12))
    header <- list("Chi-sq (df)")
    rdist <- rchisq
    variance.fn <- function(df) { 2*df }
    coverage <- simulate(rdist, param, variance.fn, delta.set, n, m, ht.list, sig.level,
                         tail, header)
    results <- c(results, list(coverage))
  }
  
  if(DEFAULT.DISTS == TRUE || any(strcmpi("gamma", DIST))) {
    #Gamma distribution.
    print("Gamma:")
    param <- list(shape.set = c(0.400, 0.600, 0.667, 0.750, 0.857, 
                                1.000,  1.200,  1.500,  6.000),
                  scale.set = 1)
    header <- list("Gamma (shape, scale)")
    rdist <- function(x, shape, scale) rgamma(x, shape, scale = scale)
    variance.fn <- function(shape, scale) { shape * scale^2 }
    coverage <- simulate(rdist, param, variance.fn, delta.set, n, m, ht.list, sig.level,
                         tail, header)
    results <- c(results, list(coverage))
  }
  
  
  if(DEFAULT.DISTS == TRUE || any(strcmpi("weibull", DIST))) {
    #Weibull distribution.
    print("Weibull:")
    # param <- list(shape.set = c(0.5, 1, 5, 10, 100),
    #               scale.set = c(0.5, 10))
    param <- list(shape.set = c(0.764, 0.859, 0.886, 0.917, 0.955, 
                                1.000,  1.056,  1.128,  1.615),
                  scale.set = 1)
    header <- list("Weibull (shape, scale)")
    rdist <- rweibull
    variance.fn <- function(shape, scale) { scale^2 * (gamma(1 + 2/shape) - (gamma(1 + 1/shape)^2)) }
    coverage <- simulate(rdist, param, variance.fn, delta.set, n, m, ht.list, sig.level,
                         tail, header)
    results <- c(results, list(coverage))
  }
  
  if(DEFAULT.DISTS == TRUE || any(strcmpi("lnorm", DIST))) {
    #Lognormal distribution.
    print("Log-Normal:")
    param <- list(lmu.set = 0,
                  lsigma.set = c(0.669, 0.595, 0.576, 0.554, 0.530, 
                                 0.503,  0.471,  0.433,  0.240))
    header <- list("Lognormal (lmu, lsigma)")
    rdist <- rlnorm
    variance.fn <- function(lmu, lsigma) { (exp(lsigma^2) - 1) * exp(2 * lmu + lsigma^2) }
    coverage <- simulate(rdist, param, variance.fn, delta.set, n, m, ht.list, sig.level,
                         tail, header)
    results <- c(results, list(coverage))
  }
  
  if(DEFAULT.DISTS == TRUE || any(strcmpi("inverse-gamma", DIST))) {
    #Inverse-gamma distribution.
    print("Inv-Gamma:")
    rdist <- function(x, shape, scale) rinvgamma(x, shape, scale = scale)
    variance.fn <- function(shape, scale) { scale^(2)/((shape - 1)^2*(shape - 2)) }
    param <- list(shape = c(6.462, 7.530, 7.880, 8.314, 8.870, 
                            9.606, 10.629, 12.155, 34.756),
                  scale = 1)
    header <- list("Inverse-gamma (shape, scale)")
    coverage <- simulate(rdist, param, variance.fn, delta.set, n, m, ht.list, sig.level,
                         tail, header)
    results <- c(results, list(coverage))
  }
  
  return(results)
}



simulate <- function(rdist, param, variance.fn, delta.set, n, m, ht.list, 
                     sig.level = 0.05, tail = "left", header = NULL) {
  
  num.param <- length(param)
  coverage.list <- header
  if(num.param == 1) {
    for(param1 in param[[1]]) {
      coverage.data <- NULL
      for(delta in delta.set) {
        rdist.wpar <- function(n) { rdist(n, param1) }
        variance <- variance.fn(param1)
        coverage <- sim.dist(rdist.wpar, variance, delta, n, m, ht.list,
                             sig.level, tail)
        coverage.data <- rbind(coverage.data, 
                               data.frame(Test = ht.list, 
                                          Delta = rep(delta, length(ht.list)), 
                                          Coverage = coverage,
                                          param1 = rep(param1, length(ht.list))))
      }
      coverage.list <- c(coverage.list, list(coverage.data))
      
    }
  } else if(num.param == 2) {
    for(param1 in param[[1]]) {
      for(param2 in param[[2]]) {
        coverage.data <- NULL
        for(delta in delta.set) {
          rdist.wpar <- function(n) { rdist(n, param1, param2) }
          variance <- variance.fn(param1, param2)
          coverage <- sim.dist(rdist.wpar, variance, delta, n, m, ht.list,
                               sig.level, tail)
          coverage.data <- rbind(coverage.data, 
                                 data.frame(Test = ht.list, 
                                            Delta = rep(delta, length(ht.list)), 
                                            Coverage = coverage,
                                            param1 = rep(param1, length(ht.list)),
                                            param2 = rep(param2, length(ht.list))))
        }
        coverage.list <- c(coverage.list, list(coverage.data))
        
      }
    }
  } else {
    return(NA)
  }
  
  return(coverage.list)
}



sim.dist <- function(rdist.wpar, variance, delta, n, m, ht.list, 
                     sig.level = 0.05, tail = "left", nuisance = NULL) {
  
  x <- matrix(rdist.wpar(n*m), nrow = n, ncol = m)
  
  if(tail == "right") {
    delta <- 1/delta
  }
  sigma.sq.0 <- delta * variance
  
  results <- NULL
  for(i in 1:length(ht.list)) {
    results <- cbind(results, apply(x, 2, get(ht.list[i]), 
                                    sigma.sq.0 = sigma.sq.0,
                                    significance.level = sig.level, 
                                    tail = tail,
                                    nuisance))
  }
  
  #Account for any samples that resulted in a NaN.
  print(apply(results, 2, function(x) sum(is.nan(x))))
  coverage <- apply(results, 2, function(x) mean(x[!is.nan(x)]))
  return(coverage)
}


