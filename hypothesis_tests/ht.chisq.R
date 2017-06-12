ht.chisq <- function(data, sigma.sq.0, significance.level = 0.05, 
                     tail = "left", nuisance = NULL) {
  n <- length(data)
  s2 <- var(data)
  t <- s2/sigma.sq.0*(n-1)
  
  
  if(tail == "left") {
    if(t < qchisq(significance.level, n-1)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else if(tail == "right") {
    if(t > qchisq(1 - significance.level, n-1)) {
      return(TRUE)
    } else {
      return(FALSE)
    }  
  } else {
    if(abs(t) > qchisq(1 - significance.level/2, n-1)) {
      return(TRUE)
    } else {
      return(FALSE)
    }    
  }
}




verify.chisq <- function(n = 100, rep = 1000, sig.level = 0.05, tail = "left") {
  
  delta.set <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
  
  ht.list <- "ht.saddle.chisq"
  df.set <- c(1, 2, 3, 4, 5, 10)
  
  rdist <- rchisq
  variance.fn <- function(df) { 2*df }
  
  plots <- vector("list", length(df.set))
  i <- 1
  for(df in df.set) {
    coverage.data <- data.frame(Test = NULL, Delta = NULL, Coverage = NULL)
    for(delta in delta.set) {
      rdist.wpar <- function(n) { rdist(n, df) }
      variance <- variance.fn(df)
      coverage <- sim.dist(rdist.wpar, variance, delta, n, rep, ht.list, sig.level, tail)
      coverage.data <- rbind(coverage.data, 
                             data.frame(Test = ht.list, 
                                        Delta = rep(delta, length(ht.list)), 
                                        Coverage = coverage))
    }
    distribution <- "chi-sq"
    parameter.names <- "df"
    parameter.vals <- df
    print(paste(distribution, ": df =", df))
    print(get.table(coverage.data))
    plots[[i]] <- get.coverage.plot(coverage.data, significance.level, 
                                    parameter.names, parameter.vals)
    i <- i + 1
  }
  get.plot(plots, layout = rbind(1:3, 4:6, rep(7, 3)))
}