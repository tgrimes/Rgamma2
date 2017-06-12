ht.robust <- function(data, sigma.sq.0, significance.level = 0.05, 
                      tail = "left", nuisance = NULL) {
  eta.hat <- kurtosis(data)
  d.hat <- 1/(1 + eta.hat/2)
  n <- length(data)
  s2 <- var(data)
  t <- s2*(n-1)*d.hat/sigma.sq.0
  
  if(tail == "left") {
    if(t < qchisq(significance.level, ceiling((n-1)*d.hat))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else if(tail == "right") {
    if(t > qchisq(1 - significance.level, ceiling((n-1)*d.hat))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else {
    if(abs(t) > qchisq(1 - significance.level/2, ceiling((n-1)*d.hat))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}