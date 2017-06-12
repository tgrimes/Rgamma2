#require("moments")

ht.z6 <- function(data, sigma.sq.0, significance.level = 0.05, 
                  tail = "left", nuisance = NULL) {
  #Center data:
  data <- data - mean(data)
  
  k <- (all.cumulants(all.moments(data, 6)))[2:7]
  s2 <- var(data)
  s4 <- s2^2
  n <- length(data)
  if(k[4] + 2 * s4 < 0) {
    k[4] <- 0
  } else if(k[4]/(n*s2) + 2*sigma.sq.0/(n - 1) < 0) {
    k[4] <- 0
  }
  b1.hat <- -sqrt(s4/(k[4] + 2*s4))
  b2.hat <- (k[6] + 12*k[4]*s2 + 4*k[3]^2 + 8*s2^3)/(k[4] + 2*s4)^(3/2)
  
  #if(is.nan(b1.hat) || is.nan(b2.hat)) {
  #  return(FALSE)
  #}
  
  temp <- (k[4]*sigma.sq.0)/(n*s2) + (2*sigma.sq.0^2)/(n-1)
  if(temp < 0) {
    #print("temp < 0 in ht.z6")
    return(NaN)
  }
  
  #Tau = sqrt(E(X-mu)^4 - sigma^4); Tau/sqrt(n) is estimated as in z6.
  z6 <- (s2 - sigma.sq.0)/(sqrt(temp))
  
  if(tail == "left") {
    z.alpha <- qnorm(significance.level) 
    right <- z.alpha + sqrt(1/n)*(b1.hat + b2.hat*((z.alpha)^2 - 1)/6)
    if(is.nan(z6) || is.nan(right)) {
      print("z6 or right is nan in z6.")
      return(NaN)
    }
    
    if(z6 < right) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else if(tail == "right") {
    z.alpha <- qnorm(1 - significance.level)
    right <- z.alpha + sqrt(1/n)*(b1.hat + b2.hat*((z.alpha)^2 - 1)/6)
    if(is.nan(z6) || is.nan(right)) {
      print("z6 or right is nan in z6.")
      return(NaN)
    }
    
    if(z6 > right) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else {
    z.alpha <- qnorm(1 - significance.level/2) 
    right <- z.alpha + sqrt(1/n)*(b1.hat + b2.hat*((z.alpha)^2 - 1)/6)
    if(is.nan(z6) || is.nan(right)) {
      print("z6 or right is nan in z6.")
      return(NaN)
    }
    
    if(abs(z6) > right) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}