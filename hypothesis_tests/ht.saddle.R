ht.saddle <- function(obs, sigma.sq.0, significance.level = 0.05,
                      tail = "left", nuisance = NULL) {
  
  #distributions <- c("normal", "gamma", "exponential", "chi-squared")
  tests <- c("ht.saddle.normal", "ht.saddle.gamma2")
  
  #Fit normal.
  fits <- tryCatch({
              (fitdistr(obs, "normal"))$loglik
            }, error = function(e) {
              return(-Inf)
            })
  
  #Fit gamma.
  fits <- c(fits, 
            tryCatch({
              (fitdistr(obs, "gamma"))$loglik
            }, error = function(e) {
              return(-Inf)
            }))
  
  test = tests[fits == max(fits)]
  if(test == tests[1]) {
    return((get(test))(obs, sigma.sq.0,  significance.level/2, tail, nuisance))
  } else {
    #Do not divide significance level by 2; the adjusted gamma test does this.
    return((get(test))(obs, sigma.sq.0,  significance.level, tail, nuisance))
  }
}

