table.single.dist <- function(data, params) {
  Tests <- unique(data[, 1])
  Deltas <- unique(data[, 2])
  
  coverage.table <- data.frame(Test = Tests)
  n <- length(Tests)
  d <- length(Deltas)
  
  i <- 1
  for(i in 1:d) {
    coverage.table <- cbind(coverage.table, data[((i - 1)*n + 1):(i*n), 3])
  }
  
  names <- NULL
  for(i in 1:d) {
    names <- c(names, paste("Delta =", Deltas[i]))
  }
  
  colnames(coverage.table) <- c("Tests", names)
  return(coverage.table)
}



table.overall <- function(data.list, delta = 1) {
  #First item in list is description of distribution. Take this and remove it.
  description <- data.list[[1]]
  if(delta == 1) {
    description <- paste("Population: ", description, 
                         ": Type-I-Error rate.", sep = "")
  } else {
    description <- paste("Population: ", description, 
                         ": Power, delta = ", delta, ".", sep = "")
  }
  
  n.params <- length(data.list) - 1
  
  #Add columns for skewness and kurtosis.
  if(data.list[[1]] == "Normal (mu, sigma)") {
    params <- sapply(data.list[2:(n.params + 1)], function(x) as.numeric(x[1, 4:5]))
    sigma <- params[2, ]
    skew <- rep(0, length(params[2, ]))
    kurt <- rep(0, length(params[2, ]))
  } else if(data.list[[1]] == "t (df)") {
    params <- sapply(data.list[2:(n.params + 1)], function(x) as.numeric(x[1, 4]))
    sigma <- ifelse(params > 2, params / (params - 2), Inf)
    skew <- ifelse(params > 3, 0, NA)
    kurt <- ifelse(params > 4, 6/(params - 4), ifelse(params > 2, Inf, NA))
  } else if(data.list[[1]] == "Chi-sq (df)") {
    params <- sapply(data.list[2:(n.params + 1)], function(x) as.numeric(x[1, 4]))
    sigma <- 2 * params
    skew <- sqrt(8 / params)
    kurt <- 12 / params
  } else if(data.list[[1]] == "Gamma (shape, scale)") {
    params <- sapply(data.list[2:(n.params + 1)], function(x) as.numeric(x[1, 4:5]))
    sigma <- params[1, ] * params[2, ]^2
    skew <- 2 / sqrt(params[1, ])
    kurt <- 6 / params[1, ]
  } else if(data.list[[1]] == "Weibull (shape, scale)") {
    params <- sapply(data.list[2:(n.params + 1)], function(x) as.numeric(x[1, 4:5]))
    mu <- params[2, ] * gamma(1 + 1 / params[1, ])
    sigma <- sqrt(params[2, ]^2 * (gamma(1 + 2 / params[1, ]) - 
                                     (gamma(1 + 1 / params[1, ]))^2))
    skew <- (params[2, ]^3 * gamma(1 + 3 / params[1, ]) - 
               3 * mu * sigma^2 - 
               mu^3
            )/sigma^3
    kurt <- (params[2, ]^4 * gamma(1 + 4 / params[1, ]) -
               4 * skew * sigma^3 * mu -
               6 * mu^2 * sigma^2 -
               mu^4
            )/sigma^4 - 3
  } else if(data.list[[1]] == "Lognormal (lmu, lsigma)") {
    params <- sapply(data.list[2:(n.params + 1)], function(x) as.numeric(x[1, 4:5]))
    sigma <- sqrt((exp(params[2, ]^2) - 1) * exp(2 * params[1, ] + params[2, ]^2))
    skew <- (exp(params[2, ]^2) + 2) * sqrt(exp(params[2, ]^2) - 1)
    kurt <- exp(4 * params[2, ]^2) + 
      2 * exp(3 * params[2, ]^2) +
      3 * exp(2 * params[2, ]^2) -
      6
  } else if(data.list[[1]] == "Inverse-gamma (shape, scale)") {
    params <- sapply(data.list[2:(n.params + 1)], function(x) as.numeric(x[1, 4:5]))
    sigma <- ifelse(params[1, ] > 2, 
                    params[2, ]^(2) / 
                      ((params[1, ] - 1)^2 * (params[1, ] - 2)),
                    NA)
    skew <- ifelse(params[1, ] > 3, 
                   4 * sqrt(params[1, ] - 2) / (params[1, ] - 3),
                   NA)
    kurt <- ifelse(params[1, ] > 4, 
                   (30 * params[1, ] - 66) / 
                     ((params[1, ] - 3) * (params[1, ] - 4)),
                   NA)
  } else {
    stop("distribution name not found.")
  }
  
  data.list[[1]] <- NULL
  
  tests <- unique((data.list[[1]])[, 1])
  deltas <- unique((data.list[[1]])[, 2])
  
  n <- length(tests)      #Number of tests
  d <- length(deltas)     #Number of delta levels.
  num.dataframes <- length(data.list)

  coverage.table <- data.frame()
  i <- match(delta, deltas)
  for(data in data.list) {
    params <- NULL
    if(ncol(data) > 4) {
      params <- paste("(", data[1, 4], ",", data[1, 5], ")", sep = "")
    } else {
      params <- paste("(", data[1, 4], ")", sep = "")
    }
    coverage.table <- rbind(coverage.table, 
                            data.frame(params, 
                                       matrix(data[((i - 1)*n + 1):(i*n), 3], 
                                              nrow = 1)))
  }
  
  coverage.table <- cbind(coverage.table[, 1], sigma, skew, 
                          kurt, coverage.table[, -1])
  
  colnames(coverage.table) <- c("params", "sigma", "skew", 
                                "ex. kurt", levels(tests)[tests])

  return(list(desc = description, table = coverage.table))
}



get.table <- function(data, delta = 1) {
  if(is.data.frame(data)) {
    table.single.dist(data)
  } else {
    table.overall(data, delta)
  }
}

