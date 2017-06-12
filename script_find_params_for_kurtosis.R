#######################################################################
#
# This script finds the parameter values required for a desired kurtosis via 
# the `uniroot` root finding function. For distributions with >1 parameters, 
# any location and scale parameters are fixed to 0 and 1, respectively. 
# Results are printed to the screen.
#
#######################################################################

kurt_set <- c(15, 10, 9, 8, 7, 6, 5, 4, 1)
vals <- matrix(0, 6, length(kurt_set))
rownames(vals) <- c("t(df)", "chisq(df)", "gamma(shape, scale = 1)",
                    "lnorm(lmu = 0, tau)", "weib(shape, scale = 1)", 
                    "invgamma(shape, scale = 1)")
colnames(vals) <- kurt_set

# Kurtosis function for gamma, weibull, lognormal, and inv.gamma distributions.
kurt_t <- function(df) {
  params <- df
  kurt <- ifelse(params > 4, 6/(params - 4), ifelse(params > 2, Inf, NA))
  as.numeric(kurt)
}
kurt_chisq <- function(df) {
  params <- df
  kurt <- 12 / params
  kurt
}
kurt_gamma <- function(shape, scale = rep(1, length(shape))) {
  params <- rbind(shape, scale)
  kurt <- 6/params[1, ]
  kurt
}
kurt_weib <- function(shape, scale = rep(1, length(shape))) {
  params <- rbind(shape, scale)
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
  as.numeric(kurt)
}
kurt_ln <- function(lsigma, lmu = rep(0, length(lsigma))) {
  params <- rbind(lmu, lsigma)
  kurt <- exp(4 * params[2, ]^2) + 
    2 * exp(3 * params[2, ]^2) +
    3 * exp(2 * params[2, ]^2) -
    6
  kurt
}
kurt_ig <- function(shape, scale = rep(1, length(shape))) {
  params <- rbind(shape, scale)
  kurt <- ifelse(params[1, ] > 4, 
         (30 * params[1, ] - 66) / 
           ((params[1, ] - 3) * (params[1, ] - 4)),
         NA)
  kurt
}

for(i in 1:length(kurt_set)) {
  vals[1, i] <- uniroot(function(x) kurt_t(x) - kurt_set[i], c(2.01, 100))$root
}
kurt_t(vals[1, ])

for(i in 1:length(kurt_set)) {
  vals[2, i] <- uniroot(function(x) kurt_chisq(x) - kurt_set[i], c(0.1, 100))$root
}
kurt_chisq(vals[2, ])

for(i in 1:length(kurt_set)) {
  vals[3, i] <- uniroot(function(x) kurt_gamma(x) - kurt_set[i], c(0.5, 5))$root
}
kurt_gamma(vals[3, ])

for(i in 1:length(kurt_set)) {
  vals[4, i] <- uniroot(function(x) kurt_weib(x) - kurt_set[i], c(0.01, 5))$root

}
kurt_weib(vals[4, ])

for(i in 1:length(kurt_set)) {
  vals[5, i] <- uniroot(function(x) kurt_ln(x) - kurt_set[i], c(4.01, 100))$root
}
kurt_ln(vals[5, ])

for(i in 1:length(kurt_set)) {
  vals[6, i] <- uniroot(function(x) kurt_ig(x) - kurt_set[i], c(0.01, 100))$root
}
kurt_ig(vals[6, ])

print("Results:")
round(vals, 3)
