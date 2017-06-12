#######################################################################
#
# This script is not used by any other file and is meant to be run
# separately. Requires that load_required_packages_and_files.R has been
# sourced.
#
# The asymptotic properties for six test statistics are verified. The tests are
# performed with their distribution assumptions and known nuisance parameter 
# assumptions satisfied. Left-tailed tests of variance are conducted to assess 
# the type-I error rate (delta <= 1) and power (delta > 1). Resulting tables are
# created and stored as .csv files in results/verify_under_assumptions.
#
# The code can be rerun with different values of n, m, or significance.level. 
#
#######################################################################

# Run this line if packages/files have not already been loaded:
#source(load_required_packages_and_files.R)

n <- 100
m <- 10^4
significance.level <- 0.05
tail <- "left"
save_directory <- "results/verify_under_assumptions"

results <- vector("list", 6)
results[[1]] <- verify.normal(n, m, significance.level, tail)
results[[2]] <- verify.chisq(n, m, significance.level, tail)
results[[3]] <- verify.exp(n, m, significance.level, tail)
results[[4]] <- verify.gamma(n, m, significance.level, tail)
results[[5]] <- verify.weibull(n, m, significance.level, tail)
results[[6]] <- verify.lognormal(n, m, significance.level, tail)

export <- function(test_name, results_table) {
  file_name <- paste(test_name, " n=", n, " m=", m, " signif=", 
                     significance.level, ".csv", sep = "")
  write.csv(results_table, paste(save_directory, file_name, sep = "/"))
}
library(reshape2)
export("normal", dcast(results[[1]], mu + sigma ~ delta, value.var = "coverage"))
export("chisq", dcast(results[[2]], df ~ delta, value.var = "coverage"))
export("exp", dcast(results[[3]], lambda ~ delta, value.var = "coverage"))
export("gamma", dcast(results[[4]], alpha + beta ~ delta, value.var = "coverage"))
export("weibull", dcast(results[[5]], shape + scale ~ delta, value.var = "coverage"))
export("lognormal", dcast(results[[6]], lmu + lsigma ~ delta, value.var = "coverage"))
