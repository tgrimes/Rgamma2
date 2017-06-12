#######################################################################
#
# Main script that runs the simulation and exports the results.
#
#######################################################################

source("init_required_packages_and_files.R")
source("init_setup_directory.R")

# Simulation parameters.
n.set = c(10, 15, 20) # Sample sizes.
m = 10^5              # Number of simulations.
signf.set = c(0.05)   # Significance levels.
tail <- "left"        # Hypothesis test (only "left" has been fully implemented).

results <- NULL
timestart <- NULL
timeend <- NULL

#------------------------------------------------------------------------
# Run simulation with the five test statistics. 
#------------------------------------------------------------------------
ht.list <- c("ht.chisq", "ht.robust", "ht.z6", "ht.saddle", "ht.saddle.gamma2")
ht.names <- c("Chisq", "Robust", "R_lh", "R_gamma2", "Z6") # Used in tables/figures.
delta.set <- c(1, 2, 3, 4)

# NOTE: Set desired parameter values for population distributions 
#       in the run.simulation() function in the simulation.R script.
for(n in n.set) {
  for(significance.level in signf.set) {
    timestart <- format(Sys.time(), "%Y_%m_%d_%H.%M.%S")
    results <- run.simulation(n, m, significance.level, ht.list, delta.set, tail)
    timeend <- format(Sys.time(), "%Y_%m_%d_%H.%M.%S")
    
    #Edit test names. These are used in graphs and tables.
    for(i in 1:length(results)) {
      for(j in 2:length(results[[i]])) {
        levels(results[[i]][[j]][, 1]) <- ht.names
      }
    }
    
    #Save environment to disk. 
    filename <- paste(tail, "-tail_n", n, "_m", m, "_signf", significance.level, 
                      sep = "")
    dir.save <- paste(dir, "/results/", filename, sep = "")
    dir.create(dir.save)
    save.image(file = paste(dir.save, "/env", sep = ""))
    
    source("script_export_figures_tables.R")
    
    print(paste("n =", n, "- start:", timestart, "     end:", timeend))
  }
}
