#######################################################################
#
# This script is used by main.R to export the simulation tables/figures. It 
# can also be used independently by first loading a saved environment (see
# Option below).
#
#######################################################################

# Optional: load a saved environment. This is used if we wish to recreate
# tables or graphs, or to check the results from previous simulations. Note,
# everything from the simulation is saved, including the functions used. 
# For example, after loading a previous environment, type "ht.saddle.gamma2"
# in the console to see the exact code used to compute that hypothesis test.
#
#load(paste(dir, "results/left-tail_n10_m1000_signf0.05/env", sep = "/"))

save_table_type_I_error <- TRUE
save_csv_result_tables <- TRUE
save_figure_power_curve <- TRUE

#Print out tables of results and export as .png image and/or .csv file.
for(delta in delta.set) {
  for(coverage.list in results) {
    tab <- get.table(coverage.list, delta = delta)
    tab$table[, -1] <- round(tab$table[, -1], 4)
    tab$table[, 2:4] <- format(tab$table[, 2:4], digits = 2)
    print(paste(tab[[1]], "significance =", significance.level))
    print((tab[[2]])[,])
    
    # Optional: save table as .png image. This is meant for type-I errors.
    # The file directory and/or file name should be adjusted if used for
    # power tables (delta > 1).
    if(save_table_type_I_error && delta == 1) {
      file_dir <- paste("tables/", coverage.list[[1]],
                        " n=", n, " m=", m, ".png", sep = "")
      if(length(coverage.list) > 7) {
        png(file_dir, 2166, 900, res = 300)
      }
      else {
        png(file_dir, 2166, 646, res = 300)
      }
      
      p <- tableGrob(tab[[2]])
      grid.arrange(p)
      dev.off()
    }
    
    #Save partial table as .csv file.
    if(save_csv_result_tables) {
      if(delta == 1) {
        file_dir <- paste("results/simulation_type_I_errors/", coverage.list[[1]], 
                          " n=", n, " m=", m, ".csv", sep = "")
      } else {
        file_dir <- paste("results/simulation_power_delta_", delta, "/", 
                          coverage.list[[1]], 
                          " n=", n, " m=", m, ".csv", sep = "")
      }
      write.csv((tab[[2]])[, -c(2, 5)], file = file_dir)
    }
  }
}

if(save_figure_power_curve) {
  #Print out graphs of results:
  for(coverage.list in results) {
    tests <- c("Robust", "R_lh", "R_gamma2", "Z6")
    file_dir <- paste("figures/", coverage.list[[1]], 
                      " n=", n, " m=", m, ".png", sep = "")
    if(length(coverage.list) > 7) {
      png(file_dir, 3000, 3000*(9/10), res = 395)
    }
    else {
      png(file_dir, 3000, 3000*(2/3), res = 375)
    }
    graph.overall(coverage.list, tests, significance.level)
    dev.off()
  }
}