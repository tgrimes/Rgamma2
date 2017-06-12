#Get working directory. Required source files should also be located here.
dir <- getwd()

#Required files:
packages <- c("distr", "ggplot2", "gridExtra", "MASS", 
              "moments", "polynom", "pracma", "actuar")
scripts <- c("graph.coverage.R", "table.coverage.R", "simulation.R", 
             "sourceDir.R")
directories <- c("hypothesis_tests")

#Install any missing packages.
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

#Load any packages that are not already loaded.
if(!all(packages %in% loadedNamespaces())) {
  missing_packages <- setdiff(packages, loadedNamespaces())
  sapply(missing_packages, require, character.only = TRUE)
}

#Load scripts from source files.
for(script in scripts) {
  source(paste(dir, "/", script, sep = ""))
}

#Load scripts from other directories.
for(directory in directories) {
  sourceDir(paste(dir, "/", directory, sep = ""))
}
