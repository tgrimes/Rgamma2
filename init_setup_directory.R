

create_dir <- function(dir_name) {
  if(dir.exists(dir_name)) {
    cat("\t", dir_name, "EXISTS.\n")
  } else {
    dir.create(dir_name)
    cat("\t", dir_name, "CREATED.\n")
  }
}

cat("Setting up directories:\n")
create_dir("results")
create_dir("results/simulation_power_delta_2")
create_dir("results/simulation_power_delta_2")
create_dir("results/simulation_power_delta_2")
create_dir("results/simulation_type_I_errors")
create_dir("results/verify_under_assumptions")
create_dir("tables")
create_dir("figures")