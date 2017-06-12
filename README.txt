To run the main simulation, see main.R. Simulation parameters (sample sizes, significance levels, delta levels, statistical tests to use, etc.) are set here and the script handles the rest.

The distribution parameters needed to obtain particular excess kurtosis values can be recalculated using the `script_find_params_for_kurtosis.R` script. These values are hand coded into the `simulation.R` file, which is used by `main.R`.

The code for each test statistic can be found in the `hypothesis_tests` folder.
Some results are saved as images:
	\figures
	 - Contains Power curves (.png files) from main simulation.
	 
	\tables
	 - Contains Type-I error rate tables (.png files) from main simulation.

Results are also saved as .csv files. These are located in the `results` folder:
	\results\verify_under_assumptions
	 - Contains the .csv files for the simulations of the six original test statistics with the test assumptions verified. See `script_verify_assumptions.R` to rerun this simulation.
 
	\results\simulation_type_I_errors
 	 - Contains .csv files of type-I errors from the main simulation.
 
	\results\simulation_power_delta_*
 	 - Contains .csv files for the power at *various delta levels from the main simulation.
