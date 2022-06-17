# Run Notes

There are several steps to replicating the flow forecasting results. 
This directory is organized as follows:

	- bayesPop: source codes and data for the bayesPop project
	- bayesPop.output: population forecasts based on migration flow model
	- data: data needed to generate migration adjusted population forecasts
	- forecast: posterior parameter samples using the most updated WPP and estimated flows
	- heldout: posterior parameter samples using data as of WPP 2015 release
	- lib: local R libraries
	- log: sampler log files
	- README.md: this file
	- src: all flow forecasting codes organized by SLURM job and function, e.g. inflow_heldout.sbatch runs inflow_heldout.R.

Codes to fit the flow model reside in src.
The heldout analysis is executed in a series of steps:

	1) Fit outflow model component. This script determines the MCMC sample size and thinning interval.
		sbatch src/outflow.sbatch
	2) Fit the inflow model component. This script executes 200 parallel jobs. 
		sbatch src/inflow_heldout.sbatch
	3) Combine the results from 1-2 and summarize the fit quality. Note that SDN is used as a proxy for SDN preferences.
		sbatch src/heldout.sbatch

Codes to generate forecasts also reside in src.
Forecasts are generated according to the following steps:
	
	1) Fit the outflow model component. This step is actually included as part of the heldout analysis.
		[executed in step 1 above]
	2) Fit the inflow model component. This script executes 200 parallel jobs.
		sbatch src/inflow_forecast.sbatch
	3) Generate probabilistic population forecasts. This step generates 1000 trajectories in 10 batches in parallel.
		sbatch src/bayespop.sbatch
	4) Combine the forecast results from the 10 batches and summarize the results. 
		sbatch src/bayespopSummary.sbatch


