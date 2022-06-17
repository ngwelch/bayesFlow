#!/bin/bash

#SBATCH --job-name forecast 
#SBATCH -o log/forecast_%a.out
#SBATCH -e log/forecast_%a.out
#SBATCH --array=1-200 #24,47,133,157
#SBATCH --mem=8000
#SBATCH --time 0-10:00          # Wall time limit in D-HH:MM
#SBATCH --partition short #faloorum #SBATCH --reservation=faloorumgroup_fall
#SBATCH --ntasks 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nwelch@uw.edu

dir="/homes/nwelch/research/raftery/welchRaftery2022revision"
#mcmcSS="20000"
#mcmcThin="10"


module load R
Rscript --no-save $dir/src/inflow_forecast.R $dir $SLURM_ARRAY_TASK_ID
#Rscript --no-save $dir/src/inflow_forecast.R $dir $mcmcSS $mcmcThin $SLURM_ARRAY_TASK_ID
