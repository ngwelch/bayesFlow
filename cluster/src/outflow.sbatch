#!/bin/bash

#SBATCH --job-name outflow
#SBATCH -o log/outflow.out
#SBATCH -e log/outflow.out
#SBATCH --mem=8000
#SBATCH --time 0-01:00          # Wall time limit in D-HH:MM
#SBATCH --partition short #faloorum #SBATCH --reservation=faloorumgroup_fall
#SBATCH --ntasks 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nwelch@uw.edu

dir="/homes/nwelch/research/raftery/welchRaftery2022revision"
mcmcSS="10000"
mcmcThin="1"


module load R
Rscript --no-save $dir/src/outflow_heldout.R $dir $mcmcSS $mcmcThin 
Rscript --no-save $dir/src/outflow_forecast.R $dir $mcmcSS $mcmcThin 

