#!/bin/bash

#SBATCH --job-name error 
#SBATCH -o log/error.out
#SBATCH -e log/error.out
#SBATCH --mem=4000
#SBATCH --time 0-00:15          # Wall time limit in D-HH:MM
#SBATCH --partition build #faloorum
#SBATCH --ntasks 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nwelch@uw.edu

dir="/homes/nwelch/research/raftery/welchRaftery2022revision"

module load R
Rscript --no-save $dir/src/inflow_fail.R "heldout"
Rscript --no-save $dir/src/inflow_fail.R "forecast"

