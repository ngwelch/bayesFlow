#!/bin/bash

#SBATCH --job-name bpsum
#SBATCH -o log/bpsum.out
#SBATCH -e log/bpsum.out
#SBATCH --mem=32000
#SBATCH --time 0-01:00          # Wall time limit in D-HH:MM
#SBATCH --partition short
#SBATCH --ntasks 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nwelch@uw.edu

dir="/homes/nwelch/research/raftery/welchRaftery2022revision"

module load R
Rscript --no-save $dir/src/bayespopSummary.R 

