#!/bin/bash

#SBATCH --job-name heldout 
#SBATCH -o log/heldout_%a.out
#SBATCH -e log/heldout_%a.out
#SBATCH --array=1-200
#SBATCH --mem=16000
#SBATCH --time 0-08:00          # Wall time limit in D-HH:MM
#SBATCH --partition short #faloorum #SBATCH --reservation=faloorumgroup_fall
#SBATCH --ntasks 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nwelch@uw.edu

dir="/homes/nwelch/research/raftery/welchRaftery2022revision"

module load R
Rscript --no-save $dir/src/inflow_heldout.R $dir $SLURM_ARRAY_TASK_ID
#Rscript --no-save $dir/src/inflow_heldout.R $dir $mcmcSS $mcmcThin $traceSize $SLURM_ARRAY_TASK_ID
