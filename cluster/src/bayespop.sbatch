#!/bin/bash

#SBATCH --job-name bayespop 
#SBATCH -o log/bayespop_%a.out
#SBATCH -e log/bayespop_%a.out
#SBATCH --array=1-10
#SBATCH --mem=8000
#SBATCH --partition faloorum
#SBATCH --ntasks 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nwelch@uw.edu

dir="/homes/nwelch/research/raftery/welchRaftery2022revision"
batchSize="100" #100

module load R
Rscript --no-save $dir/src/bayespop.R $batchSize $SLURM_ARRAY_TASK_ID 

