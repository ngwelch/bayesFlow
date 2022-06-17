#!/bin/bash

#SBATCH --job-name oos
#SBATCH -o log/oos.out
#SBATCH -e log/oos.out
#SBATCH --mem=32000
#SBATCH --time 0-10:00          # Wall time limit in D-HH:MM
#SBATCH --partition short #faloorum #SBATCH --reservation=faloorumgroup_fall
#SBATCH --ntasks 1 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nwelch@uw.edu

dir="/homes/nwelch/research/raftery/welchRaftery2022revision"
burninProportion="0"

module load R
Rscript --no-save $dir/src/heldout.R $dir $burninProportion

