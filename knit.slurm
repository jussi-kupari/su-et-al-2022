#!/bin/bash -l
#SBATCH -A snic2021-22-646
#SBATCH -p node  
#SBATCH -n 1
#SBATCH -t 12:00:00

# Load modules
module load bioinfo-tools
module load R/4.1.1

# Make sure we are in the right directory
cd /crex/proj/uppstore2017283/jussi/Jie2022  || exit

# Run R script
Rscript --vanilla ./R/knit_notebooks.R

  
