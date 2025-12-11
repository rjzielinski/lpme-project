#!/bin/bash

#SBATCH --job-name=adni_lpme
#SBATCH --partition defq
#SBATCH -n 64
#SBATCH -o adni_lpme_%j.out
#SBATCH -e adni_lpme_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert_zielinski1@brown.edu

# Run the script
Rscript code/adni_application.R

# Bonus: write information about our job in the output file
scontrol show job $SLURM_JOB_ID
