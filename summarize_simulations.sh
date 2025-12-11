#!/bin/bash
#SBATCH --job-name=simulation_summary
#SBATCH --partition defq
#SBATCH -n 1
#SBATCH -o simulation_summary_%j.out
#SBATCH -e simulation_summary_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert_zielinski1@brown.edu
R CMD BATCH code/simulation_summary.R
