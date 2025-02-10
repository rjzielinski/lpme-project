#!/bin/bash
#SBATCH --job-name=lpme_simulations
#SBATCH --partition defq
#SBATCH -n 64
#SBATCH -o lpme_simulations_%j.out
#SBATCH -e lpme_simulations_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert_zielinski1@brown.edu
R CMD BATCH code/run_simulations.R
