#!/bin/bash
#SBATCH --job-name=sim_volumes
#SBATCH --partition defq
#SBATCH -n 64
#SBATCH -o simulation_volumes_%j.out
#SBATCH -e simulation_volumes_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert_zielinski1@brown.edu
R CMD BATCH code/simulation_atrophy_analysis.R
