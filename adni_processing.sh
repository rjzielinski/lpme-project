#!/bin/bash
#SBATCH --job-name=adni_processing
#SBATCH --partition defq
#SBATCH -n 8
#SBATCH -o lpme_simulations_%j.out
#SBATCH -e lpme_simulations_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert_zielinski1@brown.edu
R CMD BATCH code/adni_mri_preprocess.R
