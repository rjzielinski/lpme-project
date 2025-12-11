#!/bin/bash
#SBATCH --job-name=adni_volumes
#SBATCH --partition defq
#SBATCH -n 64
#SBATCH -o adni_volumes_%j.out
#SBATCH -e adni_volumes_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert_zielinski1@brown.edu
R CMD BATCH code/adni_volumes.R
