#! /bin/bash

#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --job-name=fitINGARCH

module purge

module load anaconda

Rscript ./Scripts/INGARCH_models.Rmd