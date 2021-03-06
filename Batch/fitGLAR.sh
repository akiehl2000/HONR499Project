#! /bin/bash

#SBATCH --account=csu-general
#SBATCH --qos normal
#SBATCH --job-name=fitGLAR
#SBATCH --output=../Output/GLAR/%x_%a_%j.out
#SBATCH --error=../Output/GLAR/%x_%a_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=adam.kiehl@colostate.edu
#SBATCH --mail-type=ALL

source ~/../../curc/sw/anaconda/default
conda activate fitEnv2

Rscript ../Scripts/GLAR_models.R