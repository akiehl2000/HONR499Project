#! /bin/bash

#SBATCH --account=csu-general
#SBATCH --qos normal
#SBATCH --job-name=fitINGARCH
#SBATCH --output=../Output/INGARCH/%x_%a_%j.out
#SBATCH --error=../Output/INGARCH/%x_%a_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=adam.kiehl@colostate.edu
#SBATCH --mail-type=ALL

source ~/../../curc/sw/anaconda/default
conda activate fitEnv2

Rscript ../Scripts/INGARCH_models.R