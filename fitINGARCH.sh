#! /bin/bash

#SBATCH --account=csu-general
#SBATCH --qos normal
#SBATCH --job-name=fitINGARCH
#SBATCH --output=%x_%a_%j.out
#SBATCH --error=%x_%a_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --mail-user=adam.kiehl@colostate.edu
#SBATCH --mail-type=ALL

source ~/../../curc/sw/anaconda/default
conda activate fitEnv2

Rscript summitTest.R