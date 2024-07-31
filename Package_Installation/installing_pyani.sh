#!/bin/bash

#SBATCH --job-name=installing_pyANI
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10g
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

#initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# create pyani env
conda create --name pyani_env python=3.8 -y

conda activate pyani_env

conda install pyani

conda install mummer blast legacy-blast fastani -y

pip install openpyxl

# deactivate conda
conda deactivate
