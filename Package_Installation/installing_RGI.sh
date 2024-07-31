#!/bin/bash

#SBATCH --job-name=installing_rgi
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=10g
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# follow guidelines on RGI github page
conda create --name rgi -y

conda activate rgi

conda install --channel bioconda rgi -y

# deactivate conda enc
conda deactivate


