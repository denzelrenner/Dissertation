#!/bin/bash

#SBATCH --job-name=installing_R
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=100g
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# create a new conda environment for all our pipeline packages
conda create -y -n R_env r-essentials r-base

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate environment with the conda packages
conda activate R_env

# add channels
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

# install R packages
conda install -c conda-forge r-tidyverse


# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
