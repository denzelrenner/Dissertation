#!/bin/bash

#SBATCH --job-name=installing_goatools
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=70
#SBATCH --mem=40g
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# create a new conda environment for all our pipeline packages
conda create -y -n goatools_env

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate environment with the conda packages
conda activate goatools_env

# add channels
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

# install roary
conda install -c bioconda goatools

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
