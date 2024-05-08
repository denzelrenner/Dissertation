#!/bin/bash

#SBATCH --job-name=installing_packages
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20g
#SBATCH --time=3:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# create a new conda environment for all our pipeline packages
conda create -y -n pipeline_pckgs python=3.7

#initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate environment with the conda packages
conda activate pipeline_pckgs

# install prokka
conda install -y -c bioconda prokka

# install busco
conda install -y -c bioconda busco

# install roary
conda install -y -c bioconda roary=3.12.0

# install macsy finder
conda install -y -c bioconda macsyfinder

# install the ncbi datasets package
conda install -y -c conda-forge ncbi-datasets-cli

# install scoary and its dependencies
pip install scipy scoary ete3 six

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
