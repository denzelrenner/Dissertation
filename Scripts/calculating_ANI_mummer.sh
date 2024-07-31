#!/bin/bash

# This script first gets all the nucleotide and protein fasta files from the assemblies we downloaded and creates a unique directory for each of them. It then uses the assembly fasta files to calculate ANI in pairwise contrasts
#SBATCH --job-name=calculating_ANI_mummer
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=94
#SBATCH --mem=60g
#SBATCH --time=60:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate pyANI conda environment
conda activate pyani_env

# run pyANI on the folder with all of our fasta files
# using ANIm (mummer)
average_nucleotide_identity.py -i ~/all_gammaproteobacteria_data/assemblies/nucleotide_fasta_files/ \
	-o ~/all_gammaproteobacteria_data/ANIm_output -m ANIm -g -v --write_excel --gformat pdf,png --workers 94

# deactivate conda 
conda deactivate

# echo JOBID
echo "The Job ID for this job is: $SLURM_JOB_ID"
