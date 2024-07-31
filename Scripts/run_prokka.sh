#!/bin/bash

# This script is going to run roary on our final list of organisms after we had deduplicated the dataset based on the pairwise ANI values 

#SBATCH --job-name=running_prokka
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=150g
#SBATCH --time=60:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate the conda env you are using with all the envs we need
conda activate pipeline_pckgs

# run the python script
python3 ~/scripts/running_prokka_for_roary.py \
	-nd ~/all_gammaproteobacteria_data/assemblies/nucleotide_fasta_files \
	-if ~/all_gammaproteobacteria_data/busco_plotting_output_and_summary_files/busco_filtered_gammaproteobacteria.txt \
	-pd ~/all_gammaproteobacteria_data/assemblies/protein_fasta_files \
	-o ~/all_gammaproteobacteria_data/prokka_output_files_busco_filtered_input

# run roary
#roary -p 90 -f ~/all_gammaproteobacteria_data/roary_output ~/all_gammaproteobacteria_data/prokka_output_files_busco_filtered_input/*.gff

# deactivate conda
conda deactivate

# get slurm job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
