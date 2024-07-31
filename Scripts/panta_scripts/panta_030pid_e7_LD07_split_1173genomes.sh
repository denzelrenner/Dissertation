#!/bin/bash

#SBATCH --job-name=panta_030pid_e7_LD07_split_1173_genomes_withsflagsotruenosplit
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=95
#SBATCH --mem=200g
#SBATCH --time=168:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate panta

# run panta
panta main -g ~/all_gammaproteobacteria_data/prokka_output_files_busco_filtered_input/all_species_gff/*.gff \
	-e 1E-7 \
	-i 0.3 \
	-b blast \
	-t 95 \
	--LD 0.7 \
	-s \
	-o ~/all_gammaproteobacteria_data/panta_output/panta_030pid_e7_LD07_split_1173_genomes_with_sflag_sotruenosplit

# deactivate conda env
conda deactivate

# get slurm job id
echo "The Job ID for this job is: $SLURM_JOB_ID"

