#!/bin/bash

#SBATCH --job-name=running_rgi_1173_genomes_addedstrip_to_separatedfamilies
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=180g
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate rgi package
conda activate rgi

# make output directory for database and the actual species files
# creating a separate directory for the database ensures it wont affect subsequent commands
mkdir -p ~/all_gammaproteobacteria_data/rgi_output_1173_genomes/database_files

# move into the directory for the database files so they get stored there
cd ~/all_gammaproteobacteria_data/rgi_output_1173_genomes/database_files

# load card database
rgi clean --local

wget https://card.mcmaster.ca/latest/data

tar -xvf data ./card.json

rgi load --card_json ./card.json --local

# run rgi and create gene family files
python3 ~/scripts/rgi_scripts/run_RGI.py \
	-id ~/all_gammaproteobacteria_data/rgi_input_protein_fasta_1173_genomes \
        -od ~/all_gammaproteobacteria_data/rgi_output_1173_genomes/gene_family_files \
	-t 90 \
	-mj ~/scripts/rgi_scripts/all_1348_species_metrics.json \
        --remove_temp

# deactivate conda
conda deactivate

# get slurm job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
