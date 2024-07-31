#!/bin/bash

#SBATCH --job-name=running_eggnog_mapper_1173genomes
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=95
#SBATCH --mem=300g
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

#initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate eggnog env
conda activate eggnog

# make directory for eggnog output
mkdir -p ~/all_gammaproteobacteria_data/eggnog_pangenome

cd ~/all_gammaproteobacteria_data/eggnog_pangenome

# download eggnog databases
download_eggnog_data.py -y

emapper.py -i /gpfs01/home/mbydr5/all_gammaproteobacteria_data/panta_output/panta_030pid_e7_LD07_split_1173_genomes_with_sflag_sotruenosplit/representative_clusters_prot.fasta \
	--dbmem \
	--output_dir ~/all_gammaproteobacteria_data/eggnog_pangenome \
	-o whole_dataset \
	-m diamond \
	--cpu 95 \
	--tax_scope Gammaproteobacteria \
	--tax_scope_mode Bacteria

# deactivate conda
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
