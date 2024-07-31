#!/bin/bash

#SBATCH --job-name=deduplicate_dataset
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=60g
#SBATCH --time=3:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate conda env with all the packages
conda activate pyani_env

# run script
python3 ~/scripts/deduplicate_dataset.py \
	-al ~/all_gammaproteobacteria_data/ANIm_output/ANIm_alignment_lengths.tab \
	-cov ~/all_gammaproteobacteria_data/ANIm_output/ANIm_alignment_coverage.tab \
	--ani ~/all_gammaproteobacteria_data/ANIm_output/ANIm_percentage_identity.tab \
	-o ~/all_gammaproteobacteria_data/output_data/final_deduplicated_gammaproteobacteria.txt

# deactivate conda
conda deactivate
