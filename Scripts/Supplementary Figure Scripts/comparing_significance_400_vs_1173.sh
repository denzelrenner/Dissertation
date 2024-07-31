#!/bin/bash

#SBATCH --job-name=comparing_significance
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=95
#SBATCH --mem=50g
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate the conda env you are using with all the envs we need
conda activate pipeline_pckgs

python3 ~/scripts/scoary_scripts/compare_400_vs_1173.py \
	--output_directory ~/all_gammaproteobacteria_data/comparing_400_1173_genomes \
	--benjaminiH_1173genomes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families \
	--benjaminiH_400genomes ~/all_gammaproteobacteria_data/scoary_output_400genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families

# deactivate conda env
conda deactivate

# activate R plotting env
conda activate R_env

Rscript ~/all_gammaproteobacteria_data/comparing_400_1173_genomes/combined_pvalue_frequency_plot.R

conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
