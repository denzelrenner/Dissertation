#!/bin/bash

#SBATCH --job-name=COG_category_analysis
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=95
#SBATCH --mem=20g
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate pipeline_pckgs

python3 ~/scripts/cog_scripts/cog_calculations_and_plots.py \
	-od ~/all_gammaproteobacteria_data/COG_calculation_1173_genomes \
	--annotation_file ~/all_gammaproteobacteria_data/eggnog_pangenome/whole_dataset.emapper.annotations \
	--positively_correlated_genes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families/master_positive_correlated_genes.tsv \
	--negatively_correlated_genes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families/master_negative_correlated_genes.tsv \
	--Rtab_file ~/all_gammaproteobacteria_data/panta_output/panta_030pid_e7_LD07_split_1173_genomes_with_sflag_sotruenosplit/gene_presence_absence.Rtab \
	--pvalue 0.05 \
	--number_of_tests 22

# deativate conda env
conda deactivate

# activate conda env
conda activate R_env

Rscript ~/all_gammaproteobacteria_data/COG_calculation_1173_genomes/COG_plotting.R

# deactivate env
conda deactivate

