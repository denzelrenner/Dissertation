#!/bin/bash

#SBATCH --job-name=running_goatools_creating_input
#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=95
#SBATCH --mem=200g
#SBATCH --time=05:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

mkdir -p ~/all_gammaproteobacteria_data/goatools_output

conda deactivate

conda activate pipeline_pckgs

python3 ~/scripts/goatools_scripts/create_goatools_input_files.py \
	--output_directory ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes \
	--positively_correlated_genes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE_WITH_GROUPS/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families/master_positive_correlated_genes.tsv \
	--negatively_correlated_genes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE_WITH_GROUPS/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families/master_negative_correlated_genes.tsv \
	--Rtab_file ~/all_gammaproteobacteria_data/panta_output/panta_030pid_e7_LD07_split_1173_genomes_with_sflag_sotruenosplit/gene_presence_absence.Rtab \
	--annotation_file ~/all_gammaproteobacteria_data/eggnog_pangenome/whole_dataset.emapper.annotations

conda deactivate

# activate goatools env
conda activate goatools_env

# move into output dir
cd ~/all_gammaproteobacteria_data/goatools_output

# get obo file
wget http://current.geneontology.org/ontology/go-basic.obo
wget https://current.geneontology.org/ontology/subsets/goslim_prokaryote.obo

# --no_propagate_counts
# run goatools
find_enrichment.py --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --obsolete replace --goslim goslim_prokaryote.obo --obo go-basic.obo --alpha 0.05 \
	--outfile positive_results.tsv --indent ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/positively_correlated_genes_with_GOterm.txt \
	~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/pangenome_population_only.txt ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/pangenome_association_file.tsv

find_enrichment.py --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --obsolete replace --goslim goslim_prokaryote.obo --obo go-basic.obo --alpha 0.05 \
        --outfile negative_results.tsv --indent ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/negatively_correlated_genes_with_GOterm.txt \
	~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/pangenome_population_only.txt ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/pangenome_association_file.tsv

# deactivate env
conda deactivate

# activate env
conda activate pipeline_pckgs

# This script produces the files needed to then plot the enriched terms in R
python3 ~/scripts/goatools_scripts/plot_enriched_terms.py \
	--output_directory ~/all_gammaproteobacteria_data/goatools_output/stats_and_plots \
	--negative_goatools ~/all_gammaproteobacteria_data/goatools_output/negative_results.tsv \
	--positive_goatools ~/all_gammaproteobacteria_data/goatools_output/positive_results.tsv \
	--number_of_terms 15

# deactivate conda env
conda deactivate

# activate R env
conda activate R_env

Rscript ~/all_gammaproteobacteria_data/goatools_output/stats_and_plots/top_GOs_plotting.R

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"



# old run command
#find_enrichment.py --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --obsolete replace --goslim goslim_prokaryote.obo --obo go-basic.obo --alpha 0.05 \
#	--outfile positive_results.tsv --indent ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/positively_correlated_genes_with_GOterm.txt \
#	~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/pangenome_population_only.txt ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/pangenome_association_file.tsv

#find_enrichment.py --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --obsolete replace --goslim goslim_prokaryote.obo --obo go-basic.obo --alpha 0.05 \
#        --outfile negative_results.tsv --indent ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/negatively_correlated_genes_with_GOterm.txt \
#	~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/pangenome_population_only.txt ~/all_gammaproteobacteria_data/goatools_input_BH_1173genomes/pangenome_association_file.tsv
#	~/all_gammaproteobacteria_data/goatools_in_pangenome_annotation/pangenomes_pop_genes_only.txt ~/all_gammaproteobacteria_data/goatools_in_pangenome_annotation/pangenome_association_file.tsv
