#!/bin/bash

# This script is going to produce a traits file for scoary by taking in the gene absence file from roary. The trait we are interested in based on our research question has to do with the presence
# or absence of a given antibiotic gene.

#SBATCH --job-name=running_scoary_groups_more_time_1173_genomes
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=95
#SBATCH --mem=350g
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate the conda env you are using with all the envs we need
conda activate pipeline_pckgs

# make directory for the csvs scoary produces scoary output
mkdir -p ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/csv_output_files_bonferroni
mkdir -p ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/csv_output_files_benjamini_hochberg
mkdir -p ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/trait_files
mkdir -p ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/input_files

# create the traits file to use in scoary
python3 ~/scripts/scoary_scripts/create_scoarytraits_file_genefamily_batched_files.py \
	-r ~/all_gammaproteobacteria_data/panta_output/panta_030pid_e7_LD07_split_1173_genomes_with_sflag_sotruenosplit/gene_presence_absence.Rtab \
	-j ~/all_gammaproteobacteria_data/rgi_output_1173genomes/gene_family_files/all_gene_family_data.json \
	-gf ~/all_gammaproteobacteria_data/rgi_output_1173genomes/gene_family_files/unique_gene_families.txt \
	-od ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/trait_files

# convert panta file to roary csv. Takes the rtab file and an output directory as arguments
python3 ~/scripts/scoary_scripts/panta_roary_conversion.py \
	-r ~/all_gammaproteobacteria_data/panta_output/panta_030pid_e7_LD07_split_1173_genomes_with_sflag_sotruenosplit/gene_presence_absence.Rtab \
        -od ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/input_files

#echo 'conversion succesful'

# run scoary, -t is traits file, -g is presence absence csv -p is p value cut off,-s is where data for species starts
# run scoary with bonferroni
for i in ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/trait_files/*.trait; do
	echo "Running trait file $i"
	scoary \
		-t "$i" \
		-g ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/input_files/final_gene_presence_absence.Rtab.csv \
		--threads 95 \
		-o ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/csv_output_files_bonferroni \
		-s 4 \
		--no_pairwise \
		-c B -p 0.0003
done

# run scoary with benjamini correction
for i in ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/trait_files/*.trait; do
        echo "Running trait file $i"
        scoary \
                -t "$i" \
                -g ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/input_files/final_gene_presence_absence.Rtab.csv \
                --threads 95 \
                -o ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/csv_output_files_benjamini_hochberg \
                -s 4 \
                --no_pairwise \
                -c BH -p 0.0003
done

# parse the scoary output and produce files that allow us to plot in R amongst other things
python3 ~/scripts/scoary_scripts/parsing_scoary_output_NOPAIRWISE_corrected.py \
	--csv_dir ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/csv_output_files_bonferroni \
	--output_directory ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_bonferroni \
	-j ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/trait_files/trait_to_regular_name.json \
	--pvalue_correction Bonferroni_p


# parse the scoary output and produce files that allow us to plot in R amongst other things
python3 ~/scripts/scoary_scripts/parsing_scoary_output_NOPAIRWISE_corrected.py \
        --csv_dir ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/csv_output_files_benjamini_hochberg \
        --output_directory ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg \
        -j ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/trait_files/trait_to_regular_name.json \
        --pvalue_correction Benjamini_H_p

#python3 ~/scripts/scoary_scripts/cytoscape_input.py \
 #       --input_file_positive ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families/positive_correlated_genes_and_genefamily_only.json \
 #       --input_file_negative ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families/negative_correlated_genes_and_genefamily_only.json \
  #      --output_directory ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/cytoscape_input_files_benjamini_hochberg \


# creates data for plotting
python3 ~/scripts/scoary_scripts/create_significance_comparison_script.py \
	--output_directory ~/all_gammaproteobacteria_data/comparison_of_significance \
	--benjaminiH_1173genomes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families \
	--bonferroni_1173genomes ~/all_gammaproteobacteria_data/scoary_output_1173genomes_NOPAIRWISE/parsed_output_and_plotting_files_bonferroni/across_all_gene_families \
	--benjaminiH_400genomes ~/all_gammaproteobacteria_data/scoary_output_400genomes_NOPAIRWISE/parsed_output_and_plotting_files_benjamini_hochberg/across_all_gene_families \
	--bonferroni_400genomes ~/all_gammaproteobacteria_data/scoary_output_400genomes_NOPAIRWISE/parsed_output_and_plotting_files_bonferroni/across_all_gene_families

conda deactivate

# activate R plotting env
conda activate R_env

Rscript ~/all_gammaproteobacteria_data/comparison_of_significance/combined_pvalue_frequency_plot.R

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
