#!/bin/bash
# This script will run busco on all of the different assemblies we have to get an idea of their quality scores and genome completeness. This will help decide if we use
# a lot of genomes that were not from the past 4 years. You will need to have downloaded the lineage dataset before running this script90190

#SBATCH --job-name=running_busco_correct_plotting2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=250g
#SBATCH --time=80:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate the conda env you are using with all the envs we need
conda activate pipeline_pckgs

# run the script that actuall runs busco
python3 ~/scripts/busco_quality_check.py \
	-nd ~/all_gammaproteobacteria_data/assemblies/nucleotide_fasta_files \
	-od ~/all_gammaproteobacteria_data/busco_output \
	-l ~/all_gammaproteobacteria_data/busco_lineage_dataset/lineages/gammaproteobacteria_odb10

# now run the busco plotting script
python3 ~/scripts/busco_plotting.py \
	-bd ~/all_gammaproteobacteria_data/busco_output \
	-od ~/all_gammaproteobacteria_data/busco_plotting_output_and_summary_files \
	-c 90 \
	-ani ~/all_gammaproteobacteria_data/output_data/final_deduplicated_gammaproteobacteria.txt

# run R script to create figures. first argument is the working directory, second is the path to the metrics file created by the script above
#Rscript figures_for_busco.R ~/all_gammaproteobacteria_data/busco_plotting_output_and_summary_files ~/all_gammaproteobacteria_data/busco_plotting_output_and_summary_files/sorted_by_complete_all_species_busco_metrics.tsv


# deactivate conda
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
