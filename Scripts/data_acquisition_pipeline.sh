#!/bin/bash

# This script runs the whole data acquisition pipeline for the data acquisition stage of the analysis. It first gets all the different complete and refseq (GCF) assemblies for different gammaproteobcteria from ncbi
# The accession code , organism name, assembly release data, completeness, and a whole battery of other metrics are for each species is written to an output file
# Then the resulting file is filtered to remove duplicate assemblies, and between different strains of the same species, take the most recent assembly.
# All the different assemblies for our final list of gammaproteobacteria are then downloaded as zip files, and we put metrics like assembly length into a txt file.
# An R script is also ran to produce a plot showing the distribution of the GC content and assembly lengths which would become useful when doing phylogenetic trees because sometimes considering GC content is important


#SBATCH --job-name=acquiring_data
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=95g
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs01/home/mbydr5/OnE/%x.out
#SBATCH --error=/gpfs01/home/mbydr5/OnE/%x.err

# initialise conda
source ~/miniconda3/etc/profile.d/conda.sh

# activate the conda env you are using with all the envs we need
conda activate pipeline_pckgs

# make a directory to put all the output into
mkdir -p ~/all_gammaproteobacteria_data/assemblies
mkdir -p ~/all_gammaproteobacteria_data/metric_files
mkdir -p ~/all_gammaproteobacteria_data/output_data

# grab all gammaproteobacteria hits from the ncbi website. 1236 is the taxon id for gammaproteobacteria
datasets summary genome taxon 1236 --assembly-level complete --assembly-source 'RefSeq' --as-json-lines | dataformat tsv genome \
	--fields accession,organism-name,assminfo-release-date,checkm-completeness,assmstats-number-of-contigs,assmstats-number-of-scaffolds,organism-tax-id,checkm-species-tax-id,assminfo-biosample-description-organism-tax-id > \
	~/all_gammaproteobacteria_data/output_data/all_gammaprotobacteria.txt

# filter the data to remove duplicate strains by taking the most recent assembly or if that fails then taking the best checkM completeness
python3 ~/scripts/filtering_gammaprotobacteria.py \
	-i ~/all_gammaproteobacteria_data/output_data/all_gammaprotobacteria.txt \
	-d ~/all_gammaproteobacteria_data/output_data

# move into the assemblies directory so all the assembly output gets put into there
cd ~/all_gammaproteobacteria_data/assemblies

# get all the assemblies for our final list of organimsmms and then retrieve metrics like the GC content and total assembly length. This script is what actually downloads all the assemblies we need in zip format
python3 ~/scripts/assembly_and_metrics_grab.py \
	-i ~/all_gammaproteobacteria_data/output_data/filtered_gammaproteobacteria.txt \
	-d ~/all_gammaproteobacteria_data

# create the plots in R
Rscript ~/scripts/distribution_plotting.r ~/all_gammaproteobacteria_data/output_data all_species_metrics.txt

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
