#!/bin/bash

#SBATCH --job-name=acquiring_data
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=55g
#SBATCH --time=24:00:00
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
	--fields accession,organism-name,assminfo-release-date,checkm-completeness > ~/all_gammaproteobacteria_data/output_data/all_gammaprotobacteria.txt
#datasets summary genome taxon 1236 --assembly-source 'RefSeq' --as-json-lines | dataformat tsv genome \
#        --fields accession,organism-name,assminfo-release-date,checkm-completeness > ~/all_gammaproteobacteria_data/output_data/all_gammaprotobacteria.txt

# filter the data to remove duplicate strains by taking the most recent assembly or if that fails then taking the best checkM completeness
python3 ~/scripts/filtering_gammaprotobacteria.py \
	-i ~/all_gammaproteobacteria_data/output_data/all_gammaprotobacteria.txt \
	-o ~/all_gammaproteobacteria_data/output_data/filtered_all_gammaprotobacteria.txt

# move into the assemblies directory so all the assembly output gets put into there
cd ~/all_gammaproteobacteria_data/assemblies

# get all the assemblies for our final list of organimsmms and then retrieve metrics like the GC content and total assembly length
python3 ~/scripts/assembly_and_metrics_grab.py \
	-i ~/all_gammaproteobacteria_data/output_data/filtered_all_gammaprotobacteria.txt \
	-d ~/all_gammaproteobacteria_data

# create the plots in R
Rscript ~/scripts/distribution_plotting.r ~/all_gammaproteobacteria_data/output_data all_species_metrics.txt

# deactivate conda env
conda deactivate

# get job id
echo "The Job ID for this job is: $SLURM_JOB_ID"
