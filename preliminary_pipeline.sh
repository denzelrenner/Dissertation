#!/bin/bash

# This script runs the pipeline for grabbing the assemblies we need 

# SOURCE HOME PROFILE

# activate the conda environment 
conda activate ncbi_datasets

# grab all gammaproteobacteria hits from the ncbi website. 1236 is the taxon id for gammaproteobacteria
datasets summary genome taxon 1236 --assembly-level complete --assembly-source 'RefSeq' --as-json-lines | dataformat tsv genome --fields accession,organism-name > refseq_complete_gammaprotobacteria.txt

# filter the data to remove duplicate strains by taking the most recent assembly or if that fails then taking the best checkM completeness
python3 ../filtering_gammaprotobacteria.py -i refseq_complete_gammaprotobacteria.txt -o mouse.txt

# get all the assemblies for our final list of organimsmms and then retrieve metrics like the GC content and total assembly length
python3 ../test_assembly_grab.py -i ../test_file.txt -d {output_directory}

# plot the data in R to remove outliers
{R_script}