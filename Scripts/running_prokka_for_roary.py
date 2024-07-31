# This script is going to run prokka on each of our individual assemblies and then we will use the gff output of the different assemblies to run roary

# import modules
import argparse
import os
import re

my_vars = argparse.ArgumentParser(description='Taking an input and output file')

my_vars.add_argument('-nd','--nucleotide_dir',dest='nucleotide_dir', type=str, required=True, help='give the absolute path to the input directory with nucleotide fasta files')

my_vars.add_argument('-pd','--protein_dir',dest='protein_dir', type=str, required=True, help='give the absolute path to the input directory with all the protein fasta files')

my_vars.add_argument('-if','--input_file',dest='input_file', type=str, required=True, help='give the absolute path to the file containing the deduplicated set of gammoproteobacteria')

my_vars.add_argument('-o','--out',dest='output_dir', type=str, metavar='output dir', required=True, help='give the absolute path to the output directory')

my_files = my_vars.parse_args()

# create output directory if it doesnt exist
if not os.path.isdir(my_files.output_dir):
    os.makedirs(my_files.output_dir)

# make a directory to store all the species gffs
if not os.path.isdir(f"{my_files.output_dir}/all_species_gff"):
    os.makedirs(f"{my_files.output_dir}/all_species_gff")

# create log file directory if it doesnt exist 
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# open the input file
with open(my_files.input_file,'r') as input_file, open(f"{os.path.expanduser('~')}/log_files/run_prokka_for_roary.log",'w') as log_file:

    # create a string that can be searched through using re for all the different species in the final list. and replace any double underscores with a regular one
    species_data = input_file.read().replace('__','_')

    # load in the files from the input directory
    fasta_files = os.listdir(my_files.nucleotide_dir)

    # loop through each fasta file
    for file in fasta_files:

        # get rid of anything past the species name. ie the year of the assembly submission
        only_species_name = file.replace('.fna','')[:-11]
        #only_species_name = re.sub(pattern=fr'_20.*\.fna',repl=f'',string=file)

        # get the name of the protein fasta file
        #protein_fasta = re.sub(pattern=fr'\.fna',repl=f'.faa',string=file)
        
        # get the name of the protein fasta file
        protein_fasta = file.replace('.fna','.faa')

        # check if you can find the species in the input file and if you can the run prokka using that file 
        if re.search(pattern=only_species_name,string=species_data):

            # make a directory for output for that specific species
            os.makedirs(f"{my_files.output_dir}/{only_species_name}")

            log_file.write(f"Creating directory {my_files.output_dir}/{only_species_name}\n")
            
            # run prokka
            os.system(f"prokka --outdir {my_files.output_dir}/{only_species_name} --prefix {only_species_name} --proteins {my_files.protein_dir}/{protein_fasta} --cpus 95 --kingdom Bacteria --force {my_files.nucleotide_dir}/{file}")
            
            # write to log file
            log_file.write(f"the species {only_species_name} is part of the final selection of gammaproteobacteria and its protein fasta is {protein_fasta} and nucleotide fasta is {file}\n")

            # move the gff to the directory for gffs
            os.system(f"cp {my_files.output_dir}/{only_species_name}/*.gff {my_files.output_dir}/all_species_gff")

         
