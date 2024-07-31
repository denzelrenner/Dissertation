# This script will run busco on all of the assemblies we have filtered down to. This will then be plotted and ranked from the best busco score to the lowest in a different script

# import modules
import argparse
import os
import re

my_vars = argparse.ArgumentParser(description='Running busco on all of our asseblies')

my_vars.add_argument('-nd','--nucleotide_dir',dest='nucleotide_dir', type=str, required=True, help='give the absolute path to the directory containing all the nucleotide fasta files')

my_vars.add_argument('-od','--odir',dest='output_dir', type=str, required=True, help='give the absolute path to an output directory to put all the busco output into')

my_vars.add_argument('-l','--lindtst',dest='lineage_dataset', type=str, required=True, help='give the absolute path to the directory with the lineage dataset')

my_files = my_vars.parse_args()

# check if directories exist
# create output directory if it doesnt exist
if not os.path.isdir(my_files.output_dir):
    os.makedirs(my_files.output_dir)

# create log file directory if it doesnt exist 
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# store all nucloetide fasta files in a variable 
fasta_files = os.listdir(my_files.nucleotide_dir)

# open a log file to write important output
with open(f"{os.path.expanduser('~')}/log_files/running_busco.log",'w') as log_file:
    
    # loop through each fasta file
    for file in fasta_files:

        # store absolute path to file in a variable
        abs_path = f"{my_files.nucleotide_dir}/{file}"

        # get rid of anything past the species name. ie the year of the assembly submission
        #only_species_name = re.sub(pattern=fr'_20.*\.fna',repl=f'',string=file)
        only_species_name = file.replace('.fna','')[:-11]

        # run busco on the file 
        os.system(f"busco -i {abs_path} --lineage_dataset {my_files.lineage_dataset} -o {only_species_name} -m genome --out_path {my_files.output_dir} -c 90")

        log_file.write(f"Running busco on file {abs_path}, directory output name will be {only_species_name} stored inside another called {my_files.output_dir}. Lineage dataset is {my_files.lineage_dataset}\n")
