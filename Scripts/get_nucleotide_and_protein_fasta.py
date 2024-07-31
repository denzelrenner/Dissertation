# import subprocess,os


# subprocess.run(['ls','-l',f"{os.path.expanduser('~')}/Desktop/Final_Project/assembly_testing"])

# This script will unzip all of the assemblies 

# import modules
import argparse
import os

my_vars = argparse.ArgumentParser(description='Producing nucleotide and protein fasta files for filtered gammaproteobacteria')

my_vars.add_argument('-i','--input',dest='input_dir', type=str, required=True, help='give the absolute path to the directory containing all the zip files')
my_vars.add_argument('-d1','--dir_one',dest='output_dir_one', type=str, required=True, help='give the absolute path to an output directory to put all the assemblies fasta files into')
my_vars.add_argument('-d2','--dir_two',dest='output_dir_two', type=str, required=True, help='give the absolute path to an output directory to put all the protein fasta files into')

my_files = my_vars.parse_args()

# list all the zip files in the assemblies directory
all_assemblies = os.listdir(my_files.input_dir)

# check if there is a directory for log files in the home directory and if not create one 
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# check if the output directories already exist or not
if not os.path.isdir(f"{my_files.output_dir_one}"):
    os.makedirs(f"{my_files.output_dir_one}")

if not os.path.isdir(f"{my_files.output_dir_two}"):
    os.makedirs(f"{my_files.output_dir_two}")

# open the log file
with open(f"{os.path.expanduser('~')}/log_files/grabbing_protein_and_nucleotide_fasta.log",'w') as log_file:

    # loop through each file
    for file in all_assemblies:
        
        # only consider zip files. There should only be zip files in the input directory but just in case there isnt 
        if file.endswith('.zip'):
        
            #get absolute file path 
            abs_path = f"{my_files.input_dir}/{file}"

            # get the name of the species without the dataset and .zip file
            species = file.replace('_dataset.zip','')

            # create a new directory for that species to put all the zip folder contents into
            os.system(f"mkdir {my_files.input_dir}/{species}")

            # now unzip the file and send the contents to the directory you just created
            os.system(f"unzip {my_files.input_dir}/{file} -d {my_files.input_dir}/{species}")

            # move the nucleotide fasta file to the directory created for nucleotide fasta files
            os.system(f"cp {my_files.input_dir}/{species}/ncbi_dataset/data/GCF*/GCF*.fna {my_files.output_dir_one}/{species}.fna")

            # now move the protein fasta file to the directory that will be used for protein fasta files
            os.system(f"cp {my_files.input_dir}/{species}/ncbi_dataset/data/GCF*/*.faa {my_files.output_dir_two}/{species}.faa")

            # write to the log file
            log_file.write(f"Created directory called {species}. Unzipping file called {file} and copying fasta file called {species}.fna to {my_files.output_dir_one}. The absolute path to this file is {abs_path}\n")
            log_file.write(f"Created directory called {species}. Unzipping file called {file} and copying fasta file called {species}.faa to {my_files.output_dir_two}. The absolute path to this file is {abs_path}\n")



