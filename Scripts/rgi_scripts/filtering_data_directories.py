# This script is going to take a file with different gammaproteobacteria on each line. It will also take a directory with information like protein fastas.
# It will look for the protein fasta for each of those gammaproteobacteria and copy them to a new output directory

# import modules
import argparse
import os
import re
import shutil
import logging

vars = argparse.ArgumentParser(description='Creates a new directory with files for a given list of species')

vars.add_argument('-id','--input_directory',dest='input_dir', type=str, required=True, help='give the absolute path to the input directory')

vars.add_argument('--input_file',dest='input_gammaproteobacteria', type=str, required=True, help='give the absolute path to the file containing a list of gammaproteobacteria on each line')

vars.add_argument('-od','--output_directory',dest='output_dir', type=str, required=True, help='give the absolute path to the output directory')

args = vars.parse_args()

# check if there is a directory for log files in the home directory and if not create one
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# configure the logger we will be using
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/filtering_data_directories.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
logger = logging.getLogger(__name__)

# create output directory if it doesnt exist
if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)

# open the input file
with open(args.input_gammaproteobacteria,'r') as input_file:

    # create a string that can be searched through using re for all the different species in the final list. and replace any double underscores with a regular one
    species_data = input_file.read().replace('__','_')

    # load in the files from the input directory
    data_files = os.listdir(args.input_dir)

    # loop through each file
    for file in data_files:

        # get rid of anything past the species name. ie the year of the assembly submission
        species = file[:-15]

        # get absolute path
        abs_path = f"{args.input_dir}/{file}"

        # path to send the file to
        destination_path = f"{args.output_dir}/{file}"

        # check if you can find the species in the input file and if you can then copy the file over to the new directory
        if re.search(pattern=species,string=species_data):

            # copy the fasta or whatever other data file for that species to the output directory
            shutil.copy2(src=abs_path,dst=destination_path)

            # log
            logger.info(f"Species: {species}, found data file: {file}, Copying to output directory {args.output_dir}")
