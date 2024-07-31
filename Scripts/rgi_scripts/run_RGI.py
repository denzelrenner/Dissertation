# This script runs RGI on all of our protein fastas.
# This script will also use the txt files produced by RGI to identify gene families in all of our genomes.
# a json with all that information will be produced as well as a tsv that can be plotted in R

# import modules
import argparse
import os
import logging
import pandas as pd
import json
import shutil
import sys


vars = argparse.ArgumentParser(description='Running RGI main on all 1348 protein fastas and Identifying gene families for different species')

vars.add_argument('-id','--input',dest='input_dir', type=str, required=True, help='give the absolute path to the input directory containing all the protein fasta files')

vars.add_argument('-od','--output',dest='output_dir', type=str, required=True, help='give the absolute path to the output directory to host all rgi related output')

vars.add_argument('-t','--threads',dest='cores', type=int, required=True, help='number of threads for RGI to use')

vars.add_argument('-mj','--master_json',dest='master_json', type=str, required=True, help='Give the path to the file with all the json data for all the different assemblies')

vars.add_argument('-x','--exclude',dest='exclude_hits', default=[], required=False, nargs='+', help='Specify which hits to ignore when making the gene families. Enter any combination of "Perfect","Strict","Loose" if flag is used. Example "-x Loose Strict"')

vars.add_argument('--remove_temp',dest='remove_temp_files', default=False, action='store_true', help='calling this flag will remove the txt files created by RGI because they take up space and are only needed for troubleshooting')

args = vars.parse_args()

# check if there is a directory for log files in the home directory and if not create one
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# configure the logger we will be using
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/run_RGI_1173_genomes_addedstrip.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
logger = logging.getLogger(__name__)

# list all the protein fasta files in the protein fasta directory
all_protein_fasta= os.listdir(args.input_dir)

# log
logger.info(f'Getting all {len(all_protein_fasta)} protein fasta files from input directory ')

# log
logger.info(f"Hits to exclude are set to {args.exclude_hits}")

# check if their output directories already exist or not
if not os.path.isdir(f"{args.output_dir}"):
    os.makedirs(f"{args.output_dir}")
    logger.info(f"Creating output directory {args.output_dir}")

# set variable for json output data
output_json_data = {}

# variable for number of genomes processed
number_of_genomes_processed = 0

# open tsv file that will be used to plot in R, open a file that will be used to write json data to, open the master json file that has information on all assemblies like the GCF accessioncode
with open(f"{args.output_dir}/input_data_to_plot_gene_families.tsv",'w') as input_file_R_plotting, open(f"{args.output_dir}/all_gene_family_data.json",'w') as json_file, open(args.master_json,'r') as file_for_master_json:

    # log
    logger.info(f"Creating input file to plot gene families in R: {args.output_dir}/input_data_to_plot_gene_families.tsv")

    # create headers for R file
    input_file_R_plotting.write(f"file name\ttrue species name\tgene family\n")

    # load the json data from the master json
    master_json = json.load(file_for_master_json)

    # loop through each file
    for file in all_protein_fasta:

        # get absolute file path
        abs_path_protein_fasta = f"{args.input_dir}/{file}"

        # get the name of the species without the dataset and .zip file. Gets rid of the
        species = file[:-15]

        # Get true species name without the _ and stuff by looking into the master json
        true_species_name = master_json[species]['Organism Name']

        # create a new directory for that species to put all the zip folder contents into
        if not os.path.isdir(f"{args.output_dir}/{species}"):
            os.makedirs(f"{args.output_dir}/{species}")
            logger.info(f"Creating species specific directory {args.output_dir}/{species}")

        elif os.path.isdir(f"{args.output_dir}/{species}"):
            logger.error(f"Crucial error in creating directory {args.output_dir}/{species}. It already exists which means there might be a duplicated species name. Ending program")
            sys.exit(1)

        # create variable for the directory we want to remove
        dir_to_remove = f"{args.output_dir}/{species}"

        # run rgi on the protein fasta
        os.system(f"rgi main --input_sequence {abs_path_protein_fasta} --output_file {args.output_dir}/{species}/{species} --local --clean -t protein -n {args.cores}")

        # check if the txt file with all the information exists
        if os.path.isfile(f"{args.output_dir}/{species}/{species}.txt"):
            logger.info(f"Current directory: {args.output_dir}/{species}\n Found file: {args.output_dir}/{species}/{species}.txt")

        elif not os.path.isfile(f"{args.output_dir}/{species}/{species}.txt"):
            logger.error(f"Critical Error. Could not find the gene family txt file in {args.output_dir}/{species}. Ending program")
            sys.exit(1)

        # get absolute path for the input txt file
        abs_path = f"{args.output_dir}/{species}/{species}.txt"

        # create a json key/entry for a given species. we will still use the file name as the key because it allows it to be used by different files
        output_json_data[species] = {}

        # load txt file produced from RGI into a pandas dataframe
        df = pd.read_csv(filepath_or_buffer=abs_path,header=0,sep='\t')

        # store all column names in a variable
        column_headers = df.columns.to_list()

        # put all the unique values for gene families into a list. Do not separate the AMRs in the colon on the same line here because then the exact row cant be accessed
        unique_gene_families = df['AMR Gene Family'].unique()

        # log
        logger.info(f"Found directory: {species}\nSpecies name is: {species}\ngene families for species are: {unique_gene_families}")

        # go though each unique value and put that in a json
        for family in unique_gene_families:

            # for a given gene family get all the rows with that gene family
            gene_families_found = df[df['AMR Gene Family']==family]

            # there may be multiple lines that have a certain gene family. So for a given gene family we create a list of rows and put a dictionary of each row into the list
            # gene_families_found.to_dict(orient='records') does exactly that
            # rows_for_gene_family = gene_families_found.to_dict(orient='records')

            # loop through the rows and check if the Cut_Off is something that the user has excluded like loose,strict, or perfect. Using this method, gene families that have been excluded just have a [] and no dict inside of it which is true because their cutoff was part of the values the user decided to exclude
            rows_for_gene_family = [dict_of_rows for dict_of_rows in gene_families_found.to_dict(orient='records') if dict_of_rows['Cut_Off'] not in args.exclude_hits]

            # so by not including loose hits the list above is returned as empty rather than having the hit for gene , and this gets added to the json and is considered
            # as one of the gene families for a specific species. The line below ensure gene famileis are considered only if the list of dicts has a value inside of it
            if rows_for_gene_family:

                # some family lines will have a ; in them so we can split the gene family string by ; to create a list of the different families
                # even if the family doesnt have ; it will still put the gene family in a list so we can loop through it
                for separated_family in family.split(';'):

                    # put the list of dictionarys/ the rows for each gene family as the value for the gene family key
                    output_json_data[species][separated_family.strip()] = rows_for_gene_family

                    # write to the input file for the R script that will plot the different gene families that we have in our data.
                    # so we will need the species and the gene family you found for it, technically we could just loop through the final json and tally but this is probably easier
                    input_file_R_plotting.write(f"{species}\t{true_species_name}\t{separated_family.strip()}\n")

                    # log
                    logger.info(f"Writing {true_species_name} and {separated_family.strip()} to input file to plot gene families in R: {args.output_dir}/input_data_to_plot_gene_families.tsv")

            # log info if the list for the gene family is empty. will only be relevant if loose hits are included
            elif not rows_for_gene_family:

                logger.error(f"Species is {species}.Could not find any hits in dataframe for {family}. Not writing to file")

        # we are done getting the gene family information out of it so now we can remove the directory
        if args.remove_temp_files:

            # check if the path is to a directory and it does not end with my user, meaning it will not remove
            if not os.path.isdir(dir_to_remove) or dir_to_remove.split('/')[-1] in os.path.expanduser('~').split('/'):
                logger.critical("Stopping program. A main directory like the home directory is being set to be removed or the dir to remove does not exist")
                sys.exit(1)

            # if the path is to an actual dir, and it doesnt end with something critical like the users home profile
            elif os.path.isdir(dir_to_remove) and not dir_to_remove.split('/')[-1] in os.path.expanduser('~').split('/'):
                shutil.rmtree(dir_to_remove)
                logger.info(f"Directory {dir_to_remove} does not end with home directory name so we are safely removing it")

        # add to the number of genomes processed and log it
        number_of_genomes_processed += 1
        logger.info(f"Current file = {file}. Number of genomes processed = {number_of_genomes_processed}")

    # put the dictionary data for the gene families into a json file
    json.dump(obj=output_json_data,indent=4,fp=json_file)

    # log
    logger.info(f"Dumping all species and gene family data into json file: {args.output_dir}/all_gene_family_data.json.")

# now we want to create an output file with all the gene families in our data set
with open(f"{args.output_dir}/unique_gene_families.txt",'w') as output_file_gene_families:

    # load the data that will be plotted into R into a data frame
    df = pd.read_csv(filepath_or_buffer=f"{args.output_dir}/input_data_to_plot_gene_families.tsv",sep='\t',header=0)

    # get the unique gene families
    gene_families = df['gene family'].unique().tolist()

    # loop through the list on separate the hits with multiple gene families, then filter the list again to only contain unique values.
    # using lstrip and rstrip to remove the spaces at the start and end of the hits.
    # when there are multiple gene families it is in the file as genefamA; genefamB; genefamC so we split by ; and remove the spaces at the beginning og genefamB and C
    filtered_gene_families = list(set([individual_family.lstrip().rstrip() for families in gene_families for individual_family in families.split(';')]))

    # write to the output file
    output_file_gene_families.writelines([f"{unique_family}\n" for unique_family in filtered_gene_families])

    # log
    logger.info(f"Created file {args.output_dir}/unique_gene_families.txt.\nWriting unique gene families to file....")
    logger.info(f"END.")
