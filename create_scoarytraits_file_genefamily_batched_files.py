# This script will create the traits file to be used in scoary for given gene families we are interested in. Because scoary is slow when you have many traits in one file
# This script will create traits files in batches of 5 for scoary to run.
# Before running this script you should can create a file with a list of gene families you want to find traits for

# import modules
import argparse
import os
import pandas as pd
import logging
import json
import re

my_vars = argparse.ArgumentParser(description='Creates scoary traits file with gene families in mind')

my_vars.add_argument('-r','--rtab',dest='rtab_file', type=str, required=True, help='give the absolute path to the gene presence Rtab file produced by panta')

# my_vars.add_argument('-c','--csv',dest='csv_file', type=str, required=True, help='give the absolute path to the gene presence csv file produced by roary ')
# open(my_files.csv_file,'r') as csv_file,
#my_vars.add_argument('-id','--input_directory',dest='rgi_parsing_output_directory', type=str, required=True, help='give the absolute path to the directory produced from the rgi_parsing_genefamily.py script which has the tsv and json')

my_vars.add_argument('-j','--json_file',dest='gf_json', type=str, required=True, help='give the absolute path to the json file produced from the rgi_parsing_genefamily.py script')

my_vars.add_argument('-gf','--gene_families',dest='gene_families',type=str, required=True, help='absolute path to file containing the gene families of interest. Each gene family should be on a separate line')

my_vars.add_argument('-od','--out_dir',dest='output_dir', type=str, required=True, help='give the absolute path to an output directory to put all the output files into')

my_vars.add_argument('--traits_per_file',dest='batch_size', type=int, required=False,default=5, help='number of traits to put in a single file. Default:5')


my_files = my_vars.parse_args()

# create log file directory if it doesnt exist
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# configure the logger we will be using
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/creating_gene_family_traits_file.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
logger = logging.getLogger(__name__)

# create directory for scoary output
if not os.path.isdir(f"{my_files.output_dir}"):
    os.makedirs(f"{my_files.output_dir}")
    logger.info(f"Creating new directory {my_files.output_dir}")

# function to take the different gene families and
def create_blocks_of_five(gene_families:list,group_size = my_files.batch_size):

    # dict to store different groups
    group_name_gene_families = {}

    # set a variable to have as the group number
    group_id = 1

    # list to store the 5 gene families
    group_members = []

    # now go through the whole list of gene families
    for i in range(len(gene_families)):

        # add families to a group if the size is less than how many i want in a traits file
        if len(group_members) < group_size:
            group_members.append(gene_families[i].rstrip('\n'))

        # now check if the group members have been fillied to the capacity of the group size or if we are at the final element in the list.
        # if we are at the final element it wouldve been added to the list above so this is a valid method
        if len(group_members) == group_size or i == len(gene_families)-1:

            # add the group members to the dictionary with the group id
            group_name = f"group_{group_id}"

            group_name_gene_families[group_name] = group_members

            # set the group members to be empty so another five can be taken
            group_members = []

            # update the group id
            group_id += 1

    return group_name_gene_families

# open the file with the full list of gene families and
with open(my_files.gene_families,'r') as list_of_gene_families,open(my_files.rtab_file,'r') as rtab_file,open(f"{my_files.output_dir}/trait_to_regular_name.json",'w') as traits_to_regular_json_file,open(my_files.gf_json,'r') as gf_json_file:

    # load in the data matrix for Rtab file
    rtab_df = pd.read_table(rtab_file)

    # right now the data frame has the different species as column names, we store that in a variable
    different_species = rtab_df.columns.to_list()

    # delete the 'gene' column which is the column header for all the different gene names
    del different_species[0]

    # put number of species we have into a varibale
    number_of_species = len(different_species)

    # log
    logger.info(f"Number of species: {number_of_species}")

    # create json variable for the conversion from the file name to the true gene family name without underscores
    traits_to_regular_json = {}

    # load in json data with different information about the species, the gene family, then the ARO and stuff from the card database. The inital rgi output
    input_json_data_gene_family = json.load(gf_json_file)

    # store the lines in a variable
    gene_families = list_of_gene_families.readlines()

    # break them up into groups of 5
    groups_of_gene_families:dict = create_blocks_of_five(gene_families=gene_families)

    # go through each group of 5
    for group in list(groups_of_gene_families.keys()):

        # get list of gene families to write to the output file
        families_to_write = groups_of_gene_families[group]

        # we technically dont need to create a file for them, we can just log what the gene group is
        logger.info(f"Acessing Gene Group: {group}, Gene Families: {groups_of_gene_families[group]}")
        # with open(f"{my_files.output_dir}/gene_family_{group}.txt",'w') as group_out_file:
            # group_out_file.writelines([f"{family}\n" for family in families_to_write])


        # first we open traits files for that unique group of genes
        with open(f"{my_files.output_dir}/scoary_traits_file_{group}.trait",'w') as output_file:

            # store gene family data in a variable
            query_gene_family_data = groups_of_gene_families[group]

            # remove weird characters from any of the names so that scoary wont have any errors, like '',?,[,],. etc. We are basically changing anything that isnt a number or letter to an _
            query_gene_family_data_removed_weird_characters = [re.sub(pattern=r'[^a-zA-Z0-9]',repl= '_', string=family.rstrip('\n')) for family in query_gene_family_data]

            # create the json with the traits file name as they key and the proper gene family name as the value
            # traits_to_regular_json = {query_gene_family_data_removed_weird_characters[i]:query_gene_family_data[i].rstrip() for i in range(len(query_gene_family_data))}
            for i in range(len(query_gene_family_data)):
                traits_to_regular_json[query_gene_family_data_removed_weird_characters[i]] = query_gene_family_data[i].rstrip('\n')

            # add an inital empty space to the traits file, and it has to be csv separated like the gene presence file
            output_file.write("Name,")

            # # load in the data matrix for Rtab file
            # rtab_df = pd.read_table(rtab_file)

            # # right now the data frame has the different species as column names, we store that in a variable
            # different_species = rtab_df.columns.to_list()

            # # delete the 'gene' column which is the column header for all the different gene names
            # del different_species[0]

            # # put number of species we have into a varibale
            # number_of_species = len(different_species)

            # # log
            # logger.info(f"Number of species: {number_of_species}")

            # write all the gene families as headers to the traits file
            output_file.writelines([f"{query_gene_family_data_removed_weird_characters[i]}\n" if i == len(query_gene_family_data_removed_weird_characters)-1 else f"{query_gene_family_data_removed_weird_characters[i]}," for i in range(len(query_gene_family_data_removed_weird_characters))])

            # store json data with different information about the species, the gene family, then the ARO and stuff from the card database
            # input_json_data_gene_family = json.load(gf_json_file)

            # now loop through each strain/species we have and get the presence or absence for each gene family
            for species in different_species:

                # dummy vars to be able to search the json for those two species that had their names truncated
                query_json_file = species
                if species == 'Avibacterium_sp':
                    query_json_file = 'Avibacterium_sp_20_132'
                if species == 'Methylotuvimicrobium_alcaliphilum':
                    query_json_file = 'Methylotuvimicrobium_alcaliphilum_20Z'

                # first write the species to the line
                output_file.write(f'{species},')

                # get the different gene families for the species from the input json
                initial_species_specific_gene_families = list(input_json_data_gene_family[query_json_file].keys())

                # some gene families rows will have multiple gene families on the same row separated by ; so this command both helps remove duplicates, and count the multiple gene families on the line as their own distinct family
                species_specific_gene_families = list(set([family.strip() for multiple_gene_families in initial_species_specific_gene_families for family in multiple_gene_families.split(';')]))

                # go through each gene family that we are interested in in the traits file
                for i in range(len(query_gene_family_data)):

                    # create variable for gene family. Have to rstrip because it has a new line at the end
                    family = query_gene_family_data[i].rstrip()

                    # check if the species has that gene family or not and create the traits file accordingly. present = 1, absent = 0
                    if family in species_specific_gene_families:
                        presence_absence_value = 1

                    elif family not in species_specific_gene_families:
                        presence_absence_value = 0

                    # write to the output file differently if it is on the last column, or if it is one of the values leeding up to the last column
                    if i == len(query_gene_family_data) -1:
                        output_file.write(f"{presence_absence_value}\n")

                    elif i != len(query_gene_family_data) -1:
                        output_file.write(f"{presence_absence_value},")

                    # write to log file
                    logger.info(f"Species:{species}\nCurrent gene family:{family}\nPresence/absences value:{presence_absence_value}")

    # dump the json that allows you to convert the file names to their true gene family ames
    json.dump(obj=traits_to_regular_json,fp=traits_to_regular_json_file,indent=4)
