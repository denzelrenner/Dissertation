# This script will take in the distance matrix of ANI values, coverage, and alignment length produced by the ANI calculating script, and then it will remove species from the dataset that are very similar to each other

# # import modules
import argparse
import os
import pandas as pd

my_vars = argparse.ArgumentParser(description='Deduplicates the dataset based on ANI values')

my_vars.add_argument('--ani',dest='ani_input', type=str, required=True, help='give the absolute path to the file containing the distance matrix of ANI values')

my_vars.add_argument('-al','--align',dest='align_input', type=str, required=True, help='give the absolute path to the file containing the distance matrix of allignment lengths')

my_vars.add_argument('-cov','--coverage',dest='coverage_input', type=str, required=True, help='give the absolute path to the file containing the distance matrix of coverage percentage')

my_vars.add_argument('-o','--out',dest='output_file', type=str, metavar='output file', required=True,help='give the absolute path to the output file')

my_files = my_vars.parse_args()

# create log file directory if it doesnt exist 
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# create function which will return how many times a certain genus appears in the dataset. we need this because originally species which were the only ones in the genus werent 
# being included in the final calculations
def is_single_species(all_species:list,genus_of_interest:str):

    # set a counter 
    counter = 0

    # loop through each of the species we have 
    for species in all_species:

        # check if any of the species start with the genus of interest
        if species.startswith(genus_of_interest):

            counter += 1
    
    # if the counter is greater than 0 then return multiple species because it means there are others
    if counter > 1:
        return False
    
    # if the counter remains at 0 then return single species because it is the only one with that genus 
    elif counter == 1:
        return True
    
# open the input file and a log file
with open(my_files.ani_input,'r') as ani_input_file, open(my_files.output_file,'w') as output_file, open(f"{os.path.expanduser('~')}/log_files/deduplicating_dataset.log",'w') as log_file, open(my_files.align_input,'r') as alignment_input, open(my_files.coverage_input,'r') as cov_input:

    # store distance matrix of ANI values in a variable
    ANI_df = pd.read_table(ani_input_file)
    
    # store distance matrix of alignment lengths in a variable
    alignment_df = pd.read_table(alignment_input)

    # store distance matrix of alignment lengths in a variable
    coverage_df = pd.read_table(cov_input)

    # sense check data in distance matrix
    # print(df.loc[0,'secondary_endosymbiont_of_Trabutina_mannipara_2016-07-16'])

    # store all the different species names into a variable
    column_names = ANI_df.columns.to_list()
    
    # remove the first entry because it is 'unnames:0' and not one of our species
    del column_names[0]

    # put number of species we have into a varibale 
    number_of_species = len(column_names)

    # create a dictionary to host the species and the unique id/index for it 
    rownum_species = dict(enumerate(column_names))

    # store all organisms that pass the ANI checks in this variable
    final_species = []

    # store all the species that were removed for being duplicates 
    detained_species = []
    
    # store all species that were not put into detained or final species in this variable
    # they exist because within a genus there can be ANI values above 0.95 and in that case everything is detained except the last species because it cant get past the conditional statements
    error_species = []

    # now loop through each pairwise comparison and check if something is a reject or not 
    for query_species in column_names:
    
        # each species is assigned a different index so we can access that in the dictionary of species_rownum to get the correct species
        for row_number in range(number_of_species):
            
            # get the ani value
            ANI = float(ANI_df.loc[row_number,query_species])

            # get the alignment length 
            alignment_length = float(alignment_df.loc[row_number,query_species])

            # get the coverage
            coverage = float(coverage_df.loc[row_number,query_species])

            # get the name of the subject species
            subject_species = rownum_species[row_number]

            # get the genus ONLY for the query
            current_genus = query_species.split('_')[0]

            # do a check on if the current species is the only species in its genus, if it is then just add the query and break because theres no other species within the genus to compare it to
            if is_single_species(all_species=column_names,genus_of_interest=current_genus):

                final_species.append(query_species)

                log_file.write(f"This species {query_species} is the only one of its kind so it is being added to the final species list\n")
               
                break
           
            # check if the subject species is part of the detained species. In that case we shouldnt be doing comparisons with that query species anymore because it has been removed from the dataset
            # we dont need to worry about the query species here because we have already passed that species which had a high ANI similary to another subject
            # also check if we are comparing a query against itself. If we are then theres no need to consider the ANI because it will be 1
            # now if the query is being compared against a species that is not from the same genus then we just move on because this tool has issues when comparing diverse species
            elif query_species == subject_species or subject_species in detained_species or not subject_species.startswith(current_genus):
                
                continue
            
            # check whether the similarity threshold has been passed and the alignment length is sufficient to be making a conclusion about similarity added 0
            elif ANI > 0.95 and alignment_length > 1000000:
                
                # detain the current query species because it is too similar to another species
                detained_species.append(query_species)

                # remove the query from the list of final species if we have previously added it 
                if query_species in final_species:
                    
                    final_species.remove(query_species)

                log_file.write(f'Being detained... The query is {query_species} and subject is {subject_species}. The allignment length for these two is {alignment_length},The coverage is {coverage*100}. The ANI is {ANI}. The row index is {row_number}\n')

                break
            
            # check if the ANI is lower than our threshold.
            # before i also had the expression 'and query_species not in detained_species' but this doesnt make sense because when the query is added to the detained species we never see it again
            # also there is an initial if check which doesnt compare ANI with detained subject species 
            # a new query species can never be detained because it is only query species that get put into the detained list, not subjects
            elif ANI < 0.95 :

                # only add the query to the final species list if it is not currently inside there already. prevents duplicates
                if not query_species in final_species:

                    final_species.append(query_species)

                    log_file.write(f"Adding to final species... query: {query_species},\n subject: {subject_species} \n and ANI:{ANI}\n and alignment length {alignment_length}\n and coverage {coverage}\n")
            else:
                # when the allignment or coverage is very low write it to the log file 
                log_file.write(f"This is an exception with query: {query_species},\n subject: {subject_species} \n and ANI:{ANI}\n and alignment length {alignment_length}\n and coverage {coverage}\n")

    # go through all gammaproteobacteria we have and if they havent been detained or put in the final list we will save that to the error species variable
    for species in column_names:

        if species not in final_species and species not in detained_species:

            error_species.append(species)

            final_species.append(species)

            log_file.write(f"The species {species} was causing problems but has been added to the final species\n")


    # go through each of the final species and write their names to the output file
    for species in final_species:

        output_file.write(f"{species}\n")

    # write information to the log file because stdout and err can get cluttered sometimes 
    log_file.write(f'The number of species we are left with is {len(final_species)} and the number of detained species is {len(detained_species)}.\n The number of error species was {len(error_species)}\n and they were {error_species}')


