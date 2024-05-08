# This script will filter the multiple assembly entries based on one of the criteria we have set. So it looks for which assembly has the latest submission date
# and if that doesnt work out then it takes the most coverage. 

# Think about doing a json dump instead? Might make things more difficult down the line but who knows
# Think about creating a class for each species, and immediately filter it out if it doesnt meet some threshold we set

# Example run command: python3 ../filtering_gammaprotobacteria.py -i refseq_complete_gammaprotobacteria.txt -o mouse.txt


# import modules
import argparse
import os
import re

my_vars = argparse.ArgumentParser(description='Taking an input and output file')

my_vars.add_argument('-i','--input',dest='input', type=str, required=True, help='give the absolute path to the file containing the species you are interested in and their accessions')

my_vars.add_argument('-o','--out',dest='output_file', type=str, metavar='output file', required=True,
                        help='give the absolute path to the output file name with .out at the end')

my_vars.add_argument('-d','--odir',dest='output_dir', type=str, required=False, help='give an output directory')

my_files = my_vars.parse_args()


# use index to find where in the list a certain value is
# create a queue, if you arrive on a line, before adding it check whether the next line species name is the same, if not then add it and empty the queue.
# If it is then change what is in the queue depending on if what is in the queue has higher genome completeness for example

# open input file
with open(my_files.input,'r') as input_data, open(my_files.output_file,'w') as filtered_file:

    # set a default empty value for the queue. The queue should have a tuple of the species name and the metric we are using to choose which assembly to use
    queue = []

    # put all data in a variable then go back to the top of the file
    full_data = input_data.read().replace('[','').replace(']','')

    # print(full_data[0:300])
    
    input_data.seek(0)
   
    # set arbitrary value for current species
    current_species = ''

    # go through each line of the data with all the species in it
    for line in input_data.readlines()[1:]:

        # print(f"the line is {line}")
        # get the species being looked at
        full_species_name = line.split('\t')[1]

        # now we only want genus and species, we dont care about the strain
        genus_species = ' '.join(full_species_name.split()[0:2])

        # print(genus_species)
        # If we have already done summary calculations for that species, move on to the next line
        if genus_species == current_species:
            continue
        
        # otherwise if the species in the line isnt in the current species then do the re thing
        elif genus_species != current_species:
            
            # add the organism to the current species variable
            current_species = genus_species

            # set variable for what will be the representative species
            representative_species = ''

            # find the organism in the data
            all_organisms = dict(enumerate(re.findall(rf"GC.*{genus_species}.*\n",full_data)))

            # print(all_organisms)
            # set arbitrary default value to contest against
            latest_assembly_best_covg = [-9,(1999,1,1),55]

            # loop through all the re search hits. So all the different strains of a given species
            for tracker,data in all_organisms.items():
                
                # to get the year month day first split the line by tabs then split the date by dashes. its fine to have 08 string because python int will make it 8
                year,month,day = data.split('\t')[2].split('-')

                coverage = float(data.split('\t')[3])

                # now do checks to see which is the latest assembly
                if int(year) > latest_assembly_best_covg[1][0]:
                    
                    # change the latest assembly to be the one we just looked at
                    latest_assembly_best_covg = [tracker,(int(year),int(month),int(day)),coverage]

                # if the year in the strain/organism you are looking at is lower than our current best then forget about it entirely and move on
                elif int(year) < latest_assembly_best_covg[1][0]:

                    continue

                # check if the year is the same
                elif int(year) == latest_assembly_best_covg[1][0]:

                    # check if the month is greater
                    if int(month) > latest_assembly_best_covg[1][1]:

                        # change the latest assembly to be the one we just looked at
                        latest_assembly_best_covg = [tracker,(int(year),int(month),int(day)),coverage]

                    # if the month is lower theres no need to consider this organism again
                    elif int(month) < latest_assembly_best_covg[1][1]:

                        continue
                    
                    # check if the months are the same
                    elif int(month) == latest_assembly_best_covg[1][1]:

                        # if they are check the day. then if day is the same use coverage
                        if int(day) > latest_assembly_best_covg[1][2]:

                            # change the latest assembly to be the one we just looked at
                            latest_assembly_best_covg = [tracker,(int(year),int(month),int(day)),coverage]

                        # if the day is the earlier than what we currently have just move on
                        elif int(day) < latest_assembly_best_covg[1][2]:

                            continue

                        # if the days are the same then resort to using coverage
                        elif int(day) == latest_assembly_best_covg[1][2]:

                            # check if coverage is better with the current organism
                            if coverage >= latest_assembly_best_covg[2]:

                                # replace the current best assembly
                                latest_assembly_best_covg = [tracker,(int(year),int(month),int(day)),coverage]

            # write to the output file using the tracker (dictionary id) for the assembly that passed all the criteria i gave
            # print(f"the best assembly is {all_organisms[latest_assembly_best_covg[0]]}")
            filtered_file.write(all_organisms[latest_assembly_best_covg[0]])
