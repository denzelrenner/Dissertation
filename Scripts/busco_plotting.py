# This script will take the different output from running busco and plot the metrics for our report
# Produces a file called all_species_busco_metrics in the output directory that is specified

# import modules
import argparse
import os
import json
import shutil

my_vars = argparse.ArgumentParser(description='Plotting the data from running busco on all of our assemblies')

my_vars.add_argument('-bd','--busco_dir',dest='busco_directory', type=str, required=True, help='give the absolute path to the directory with all the different busco output')

my_vars.add_argument('-od','--out_dir',dest='output_dir', type=str, required=True, help='give the absolute path to the directory to store the plots and any relevant files in')

my_vars.add_argument('-c','--cutoff',dest='completeness_cutoff', type=float, required=True, help='give the cutoff you want to use for completeness [i.e 91.2] if all assemblies with a busco completeness of less than 91.2 should not be considered ')

my_vars.add_argument('-ani',dest='ani_filtered_file', type=str, required=True, help='path to the deduplicated file produced by using the ANI')

my_files = my_vars.parse_args()

# store all nucloetide fasta files in a variable. important to do this before making a new directories so they are not counted as one of the files
busco_files = os.listdir(my_files.busco_directory)

# check if directories exist
# create log file directory if it doesnt exist 
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# make output directory to put plots and busco summaries into
if not os.path.isdir(f"{my_files.output_dir}"):
    os.makedirs(f"{my_files.output_dir}")

# make temporary output directory for json files from busco 
if not os.path.isdir(f"{my_files.busco_directory}/json_files"):
    os.makedirs(f"{my_files.busco_directory}/json_files")


## Copying busco txt and json files to different directories. And extracting info about the busco scores
# open a log file to write important output
with open(f"{os.path.expanduser('~')}/log_files/plotting_busco_output2.log",'w') as log_file, open(f"{my_files.output_dir}/all_species_busco_metrics.tsv",'w') as output_file:
    
    # write headers to the busco metrics tsv file
    output_file.write(f"Organism\tComplete\tSingle_Copy\tMulti_Copy\tFragmented\tMissing\tn_markers\tdomain\tNumber_of_Scaffolds\tNumber_of_contigs\tTotal_Length\tPercent_gaps\tScaffold_N50\tContigs_N50\n")
    
    # loop through each busco output file/directory
    for file in busco_files:

        # write the name of the current species/organism to the output file
        output_file.write(f"{file}\t")

        # store absolute path to file in a variable
        abs_path = f"{my_files.busco_directory}/{file}"

        # copy the short summary txt file to the output directory which will be used as input for the plotting tool
        os.system(f"cp {abs_path}/short_summary*.txt {my_files.output_dir}")

         # write to log file to track which short summary txt files have been copied over
        log_file.write(f"Copied the short summary txt file for {file} to {my_files.output_dir}\n")
        
        # copy the json files to their own directory and change the names 
        os.system(f"cp {abs_path}/short_summary*.json {my_files.busco_directory}/json_files/{file}.json")

        # open the json file
        temp_file = open(f"{my_files.busco_directory}/json_files/{file}.json",'r')

        # load the json file 
        metrics_from_json = json.load(temp_file)

        # write each of the metrics to the tsv file except for the one line summary one
        for metric,value in metrics_from_json['results'].items():

            # separate values by tabs if it is not the final value to be output
            if metric != 'one_line_summary' and metric != 'Contigs N50' and metric != 'Percent gaps':

                output_file.write(f"{value}\t")

            # remove the % from percentage gap output
            elif metric == 'Percent gaps':

                temp_value = value.replace('%','')

                output_file.write(f"{temp_value}\t")

            # write a new line if it is the last metric
            elif metric == 'Contigs N50':

                output_file.write(f"{value}\n")
        
        # close the temp file
        temp_file.close()

        # check whether json file we just used exists and if it does, remove it
        if os.path.isfile(f"{my_files.busco_directory}/json_files/{file}.json"):

            os.remove(f"{my_files.busco_directory}/json_files/{file}.json")

            log_file.write(f"Found file {my_files.busco_directory}/json_files/{file}.json. Removing...\n")

        else:

            log_file.write(f"Error removing {my_files.busco_directory}/json_files/{file}.json. Not a real file.\n")
    
       

    # remove empty directory that was used to store json files
    # first check if directory actually exists 
    if os.path.isdir(f"{my_files.busco_directory}/json_files"):

        # exception handling in case the directory is not empty so os.rmdir cannot remove it
        try:
            os.rmdir(f"{my_files.busco_directory}/json_files")
            log_file.write(f"Removing empty directory {my_files.busco_directory}/json_files")

        except OSError as error:

            log_file.write(f"Cannot remove the directory {my_files.busco_directory}/json_files because it is not an empty directory")

    # if the directory doesnt exist write to log file
    else:

        log_file.write(f"Cannot remove the directory {my_files.busco_directory}/json_files because it doesnt exist")

# now run the generate plot command to get plots for the busco output
os.system(f"generate_plot.py -wd {my_files.output_dir}")

# sort the metrics file produced by the completeness
os.system(f"sort -k 2 -n {my_files.output_dir}/all_species_busco_metrics.tsv > {my_files.output_dir}/sorted_by_complete_all_species_busco_metrics.tsv")

# filter final gammaproteobacteria based on busco scores
# open the sorted busco metrics file and the output file for all the species that passed the busco completeness cutoff
with open(f"{my_files.output_dir}/sorted_by_complete_all_species_busco_metrics.tsv",'r') as busco_stats_file, open(f"{my_files.output_dir}/busco_filtered_gammaproteobacteria.txt",'w') as busco_filtered_file, open(f"{os.path.expanduser('~')}/log_files/plotting_busco_output2.log",'a') as log_file, open(f"{my_files.ani_filtered_file}",'r') as ani_file:

    # store data from ani file in a variable
    ani_file_data = ani_file.read()

    # store busco species metrics data in a variable
    busco_stats = busco_stats_file.readlines()[1:]

    # loop through all stats and write an organism to the output file if it has a higher busco completeness than the cutoff
    for organism in busco_stats:

        # get the different data for the organism
        name,completeness = organism.split('\t')[0],float(organism.split('\t')[1])

        # check if the name from the all species metric file is in the deduplicated dataset and hasnt been removed. alternative is to do this at the prokka stage but for now we do it here
        if name not in ani_file_data:

            log_file.write(f"The species {name}, is not in the deduplicated dataset so it is not being considered\n")

            continue
        
        # check if it has a higher cut off value
        elif completeness >= my_files.completeness_cutoff:

            # write to output file
            busco_filtered_file.write(f"{name}\n")

            # write to log file
            log_file.write(f"The species {name} has a busco completeness of {completeness} and so it is being written to {my_files.output_dir}/busco_filtered_gammaproteobacteria.txt\n")

        elif completeness < my_files.completeness_cutoff:

            # write to log file
            log_file.write(f"The species {name} has a busco completeness of {completeness} and so it is NOT being written to {my_files.output_dir}/busco_filtered_gammaproteobacteria.txt\n")


    # remove the unsorted metrics file
    if os.path.isfile(f"{my_files.output_dir}/all_species_busco_metrics.tsv"):

        os.remove(f"{my_files.output_dir}/all_species_busco_metrics.tsv")

        log_file.write(f"Found file {my_files.output_dir}/all_species_busco_metrics.tsv   Removing...\n")

    else:

        log_file.write(f"Error removing {my_files.output_dir}/all_species_busco_metrics.tsv Not a real file.\n")
