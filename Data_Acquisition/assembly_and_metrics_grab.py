# This script is meant to test the pipeline idea i have where it takes in the different species and grabs their assemblies. 
# You have to be in the right directory and have the datasets conda environment downloaded. Move into the right output directory before running this 
# all directories are made in the shell script prior to running anything

# Example run command:
#                      python3 ../test_assembly_grab.py -i ../test_file.txt -d /Users/p3/Desktop/Final_Project/assembly_testing
# import modules
import argparse
import os

my_vars = argparse.ArgumentParser(description='Taking an input and output file')

my_vars.add_argument('-i','--input',dest='array_of_assembly_species', type=str, required=True, help='give the absolute path to the file containing the species you are interested in and their accessions')

# my_vars.add_argument('-o','--out',dest='output_file', type=str, metavar='output file', required=True,
#                         help='give the absolute path to the output file name with .out at the end')

my_vars.add_argument('-d','--odir',dest='output_dir', type=str, required=True, help='give an output directory to put all the assemblies information into')

my_files = my_vars.parse_args()

# put the value for ~ in a variable
# user = {os.path.expanduser('~')}

# make a new directory to put the individual metrics for each species into 
# os.system(f"mkdir {my_files.output_dir}/metric_files")

# make a new directory to put all the assembly data into 
# os.system(f"mkdir {my_files.output_dir}/assemblies")

# open the data file with assembly accessuin code and the species
with open(my_files.array_of_assembly_species,'r') as input_data, open(f'{my_files.output_dir}/output_data/all_species_metrics.txt','w+') as metrics_file :

    # write a header to the metrics file
    metrics_file.write('Species_Name\tGC_Content\tGC_Percentage\tAssembly_Length\n')
    # go through each line of data. start after the header
    for line in input_data.readlines()[1:]:

        # first split the line
        line_data = line.split()

        # get the assembly accession
        assembly_accession = line_data[0]

        # get the full species name. Join everything except the accessions
        species = ' '.join(line_data[1:])

        # naming with underscores for file naes
        file_name_species = '_'.join(line_data[1:]).replace("[","").replace("]","").replace("(","").replace(")","").replace("=","").replace("/","").replace(".","_")

        # print useful info
        print(f"This species is {species} and for files it is {file_name_species}")
        print(f"The accession code for the assembly is {assembly_accession}")

        # download the assembly, the protein fasta, the genbank format file, annoation file, the sequence report file
        os.system(f"datasets download genome accession {assembly_accession} --include genome,protein,gff3,gbff,seq-report --filename {file_name_species}_dataset.zip")

        # download statistics about the assembly like N50, L50, etc
        os.system(f"datasets summary genome accession {assembly_accession} --assembly-level complete --assembly-source 'RefSeq' --as-json-lines | dataformat tsv genome --fields accession,organism-name,assmstats-total-sequence-len,assmstats-contig-l50,assmstats-contig-n50,assmstats-gc-count,assmstats-gc-percent,assmstats-genome-coverage,assminfo-assembly-method,checkm-completeness > {my_files.output_dir}/metric_files/{file_name_species}.txt")

        # get the species name, GC content,GC precentage,and assembly length from the statistics file that was just downloaded (code above) and add it to the metrics file
        temp_file = open(f"{my_files.output_dir}/metric_files/{file_name_species}.txt",'r')

        temp_file_data = temp_file.readlines()[1].split('\t')

        GC_content,GC_percentage,assembly_length = temp_file_data[5],temp_file_data[6],temp_file_data[2]

        # print(temp_file_data.split('\t'))
        # write to metrics file
        metrics_file.write(f"{file_name_species}\t{GC_content}\t{GC_percentage}\t{assembly_length}\n")

        temp_file.close()

