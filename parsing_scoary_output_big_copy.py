# This script takes the csv files from scoary and reformats them so that it is easier to read and subset data from it

# # import modules
import argparse
import os
import json
import pandas as pd
import numpy as np
import logging
import csv
import sys
from io import StringIO

my_vars = argparse.ArgumentParser(description='Interactive viewer of scoary output')

my_vars.add_argument('-c','--csv',dest='csv_file', type=str,nargs = '+', default=None, required=False, help='give the absolute path to the file containing the data matrix')

my_vars.add_argument('-d','--cdir',dest='csv_directory', type=str, default=None, required=False, help='give the absolute path to the directory of scoary output genes you are interested in')

my_vars.add_argument('-iv','--interactive',dest='interactive_usage', default=False, required=False,action='store_true',help='True or False on if you want to use an interactive session')

my_vars.add_argument('--output_directory',dest='output_dir', type=str, required=True, help='absolute path to output directory must be given')

my_vars.add_argument('--exclude_group',dest='exclude_group_genes', default=False,action='store_true', required=False, help='set flag if group genes (i.e group123) should not be considered when writing to out file')

my_vars.add_argument('-j',dest='json_file',type=str, required=True, help='absolute path to the json file that convert the trait names back to their regular one')

# my_vars.add_argument('-pvalue_correction',dest='pvalue_correction',type=str, required=True,options=['Bonferonni','Benjamini'], help='p value corr')



my_files = my_vars.parse_args()


# increase line limit
maxInt = sys.maxsize

while True:

    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

# create log file directory if it doesnt exist
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# create output directories if they dont exist
if not os.path.isdir(f"{my_files.output_dir}/specific_filtered_files"):
    os.makedirs(f"{my_files.output_dir}/specific_filtered_files")

# configure the logger we will be using
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/parsing_scoary_output_big_copy.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
logger = logging.getLogger(__name__)

# create output directory if it doesnt exist
# if not os.path.isdir(my_files.output_dir"):
#     os.makedirs(my_files.output_dir)

# log line limit
logger.info(f"Csv val confirmation:{csv.field_size_limit()}")

# transform number of strains with trait positivity/negativity and gene absence/presence into a percentage
#def create_percentage(number,total_number_of_strains=400):

    # do divison
#    percentage = float(int(number)/total_number_of_strains) * 100

    # return the rounded percentage
#    return round(percentage,2)

# function to process the assiocations and pvalues for a given arg trait file
def parse_arg_results_csv(files:list):

    # define significance level for bonferonni correction
    significance_threshold = 0.05
    no_of_tests_carried_out = 70
    bonferonni_corrected_signficance_threshold = significance_threshold/no_of_tests_carried_out

    # open a main file to put data from all the different csvs into one place, create two extra files to have the negative correlated genes, and the positively correlated genes
    with \
        open(f"{my_files.output_dir}/specific_filtered_files/gene_families_associated_genes_and_annotations.tsv",'w') as master_output_file, \
        open(my_files.json_file,'r') as json_file, \
        open(f"{my_files.output_dir}/specific_filtered_files/bonferonni_negative_correlated_genes.tsv",'w') as negative_correlated_genes, \
        open(f"{my_files.output_dir}/specific_filtered_files/bonferonni_positive_correlated_genes.tsv",'w') as positive_correlated_genes, \
        open(f"{my_files.output_dir}/specific_filtered_files/bonferonni_no_correlation_genes.tsv",'w') as no_correlated_genes:

        # write column header to output files
        master_output_file.writelines(["Gene family\t","Gene\t","Annotation\t","Odds_ratio\t","Naive_p\t","Bonferroni_p\t","Benjamini_H_p\t","Percentage_pos_present_in\t","Percentage_neg_present_in\t","Percentage_pos_not_present_in\t","Percentage_neg_not_present_in\t","Correlation Type\n"])
        positive_correlated_genes.writelines(["Gene family\t","Gene\t","Annotation\t","Odds_ratio\t","Bonferroni_p\n"])
        negative_correlated_genes.writelines(["Gene family\t","Gene\t","Annotation\t","Odds_ratio\t","Bonferroni_p\n"])
        no_correlated_genes.writelines(["Gene family\t","Gene\t","Annotation\t","Odds_ratio\t","Bonferroni_p\n"])

        # store dictionary which has trait name as key and the actual gene family name as value
        traits_to_scientific_name = json.load(json_file)

        # create dictionary to store the gene family from the traits file as the key, and all its associated genes/statistics as the value
        # all_association_data = {}

        # loop through each file that has been given. When a single csv file is given, it is given as a list so this method still works
        for file in files:

            # check if it is actually a scoary output file
            if file.endswith('.results.csv'):

                # log
                logger.info(f"Opening {file}...")

                # get the name of the main gene family. we need to split ABSOLUTE the file path, remove the results.csv from the end, and then the data and random 4 digits
                query_gene_family = traits_to_scientific_name[file.split('/')[-1].replace('.results.csv','')[:-16]]
                out_file_name_for_family = file.split('/')[-1].replace('.results.csv','')[:-16]

                # create a new output file for the better formatted output file
                output_file_path = f'{my_files.output_dir}/{out_file_name_for_family}.tsv'

                # open the csv file and the nnew output file
                with open(file,'r') as ARG_file, open(output_file_path,'w') as out_file:

                    # write column header for output file
                    out_file.writelines(["Gene\t","Annotation\t","Odds_ratio\t","Naive_p\t","Bonferroni_p\t","Benjamini_H_p\t","Percentage_pos_present_in\t","Percentage_neg_present_in\t","Percentage_pos_not_present_in\t","Percentage_neg_not_present_in\t","Correlation Type\n"])

                    # get list of headers
                    # csv_column_headers = ["Gene","Non-unique Gene name","Annotation","Number_pos_present_in","Number_neg_present_in","Number_pos_not_present_in","Number_neg_not_present_in","Sensitivity","Specificity","Odds_ratio","Naive_p","Bonferroni_p","Benjamini_H_p","Max_Pairwise_comparisons","Max_supporting_pairs","Max_opposing_pairs","Best_pairwise_comp_p","Worst_pairwise_comp_p"]

                    # load in the csv file from scoary
                    ARG_associations = csv.reader(ARG_file,delimiter=',',skipinitialspace=True)

                    # get the header line in the csv. This also moves from the header line to the start of the actual data
                    header = next(ARG_associations)

                    # get the index for all the different columns
                    gene_col = header.index("Gene")
                    annotation_col = header.index("Annotation")
                    odds_ratio = header.index("Odds_ratio")
                    naive_pvalue = header.index("Naive_p")
                    benferonni_pvalue = header.index("Bonferroni_p")
                    benjamini_hoch_pvalue = header.index("Benjamini_H_p")

                    # finding the percentage of strains with these combinations. change number of strains
                    numberofstrains_traitpositive_and_genepresent = header.index("Number_pos_present_in")
                    numberofstrains_traitnegative_and_genepresent = header.index("Number_neg_present_in")
                    numberofstrains_traitpositive_and_geneabsent = header.index("Number_pos_not_present_in")
                    numberofstrains_traitnegative_and_geneabsent = header.index("Number_neg_not_present_in") # when gene was absent, trait was negative

                    # determine the association type. whether the gene presence or absence is correlated with the trait
                    correlation_type = None

                    # got through the csv
                    for line in ARG_associations:

                        if line[odds_ratio] == 'inf':
                            correlation_type = 'positive'
                         # if odds ratio is 1 the odds of gene family present is the same in the presence or absence of a gene
                        elif float(line[odds_ratio]) == 1:
                            correlation_type = 'none'

                        # absence of the gene is negatively associated with the trait. gene family present, gene absent
                        elif float(line[odds_ratio]) < 1:
                            correlation_type = 'negative'

                        # presence of the gene is positively associated with the trait. gene family present, gene present
                        elif float(line[odds_ratio]) > 1:
                            correlation_type = 'positive'
                        # plotting file variables. Unused varibale but allows me to see everything i want to write to output
                        # plotting_variables = [line[gene_col],line[annotation_col],line[odds_ratio],line[naive_pvalue],line[benferonni_pvalue],line[benjamini_hoch_pvalue]

                        # take the associated gene name and add it as the key for the dictionary within the dictionary for the main query gene family.
                        associated_gene = line[gene_col].replace('"','')

                        # handle excluding group argument
                        if my_files.exclude_group_genes:

                            # only get the data for the associated gene if the gene name is a real gene and not a group_
                            if not associated_gene.startswith('group'):

                                # I just realised when writing to the out files there is no reason to keep columns i wont be using and i coudl come back and change it at any time
                                # write lines writes a list to the output file. I am using list comprehension rather than just giving the list of values because i need to add tabs
                                # and a new line if it is the last eleement of the list
                                # out_file.writelines([f"{new_line_data[i]}\n" if i == len(new_line_data)-1 else f"{new_line_data[i]}\t" for i in range(len(new_line_data))])
                                out_file.writelines([f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[naive_pvalue]}\t",f"{line[benferonni_pvalue]}\t",f"{line[benjamini_hoch_pvalue]}\t",f"{line[numberofstrains_traitpositive_and_genepresent]}\t",f"{line[numberofstrains_traitnegative_and_genepresent]}\t",f"{line[numberofstrains_traitpositive_and_geneabsent]}\t",f"{line[numberofstrains_traitnegative_and_geneabsent]}\t",f"{correlation_type}\n"])

                                # write to master output file with just gene fmaily, associate gene and the first annotation
                                master_output_file.writelines([f"{query_gene_family}\t",f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[naive_pvalue]}\t",f"{line[benferonni_pvalue]}\t",f"{line[benjamini_hoch_pvalue]}\t",f"{line[numberofstrains_traitpositive_and_genepresent]}\t",f"{line[numberofstrains_traitnegative_and_genepresent]}\t",f"{line[numberofstrains_traitpositive_and_geneabsent]}\t",f"{line[numberofstrains_traitnegative_and_geneabsent]}\t",f"{correlation_type}\n"])

                                # handle if a significant p value is found for bonferonni
                                if float(line[benferonni_pvalue]) < bonferonni_corrected_signficance_threshold:

                                    # write to a different file depending on what the correlation type is
                                    if correlation_type == 'positive':
                                        positive_correlated_genes.writelines([f"{query_gene_family}\t",f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[benferonni_pvalue]}\n"])

                                    elif correlation_type == 'negative':
                                        negative_correlated_genes.writelines([f"{query_gene_family}\t",f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[benferonni_pvalue]}\n"])

                                    elif correlation_type == 'none':
                                        no_correlated_genes.writelines([f"{query_gene_family}\t",f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[benferonni_pvalue]}\n"])


                        # if the exclude group genes isnt flagged
                        elif not my_files.exclude_group_genes:

                            # write to R file that will be taken in by R
                            out_file.writelines([f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[naive_pvalue]}\t",f"{line[benferonni_pvalue]}\t",f"{line[benjamini_hoch_pvalue]}\t",f"{line[numberofstrains_traitpositive_and_genepresent]}\t",f"{line[numberofstrains_traitnegative_and_genepresent]}\t",f"{line[numberofstrains_traitpositive_and_geneabsent]}\t",f"{line[numberofstrains_traitnegative_and_geneabsent]}\t",f"{correlation_type}\n"])

                            # write to master output file with just gene fmaily, associate gene and, annotation
                            # write to master output file with just gene fmaily, associate gene and the first annotation
                            master_output_file.writelines([f"{query_gene_family}\t",f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[naive_pvalue]}\t",f"{line[benferonni_pvalue]}\t",f"{line[benjamini_hoch_pvalue]}\t",f"{line[numberofstrains_traitpositive_and_genepresent]}\t",f"{line[numberofstrains_traitnegative_and_genepresent]}\t",f"{line[numberofstrains_traitpositive_and_geneabsent]}\t",f"{line[numberofstrains_traitnegative_and_geneabsent]}\t",f"{correlation_type}\n"])

                            # handle if a significant p value is found for bonferonni
                            if float(line[benferonni_pvalue]) < bonferonni_corrected_signficance_threshold:

                                # write to a different file depending on what the correlation type is
                                if correlation_type == 'positive':
                                    positive_correlated_genes.writelines([f"{query_gene_family}\t",f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[benferonni_pvalue]}\n"])

                                elif correlation_type == 'negative':
                                    negative_correlated_genes.writelines([f"{query_gene_family}\t",f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[benferonni_pvalue]}\n"])

                                elif correlation_type == 'none':
                                    no_correlated_genes.writelines([f"{query_gene_family}\t",f"{line[gene_col]}\t",f"{line[annotation_col].split(',')[0]}\t",f"{line[odds_ratio]}\t",f"{line[benferonni_pvalue]}\n"])

###############################################################
# PARSE INFORMATION FROM THE DIFFERENT CSVs PRODUCED BY SCOARY#
###############################################################

# if we are only dealing with a single csv file
if my_files.csv_file:

    # log
    logger.info(f"{len(my_files.csv_file)} files provided. running operation on one csv file only")

    # run function
    parse_arg_results_csv(files = my_files.csv_file)

# if a direcotry with the scoary output for the different traits are given
if my_files.csv_directory:

    # log
    logger.info(f"full directory given. listing directory to get file names. Beginning parsing...")

    # store all files in a variable
    file_names = os.listdir(my_files.csv_directory)

    # add the full path for each file int he csv directory
    all_files = [f"{my_files.csv_directory}/{file}" for file in file_names]

    # if the files end with .results.csv then open it and process it like with if it was only one file
    # run function
    parse_arg_results_csv(files = all_files)


########################################################
# PUTTING SIGNIFICANT GENES IN THEIR OWN FILE AS A JSON#
########################################################

# file to write output to
with open(f"{my_files.output_dir}/specific_filtered_files/significant_genes.json",'w') as json_output:

    # create a json to store the gene as key and the list of gene families as the value. These are only for signficant associations using benferonni
    outfile_data = {}

    # open the master file with all the arg gene families and gene associations
    # we are going to find all the unqiue gene, and see which species they are in
    df = pd.read_csv(filepath_or_buffer=f"{my_files.output_dir}/specific_filtered_files/gene_families_associated_genes_and_annotations.tsv",header=0,sep='\t')

    # put unique genes in a variable
    unique_genes = df['Gene'].unique().tolist()
    #print(unique_genes)

    # benferonni cut off is the p value threshold divided by pangenome size
    significance_threshold = 0.05
    no_of_tests_carried_out = 70
    bonferonni_corrected_signficance_threshold = significance_threshold/no_of_tests_carried_out

    # go through each line and find rows of data that have the gene and a significant p value
    for gene in unique_genes:

        # extract all lines with the gene name, and a significant p value. This way we can find
        all_hits = df[(df['Gene'] == gene) & (df['Bonferroni_p'] < bonferonni_corrected_signficance_threshold)]

        # then for a given gene, i want to print out all the gene families it was significant in
        rows = all_hits.values.tolist()

        # get all the genes families the gene had significant associations for, and also get the direction of correlation. last value in row is correlation type so we can use -1
        gene_families_only = [row[0] for row in rows]
        gene_families_and_correlation = [f"{row[0]}:{row[-1]}" for row in rows]

        # only store if gene families were found,i.e there were significnat associations with a gene family
        if gene_families_only:

            # store gene name and list of gene family traits where this gene was significant to the json
            outfile_data[gene] = gene_families_and_correlation

    # dump dict to json
    json.dump(obj=outfile_data,indent=4,fp=json_output)
