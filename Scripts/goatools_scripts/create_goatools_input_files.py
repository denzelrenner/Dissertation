# This script produces the study genes (the list of genes which we are loking for enrichment of processes in), the population genes (background list), and an associations
# file (all genes in pan-genome with GO terms). These 3 files serve as input for goatools

# import modules
import argparse
import os
import logging
import pandas as pd
import json
import shutil
import subprocess
import sys


vars = argparse.ArgumentParser(description='Creating a database to store all the genes, their cog catgories, and the species they were found in')

vars.add_argument('--output_directory',dest='output_dir', type=str, required=True, help='give the absolute path to the output directory to host all parsed output')

vars.add_argument('--positively_correlated_genes',dest='positively_correlated_genes',required=True ,type=str, help='file with only positively correlated genes from scoary output')

vars.add_argument('--negatively_correlated_genes',dest='negatively_correlated_genes',required=True ,type=str, help='file with only negatively correlated genes from scoary output')

vars.add_argument('--Rtab_file',dest='Rtab_file',required=True ,type=str, help='Rtab file from scoary')

vars.add_argument('--exclude_groups',dest='exclude_groups', default=False, action = 'store_true', help='set flag to ignore genes with group_10101 when writing to association file')

vars.add_argument('--annotation_file',dest='annotation_file', type=str, required=True, help='give the absolute path to the annotated pangenome file produced by emapper')


args = vars.parse_args()

# check if there is a directory for log files in the home directory and if not create one
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# configure the logger we will be using
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/goatools_prelim.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
logger = logging.getLogger(__name__)

# check if their output directories already exist or not
if not os.path.isdir(f"{args.output_dir}"):
    os.makedirs(f"{args.output_dir}")
    logger.info(f"Creating output directory {args.output_dir}")


# open the file with genes and their GOs, open the master file for positive genes, open an outfile to postively correlated genes that have no GO
with \
    open(f"{args.output_dir}/positively_correlated_genes_with_GOterm.txt",'w') as out_file_pos_genes, \
    open(f"{args.output_dir}/negatively_correlated_genes_with_GOterm.txt",'w') as out_file_neg_genes, \
    open(f"{args.output_dir}/pangenome_population_only.txt",'w') as outfile_pangenome_pops:
    # open(f"{args.output_dir}/pangenome_association_file.tsv",'w') as outfile_pangenome_association:

    # open file with genes and their GOs
    #gene_ontology_df = pd.read_csv(filepath_or_buffer=f"{args.gene_ontology}",sep='\t')

    # store only the positively correlated genes into a df
    positively_correlated_genes_df = pd.read_csv(filepath_or_buffer=f"{args.positively_correlated_genes}",sep='\t')

    # store all of those genes into a list, do the same for the whole dataset
    pos_correlated_genes = positively_correlated_genes_df['Gene'].unique().tolist()

    # store only the negatively correlated genes into a df
    negatively_correlated_genes_df = pd.read_csv(filepath_or_buffer=f"{args.negatively_correlated_genes}",sep='\t')

    # store all of those genes into a list, do the same for the whole dataset
    neg_correlated_genes = negatively_correlated_genes_df['Gene'].unique().tolist()

    # loop through all the pos correlated genes and check if they
    for gene in pos_correlated_genes:

        # # write name to out file if it has a GO
        # elif not temp_df.empty:
        out_file_pos_genes.write(f"{gene}\n")

    # loop through all the neg correlated genes and check if they
    for gene in neg_correlated_genes:

        # # write name to out file if it has a GO
        # elif not temp_df.empty:
        out_file_neg_genes.write(f"{gene}\n")


    #######################################################
    # CREATE POPULATION FILE WITH ALL GENES IN PAN-GENOME #
    #######################################################

    # take in Rtab, get all gene names, find if they have cogs
    pangenome_df = pd.read_csv(filepath_or_buffer=args.Rtab_file,sep='\t')

    # write all genes in the pangenome
    pangenome_genes = pangenome_df['Gene'].tolist()

    # sense check
    print(f"Number of genes in pangenome is {len(pangenome_genes)}")

    # loop through pangenome genes and write to the association file
    for gene in pangenome_genes:

        # write the gene to the population file regardless of if it has a GO term or not
        outfile_pangenome_pops.write(f"{gene}\n")

    ############################
    # CREATE ASSOCIATIONS FILE #
    ############################

    # code to do with writing to the population file has been removed from this file because it is not needed
    # outfile path
    associations_outfile_path = f"{args.output_dir}/pangenome_association_file.tsv"
    # pop_genes_outfile_path = f"{args.output_dir}/pangenomes_pop_genes_only.txt"
    # open(pop_genes_outfile_path,'w') as pop_genes_file

    # open a main output file
    with open(associations_outfile_path,'w') as assoc_file:

        # load in the emapper annotations file
        annotations_df = pd.read_csv(args.annotation_file,sep='\t',header=4,skipfooter=3,engine='python')
        print(annotations_df.head())
        # loop through each line of data
        for row in zip(annotations_df['#query'],annotations_df['GOs']):

            # set vars
            gene = row[0]
            # print(row[1])
            GO_terms = ';'.join(row[1].split(','))
            # print(row)
            if not args.exclude_groups:

                # add the gene to the total pop genes
                # pop_genes_file.write(f"{gene}\n")

                # now check if the GO term is empty or not. If it is empty dont write to the assocation file
                if GO_terms != '-':

                    # add gene and GO term to association file
                    assoc_file.write(f"{gene}\t{GO_terms}\n")

                elif GO_terms == '-':

                    print(f"{gene} doesnt have a GO term")

            # if i dont want any of the group genes
            elif args.exclude_groups:

                # if the gene has group in it then dont add it regardless of
                if 'group' in gene:
                    print(f"Excluding {gene} from population list")
                    continue

                elif not 'group' in gene:

                    # add the gene to the total pop genes
                    # pop_genes_file.write(f"{gene}\n")

                    # now check if the GO term is empty or not. If it is empty dont write to the assocation file
                    if GO_terms != '-':

                        # add gene and GO term to association file
                        assoc_file.write(f"{gene}\t{GO_terms}\n")

                    elif GO_terms == '-':

                        #
                        print(f"{gene} doesnt have a GO term")


