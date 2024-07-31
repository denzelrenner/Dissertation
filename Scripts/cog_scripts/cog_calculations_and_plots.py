# This script will take the emapper output from all the individual genomes, parse it to remove fields where theres no gene name.
# The script then filters the file created from the parameters mentioned above to combine the GO terms for genes. A gene may appear in different genomes so it will have multiple
# rows. So we are effectively combining the data from those rows

# import modules
import argparse
import os
import logging
import pandas as pd
import json
import shutil
import subprocess
import sys
from scipy.stats import fisher_exact

vars = argparse.ArgumentParser(description='Creating a database to store all the genes, their cog catgories, and the species they were found in')

vars.add_argument('-od','--output',dest='output_dir', type=str, required=True, help='give the absolute path to the output directory to host all parsed output')

vars.add_argument('--exclude_groups',dest='exclude_groups', default=False, action = 'store_true', help='set flag to ignore genes with group_10101 when writing to association file')

vars.add_argument('--annotation_file',dest='annotation_file', type=str, required=True, help='give the absolute path to the annotated pangenome file produced by emapper')

vars.add_argument('--Rtab_file',dest='Rtab_file',required=True ,type=str, help='Rtab file from scoary')

vars.add_argument('--positively_correlated_genes',dest='positively_correlated_genes',required=True ,type=str, help='file with only positively correlated genes from scoary output')

vars.add_argument('--negatively_correlated_genes',dest='negatively_correlated_genes',required=True ,type=str, help='file with only negatively correlated genes from scoary output')

vars.add_argument('--number_of_tests',dest='number_of_tests',required=True ,type=int, help='input number of tests being conducted so we can correct the pvalue threshold')

vars.add_argument('--pvalue',dest='pvalue_threshold',required=True ,type=float, help='significance threshold to be used')



args = vars.parse_args()


# check if there is a directory for log files in the home directory and if not create one
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# configure the logger we will be using
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/cog_category_calculation.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
logger = logging.getLogger(__name__)

# check if their output directories already exist or not
if not os.path.isdir(f"{args.output_dir}"):
    os.makedirs(f"{args.output_dir}")
    logger.info(f"Creating output directory {args.output_dir}")

# varibale for number of genomes processed
# number_of_genomes_processed = 0

# # set variable for json output data
# output_json_data = {}

# create json file that can be used
# cog_pairs = [('J', 'translation, including ribosome structure and biogenesis'), ('L','replication, recombination and repair'), ('K','transcription'), ('O','molecular chaperones and related functions'), ('M','cell wall structure and biogenesis and outer membrane'), ('N','secretion, motility and chemotaxis'),('T','signal transduction'),('P','inorganic ion transport and metabolism'),('C','energy production and conversion'),('G','carbohydrate metabolism and transport'),('E','amino acid metabolism and transport'),('F','nucleotide metabolism and transport'),('H','coenzyme metabolism'),('I','lipid metabolism'),('D','cell division and chromosome partitioning'),('R','general functional prediction only'),('S','no functional prediction'),('-','MISSING DATA')]

cog_pairs = [('A','RNA processing and modification'),
('B','Chromatin Structure and dynamics'),
('C','Energy production and conversion'),
('D','Cell cycle control and mitosis'),
('E','Amino Acid metabolis and transport'),
('F','Nucleotide metabolism and transport'),
('G','Carbohydrate metabolism and transport'),
('H','Coenzyme metabolism'),
('I','Lipid metabolism'),
('J','Tranlsation'),
('K','Transcription'),
('L','Replication and repair'),
('M','Cell wall/membrane/envelop biogenesis'),
('N','Cell motility'),
('O','Post-translational modification, protein turnover, chaperone functions'),
('P','Inorganic ion transport and metabolism'),
('Q','Secondary Structure'),
('T','Signal Transduction'),
('U','Intracellular trafficing and secretion'),
('Z','Cytoskeleton'),
('S','Function Unknown'),
('-','MISSING DATA'),
('V','Defence Mechanism'),
('W','Extracellular Structures')]

eggnog_json = {category:description for category,description in cog_pairs}

# log
logger.info(f'COG Descriptions: {eggnog_json}')

# define outfile paths
cog_genes_outfile_path = f"{args.output_dir}/stats_output.txt"
pos_cog_plotting_outfile_path = f"{args.output_dir}/plotting_input_cog_frequency_positive.txt"
neg_cog_plotting_outfile_path = f"{args.output_dir}/plotting_input_cog_frequency_negative.txt"


# store only the positively correlated genes into a df
positively_correlated_genes_df = pd.read_csv(filepath_or_buffer=f"{args.positively_correlated_genes}",sep='\t')

# store all of those genes into a list, do the same for the whole dataset
pos_cor_genes = positively_correlated_genes_df['Gene'].unique().tolist()

# store only the negatively correlated genes into a df
negatively_correlated_genes_df = pd.read_csv(filepath_or_buffer=f"{args.negatively_correlated_genes}",sep='\t')

# store all of those genes into a list, do the same for the whole dataset
neg_cor_genes = negatively_correlated_genes_df['Gene'].unique().tolist()

# for cog in cog_pairs:
#     print(cog[0])

# open a main output file
with open(cog_genes_outfile_path,'w') as cog_genes_file, open(pos_cog_plotting_outfile_path,'w') as pos_cog_plotting_file, open(neg_cog_plotting_outfile_path,'w') as neg_cog_plotting_file:

    # take in Rtab, get all gene names, find if they have cogs
    pangenome_df = pd.read_csv(filepath_or_buffer=args.Rtab_file,sep='\t')

    # create var for all genes in the pangenome
    pangenome_genes = pangenome_df['Gene'].tolist()

    # load in the emapper annotations file
    annotations_df = pd.read_csv(args.annotation_file,sep='\t',header=4,skipfooter=3,engine='python')

    # var for genes annotated by emapper
    all_emapper_genes = annotations_df['#query']

    print(annotations_df.head())

    # get frequency of COG categories
    cog_category_frequency_tracker_pangenome = {cog[0]:0 for cog in cog_pairs}
    cog_category_frequency_tracker_pos_genes = {cog[0]:0 for cog in cog_pairs}
    cog_category_frequency_tracker_neg_genes = {cog[0]:0 for cog in cog_pairs}

    # create vars for genes with no COG category for the different datasets
    no_cog_pangenome = []
    no_cog_pos = []
    no_cog_neg = []

    # now do the same but for COG categories
    for row in zip(annotations_df['#query'],annotations_df['COG_category']):

        # set vars
        gene = row[0]
        #print(row[1])
        multiple_cogs = list(row[1])

        # print(row)
        if not args.exclude_groups:

            # go through each COG and add it to the tally
            for COG in multiple_cogs:

                # increase cog tracker if the gene is in the gene set
                if gene in pangenome_genes:

                    # if COG == '-':
                    #     no_cog_pangenome.append(gene)
                    #     continue

                    cog_category_frequency_tracker_pangenome[COG] += 1

                if gene in pos_cor_genes:

                    # if COG == '-':
                    #     no_cog_pos.append(gene)
                    #     continue

                    cog_category_frequency_tracker_pos_genes[COG] += 1

                if gene in neg_cor_genes:

                    # if COG == '-':
                    #     no_cog_neg.append(gene)

                    cog_category_frequency_tracker_neg_genes[COG] += 1


            # if i dont want any of the group genes
        elif args.exclude_groups:

            # if the gene has group in it then dont add it regardless of
            if 'group' in gene:
                print(f"Excluding {gene} from population list")
                continue

            elif not 'group' in gene:

                # go through each COG and add it to the tally
                for COG in multiple_cogs:

                    # increase cog tracker if the gene is in the gene set
                    if gene in pangenome_genes:
                        cog_category_frequency_tracker_pangenome[COG] += 1

                    if gene in pos_cor_genes:
                        cog_category_frequency_tracker_pos_genes[COG] += 1

                    if gene in neg_cor_genes:
                        cog_category_frequency_tracker_neg_genes[COG] += 1


    # now do statistics
    # the contigency table will be
    #               Genes with COG  Genes without COG
    # pos correlated    i               j
    # pangenome         k               l

    # get total genes in all different studies
    total_pos_genes = len(pos_cor_genes)
    total_neg_genes = len(neg_cor_genes)
    total_pangenome_genes = len(pangenome_genes)

    # now that weve added the genes that have no COGs for a particular dataset, we can add the ones we counted before
    cog_category_frequency_tracker_pangenome['-'] = total_pangenome_genes - sum(cog_category_frequency_tracker_pangenome.values())
    cog_category_frequency_tracker_pos_genes['-'] = total_pos_genes - sum(cog_category_frequency_tracker_pos_genes.values())
    cog_category_frequency_tracker_neg_genes['-'] = total_neg_genes - sum(cog_category_frequency_tracker_neg_genes.values())

    logger.info(f"Total pos genes : {total_pos_genes} Total neg genes : {total_neg_genes} Total pangenome genes : {total_pangenome_genes}")

    # to get the number of missing cogs with -, count the number of COGs for everythign else and subtract from the total genes

    # write to main input file for plotting
       # tsv format is gene set type , then all the cogs, as header
    header_line = 'gene_set_type\tCOG_Category\tFrequency\tPercentage\n'
    pos_cog_plotting_file.write(header_line)
    neg_cog_plotting_file.write(header_line)

    # go through each cog
    for cog in cog_category_frequency_tracker_pangenome.keys():

        # create a line to write
        line_to_write = f"Pan-genome\t{cog}\t{cog_category_frequency_tracker_pangenome[cog]}\t{(cog_category_frequency_tracker_pangenome[cog]/total_pangenome_genes)*100}\n"

        # write to out file for both neg and pos
        pos_cog_plotting_file.write(line_to_write)
        neg_cog_plotting_file.write(line_to_write)

    for cog in cog_category_frequency_tracker_pos_genes.keys():

        # create a line to write
        line_to_write = f"Positive\t{cog}\t{cog_category_frequency_tracker_pos_genes[cog]}\t{(cog_category_frequency_tracker_pos_genes[cog]/total_pos_genes)*100}\n"

        # write to out
        pos_cog_plotting_file.write(line_to_write)

    for cog in cog_category_frequency_tracker_neg_genes.keys():

        # create a line to write
        line_to_write = f"Negative\t{cog}\t{cog_category_frequency_tracker_neg_genes[cog]}\t{(cog_category_frequency_tracker_neg_genes[cog]/total_neg_genes)*100}\n"

        # write to out
        neg_cog_plotting_file.write(line_to_write)

    # define what is considered significant
    significance_threshold = args.pvalue_threshold/args.number_of_tests

    logger.info(f"Original signficance threshold was {args.pvalue_threshold}. Corrected significance threshold is {significance_threshold}")

    # pangenome vs positively correlated genes
    # loop through each cog category
    for COG_Definition in cog_pairs:

        # store just cog in var
        cog = COG_Definition[0]

        # FOR POSITIVE GENES
        # do it for everything thats not -
        if cog != '-':
            pos_genes_with_cog = cog_category_frequency_tracker_pos_genes[cog]
            pos_genes_without_cog = total_pos_genes - cog_category_frequency_tracker_pos_genes[cog]
            pangenome_genes_with_cog = cog_category_frequency_tracker_pangenome[cog]
            pangenome_genes_without_cog = total_pangenome_genes - cog_category_frequency_tracker_pangenome[cog]

            print(cog,pos_genes_with_cog,pos_genes_without_cog,pangenome_genes_with_cog,pangenome_genes_without_cog)

            # create contingency table
            contingency_table = [[pos_genes_with_cog,pos_genes_without_cog],[pangenome_genes_with_cog,pangenome_genes_without_cog]]

            # do fishers exact test
            odds_ratio, pvalue = fisher_exact(contingency_table,alternative='two-sided')

            # output line to write, add the pangenome frequency as well
            if pvalue >= significance_threshold:
                line_to_write = f"Positive\t{cog}\t{pvalue}\t{odds_ratio}\tno_significance\t{(cog_category_frequency_tracker_pos_genes[cog]/total_pos_genes)*100}\n"

            elif pvalue < significance_threshold:
                line_to_write = f"Positive\t{cog}\t{pvalue}\t{odds_ratio}\tsignificant\t{(cog_category_frequency_tracker_pos_genes[cog]/total_pos_genes)*100}\n"

            cog_genes_file.write(line_to_write)

        # FOR NEGATIVE GENES
        # do it for everything thats not -
        if cog != '-':
            neg_genes_with_cog = cog_category_frequency_tracker_neg_genes[cog]
            neg_genes_without_cog = total_neg_genes - cog_category_frequency_tracker_neg_genes[cog]
            pangenome_genes_with_cog = cog_category_frequency_tracker_pangenome[cog]
            pangenome_genes_without_cog = total_pangenome_genes - cog_category_frequency_tracker_pangenome[cog]

            print(cog,neg_genes_with_cog,neg_genes_without_cog,pangenome_genes_with_cog,pangenome_genes_without_cog)

            # create contingency table
            contingency_table = [[neg_genes_with_cog,neg_genes_without_cog],[pangenome_genes_with_cog,pangenome_genes_without_cog]]

            # do fishers exact test
            odds_ratio, pvalue = fisher_exact(contingency_table,alternative='two-sided')

            # output line to write
            if pvalue >= significance_threshold:
                line_to_write = f"Negative\t{cog}\t{pvalue}\t{odds_ratio}\tno_significance\t{(cog_category_frequency_tracker_neg_genes[cog]/total_neg_genes)*100}\n"

            elif pvalue < significance_threshold:
                line_to_write = f"Negative\t{cog}\t{pvalue}\t{odds_ratio}\tsignificant\t{(cog_category_frequency_tracker_neg_genes[cog]/total_neg_genes)*100}\n"

            cog_genes_file.write(line_to_write)

    pan_genes = ','.join([str(cog) for cog in cog_category_frequency_tracker_pangenome.values()])
    pos_genes = ','.join([str(cog) for cog in cog_category_frequency_tracker_pos_genes.values()])
    neg_genes = ','.join([str(cog) for cog in cog_category_frequency_tracker_neg_genes.values()])

    print('Pangenome,','271220,',f"[{pan_genes}]")
    print('Positively Correlated,','27686,',f"[{pos_genes}]")
    print('Negatively Correlated,','16751,',f"[{neg_genes}]")



# create code to be written to add significance stars to either the positive vs pangenome plot or the negative vs pangenome plot
# open positive file
with open(cog_genes_outfile_path,'r') as stats_file,open(f"{args.output_dir}/COG_plotting.R",'w') as Rplotter_file:

    # check if it is significant,
    # create a subset
    # add a label
    # Set y value
    # set size
    # make vjust 0

    # create variables for the commands to write for plotting the positive and negative
    positive_plot_command_chunk = ''''''
    negative_plot_command_chunk = ''''''

    # store data in variable
    stats_data = stats_file.readlines()

    # loop through each line of data
    for line in stats_data:

        # set vars
        gene_set,cog,_,empty,significance,cog_freq = line.rstrip('\n').split('\t')

        # check if significant
        if significance == 'significant' and gene_set == 'Positive':

            # write command. create subset of main data file to get that specific COG for the positive
            command = f"geom_text(data=subset(df_pos,gene_set_type=='Pan-genome' & COG_Category=='{cog}'),vjust=0,size=17,aes(y={cog_freq}+2,label='*'))+\n"

            # add to the main command chunk
            positive_plot_command_chunk += command

        elif significance == 'significant' and gene_set == 'Negative':

            # write command. create subset of main data file to get that specific COG for the positive
            command = f"geom_text(data=subset(df_neg,gene_set_type=='Pan-genome' & COG_Category=='{cog}'),vjust=0,size=17,aes(y={cog_freq}+2,label='*'))+\n"

            # add to the main command chunk
            negative_plot_command_chunk += command

# The style of this plot was inspired by the Vibrio plot

# open file for R plotting


    ############################
    # START OF R PLOTTING FILE #
    ############################


    # store installation code in its own chunk because it is hard to write
    installation_chunk = '''if (!require("ggplot2", quietly = TRUE)) {
  print("The package ggplot2 is not installed. Beginning installation...")
  install.packages("ggplot2",dependencies=TRUE,repos='https://cloud.r-project.org')
} else {
  print("The package ggplot2 is installed.")
}

if (!require("gridExtra", quietly = TRUE)) {
  print("The package gridExtra is not installed. Beginning installation...")
  install.packages("gridExtra",dependencies=TRUE,repos='https://cloud.r-project.org')
} else {
  print("The package gridExtra is installed.")
}

if (!require("ggsci", quietly = TRUE)) {
  print("The package ggsci is not installed. Beginning installation...")
  install.packages("ggsci",repos='https://cloud.r-project.org')
} else {
  print("The package ggsci is installed.")
}

if (!require("colorspace", quietly = TRUE)) {
  print("The package colorspace is not installed. Beginning installation...")
  install.packages("colorspace", repos = 'https://cloud.r-project.org')
} else {
  print("The package colorspace is installed.")
}

if (!require("remotes", quietly = TRUE)) {
print("The package remotes is not installed. Beginning installation...")
install.packages("remotes", repos = 'https://cloud.r-project.org')
} else {
print("The package remotes is installed.")
}

remotes::install_github("wilkelab/cowplot")
install.packages("colorspace", repos = "http://R-Forge.R-project.org")
remotes::install_github("clauswilke/colorblindr")\n'''

    # write installation code to file
    Rplotter_file.write(installation_chunk)

    # library packages we will be using
    Rplotter_file.write('library(ggplot2)\n')
    Rplotter_file.write('library(gridExtra)\n')
    Rplotter_file.write('library(colorspace)\n')
    Rplotter_file.write('library(colorblindr)\n')
#    Rplotter_file.write('library(ggpubr)\n')

    # COG Categories (%)
    # set a working directory
    Rplotter_file.write(f"setwd('{args.output_dir}')\n")

    # load in the different files
    Rplotter_file.write(f"df_pos <- read.table(file='{pos_cog_plotting_outfile_path}',header=TRUE,sep='\\t')\n")
    Rplotter_file.write(f"df_neg <- read.table(file='{neg_cog_plotting_outfile_path}',header=TRUE,sep='\\t')\n")


    # start graphics driver gene_set_type\tCOG_Category\tFrequency\tPercentage\n
    # load graphics driver
#    Rplotter_file.write(f"pdf('{args.output_dir}/positive_top_GOterms.pdf', width = 12, height = 12)\n")
    Rplotter_file.write(f"pdf('{args.output_dir}/COG_category_plot_positive.pdf', width = 25, height = 25)\n")

    # pos genes plot
    Rplotter_file.write("ggplot(df_pos, aes(x = COG_Category, y = Percentage,fill=gene_set_type)) +\n")
    Rplotter_file.write("geom_bar(stat = 'identity',position = 'dodge',color='black') +\n")
    Rplotter_file.write("labs(x = 'COG Category',y = 'Percentage of COG Category (%)',fill = 'Dataset') +\n")
    Rplotter_file.write(f"{positive_plot_command_chunk}")
    # remove grid from background of plot
     # remove grid from background of plot
    Rplotter_file.write('theme(panel.background = element_rect(color="#FFFFFF", fill="white")) +\n')
    Rplotter_file.write('theme(panel.grid.minor = element_blank()) + \n')
    Rplotter_file.write('theme(panel.grid.major = element_blank()) + \n')
    Rplotter_file.write('''theme(axis.title.y = element_text(),axis.text = element_text(size = 40),
        axis.line.y = element_line(),
        axis.title=element_text(size = 50),legend.title = element_text(size=40),legend.key.size = unit(2, 'cm'),
        axis.ticks.y = element_line(),legend.text=element_text(size=40),
        axis.ticks.x = element_line(),
        axis.line.x = element_line())+\n''')
    Rplotter_file.write("scale_fill_manual(labels = c('Pan-genome', 'Positively Correlated Genes'),values=c('darkorange','purple')) +\n")
    Rplotter_file.write("scale_y_continuous(expand = c(0, 0), limits = c(0, 75))+\n")
    # flip axis
    #Rplotter_file.write("coord_flip()+\n")

    # centre plot tile
    Rplotter_file.write("theme(plot.title=element_text(hjust=0.5))\n")

    Rplotter_file.write("dev.off()\n")

    # neg genes plot
    Rplotter_file.write(f"pdf('{args.output_dir}/COG_category_plot_negative.pdf', width = 25, height = 25)\n")
    Rplotter_file.write("ggplot(df_neg, aes(x = COG_Category, y = Percentage,fill=gene_set_type)) +\n")
    Rplotter_file.write("geom_bar(stat = 'identity',position = 'dodge',color='black') +\n")
    Rplotter_file.write("labs(x = 'COG Category',y = 'Percentage of COG Category (%)',fill = 'Dataset') +\n")
    Rplotter_file.write(f"{negative_plot_command_chunk}")

    # remove grid from background of plot
     # remove grid from background of plot
    Rplotter_file.write('theme(panel.background = element_rect(color="#FFFFFF", fill="white")) +\n')
    Rplotter_file.write('theme(panel.grid.minor = element_blank()) + \n')
    Rplotter_file.write('theme(panel.grid.major = element_blank()) + \n')
    Rplotter_file.write('''theme(axis.title.y = element_text(),axis.text = element_text(size = 40),
        axis.line.y = element_line(),legend.title = element_text(size=40),legend.key.size = unit(2, 'cm'),
        axis.title=element_text(size = 50),legend.text=element_text(size=40),
        axis.ticks.y = element_line(),
        axis.ticks.x = element_line(),
        axis.line.x = element_line())+\n''')
    Rplotter_file.write("scale_fill_manual(labels = c('Negatively Correlated Genes','Pan-genome'),values=c('darkblue','darkorange')) +\n")
    Rplotter_file.write("scale_y_continuous(expand = c(0, 0), limits = c(0, 75))+\n")
    #Rplotter_file.write("scale_fill_OkabeIto(palette = 'contrast') +\n")

    # flip axis
    #Rplotter_file.write("coord_flip()+\n")

    # centre plot tile
    Rplotter_file.write("theme(plot.title=element_text(hjust=0.5))\n")

    Rplotter_file.write("dev.off()\n")
