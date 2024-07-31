# This script will take the GOATOOLS output tsv and separate it into Molecular Function, Biological Process, and Cellular Component
# The top 20 or 10 terms will then be taken from each individual Ontology and plotted as a rotated Barplot in R

# import modules
import argparse
import os
import logging
import pandas as pd
import numpy as np
import json
import shutil
import subprocess
import sys



vars = argparse.ArgumentParser(description='Creating a database to store all the genes, their cog catgories, and the species they were found in')

vars.add_argument('--output_directory',dest='output_dir', type=str, required=True, help='give the absolute path to the output directory to host all parsed output')

vars.add_argument('--negative_goatools',dest='negative_goatools', type=str, required=True, help='give the absolute path to the goatools output file with the negatively correlated genes as input')

vars.add_argument('--positive_goatools',dest='positive_goatools', type=str, required=True, help='give the absolute path to the goatools output file with the positively correlated genes as input')

vars.add_argument('--number_of_terms',dest='number_of_terms', type=int, required=True, help='number of significant terms to include in each plot')

args = vars.parse_args()

# check if there is a directory for log files in the home directory and if not create one
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# configure the logger we will be using
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/plotting_top_significant_GOs.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
logger = logging.getLogger(__name__)

# check if their output directories already exist or not
if not os.path.isdir(f"{args.output_dir}"):
    os.makedirs(f"{args.output_dir}")
    logger.info(f"Creating output directory {args.output_dir}")

# create file paths for the files with top15 terms for both pos and neg
combined_filepath_BP = f"{args.output_dir}/combined_terms_BP.tsv"
combined_filepath_MF = f"{args.output_dir}/combined_terms_MF.tsv"
combined_filepath_CC = f"{args.output_dir}/combined_terms_CC.tsv"

with open(combined_filepath_BP,'w') as combined_BP_file, open(combined_filepath_MF,'w') as combined_MF_file, open(combined_filepath_CC,'w') as combined_CC_file:

    # write headers for combined out file
    combined_BP_file.write(f"name\tp_fdr_bh\tcorrelation_type\n")
    combined_MF_file.write(f"name\tp_fdr_bh\tcorrelation_type\n")
    combined_CC_file.write(f"name\tp_fdr_bh\tcorrelation_type\n")

    # define output files for the pos and neg separately
    outfile_pos_BP = f"{args.output_dir}/positive_enriched_BP.tsv"
    outfile_pos_MF = f"{args.output_dir}/positive_enriched_MF.tsv"
    outfile_pos_CC = f"{args.output_dir}/positive_enriched_CC.tsv"

    outfile_neg_BP = f"{args.output_dir}/negative_enriched_BP.tsv"
    outfile_neg_MF = f"{args.output_dir}/negative_enriched_MF.tsv"
    outfile_neg_CC = f"{args.output_dir}/negative_enriched_CC.tsv"

    # load in the positive and negative dfs for GOATOOLS output
    df_pos = pd.read_csv(filepath_or_buffer=args.positive_goatools,sep='\t')
    df_neg = pd.read_csv(filepath_or_buffer=args.negative_goatools,sep='\t')

    print(df_pos[(df_pos['enrichment']=='p')])

    # separate the positive dfs into one for BP, MF, and CC, and arrange each df from smallest to largets pvalue. we dont include p values of 0 because that is the big terms like biological process which all terms have
    df_pos_BP = df_pos[(df_pos['NS']=='BP') & (df_pos['p_fdr_bh']!=0.0) & (df_pos['enrichment']=='e')].sort_values(by=['p_fdr_bh'],ascending=True).head(args.number_of_terms).rename(columns={'# GO':'GO'})
    df_pos_MF = df_pos[(df_pos['NS']=='MF') & (df_pos['p_fdr_bh']!=0.0) & (df_pos['enrichment']=='e')].sort_values(by=['p_fdr_bh'],ascending=True).head(args.number_of_terms).rename(columns={'# GO':'GO'})
    df_pos_CC = df_pos[(df_pos['NS']=='CC') & (df_pos['p_fdr_bh']!=0.0) & (df_pos['enrichment']=='e')].sort_values(by=['p_fdr_bh'],ascending=True).head(args.number_of_terms).rename(columns={'# GO':'GO'})

    print(df_pos_BP)
    print(df_pos_MF)

    # separate the negative dfs into one for BP, MF, and CC, and arrange each df from smallest to largets pvalue
    df_neg_BP = df_neg[(df_neg['NS']=='BP') & (df_neg['p_fdr_bh']!=0.0) & (df_neg['enrichment']=='e')].sort_values(by=['p_fdr_bh'],ascending=True).head(args.number_of_terms).rename(columns={'# GO':'GO'})
    df_neg_MF = df_neg[(df_neg['NS']=='MF') & (df_neg['p_fdr_bh']!=0.0) & (df_neg['enrichment']=='e')].sort_values(by=['p_fdr_bh'],ascending=True).head(args.number_of_terms).rename(columns={'# GO':'GO'})
    df_neg_CC = df_neg[(df_neg['NS']=='CC') & (df_neg['p_fdr_bh']!=0.0) & (df_neg['enrichment']=='e')].sort_values(by=['p_fdr_bh'],ascending=True).head(args.number_of_terms).rename(columns={'# GO':'GO'})

    print(df_neg_BP)
    print(df_neg_MF)


    #####################################
    # FOR BIOLOGICAL PROCESSES ONTOLOGY #
    #####################################
    # when combining dfs for plotting, by default all positives will be to the right, and negatives to the left
    # combine the dfs for bp

    # first get only the GO terms from each list
    pos_BP_GOterms = df_neg_BP['name'].tolist()
    neg_BP_GOterms = df_neg_BP['name'].tolist()

    # put the terms for the pos and neg to a list
    combined_GOterms_BP = df_neg_BP['name'].tolist() + df_pos_BP['name'].tolist()

    # Positive Correlated Genes #

    # for term in the list of top15 go terms for positively correlated genes
    for row in zip(df_pos_BP['name'],df_pos_BP['p_fdr_bh']):

        GOterm = row[0]

        # log transform the pvalue
        pvalue = np.log10(float(row[1]))

        # write to out file
        combined_BP_file.write(f"{GOterm}\t{pvalue}\tpositive\n")

    # now go through each negative term, and if it is not in the top15 positive terms, add a row for the positives with that term with a signficance of 0
    for term in neg_BP_GOterms:

        if term not in pos_BP_GOterms:

            # write to outfile
            combined_BP_file.write(f"{term}\t0\tpositive\n")

    # Negative Correlated Genes #

    # for term in the list of top15 go terms for negatively correlated genes
    for row in zip(df_neg_BP['name'],df_neg_BP['p_fdr_bh']):

        GOterm = row[0]

        # log transform the pvalue
        pvalue = np.log10(float(row[1])) * -1

        # write to out file
        combined_BP_file.write(f"{GOterm}\t{pvalue}\tnegative\n")

    # now go through each positive term, and if it is not in the top15 positive terms, add a row for the positives with that term with a signficance of 0
    for term in pos_BP_GOterms:

        if term not in neg_BP_GOterms:

            # write to outfile
            combined_BP_file.write(f"{term}\t0\tnegative\n")

    ##################################
    # FOR MOLECULAR FUNCTION ONTOLOGY#
    ##################################

    # first get only the GO terms from each list
    pos_MF_GOterms = df_neg_MF['name'].tolist()
    neg_MF_GOterms = df_neg_MF['name'].tolist()

    # put the terms for the pos and neg to a list
    combined_GOterms_MF = df_neg_MF['name'].tolist() + df_pos_MF['name'].tolist()

    # Positive Correlated Genes #

    # for term in the list of top15 go terms for positively correlated genes
    for row in zip(df_pos_MF['name'],df_pos_MF['p_fdr_bh']):

        GOterm = row[0]

        # log transform the pvalue
        pvalue = np.log10(float(row[1]))

        # write to out file
        combined_MF_file.write(f"{GOterm}\t{pvalue}\tpositive\n")

    # now go through each negative term, and if it is not in the top15 positive terms, add a row for the positives with that term with a signficance of 0
    for term in neg_MF_GOterms:

        if term not in pos_MF_GOterms:

            # write to outfile
            combined_MF_file.write(f"{term}\t0\tpositive\n")

    # Negative Correlated Genes #

    # for term in the list of top15 go terms for negatively correlated genes
    for row in zip(df_neg_MF['name'],df_neg_MF['p_fdr_bh']):

        GOterm = row[0]

        # log transform the pvalue
        pvalue = np.log10(float(row[1])) * -1

        # write to out file
        combined_MF_file.write(f"{GOterm}\t{pvalue}\tnegative\n")

    # now go through each positive term, and if it is not in the top15 positive terms, add a row for the positives with that term with a signficance of 0
    for term in pos_MF_GOterms:

        if term not in neg_MF_GOterms:

            # write to outfile
            combined_MF_file.write(f"{term}\t0\tnegative\n")

    ###################################
    # FOR CELLULAR COMPONENT ONTOLOGY #
    ###################################

    # first get only the GO terms from each list
    pos_CC_GOterms = df_neg_CC['name'].tolist()
    neg_CC_GOterms = df_neg_CC['name'].tolist()

    # put the terms for the pos and neg to a list
    combined_GOterms_MF = df_neg_CC['name'].tolist() + df_pos_CC['name'].tolist()

    # Positive Correlated Genes #

    # for term in the list of top15 go terms for positively correlated genes
    for row in zip(df_pos_CC['name'],df_pos_CC['p_fdr_bh']):

        GOterm = row[0]

        # log transform the pvalue
        pvalue = np.log10(float(row[1]))

        # write to out file
        combined_CC_file.write(f"{GOterm}\t{pvalue}\tpositive\n")

    # now go through each negative term, and if it is not in the top15 positive terms, add a row for the positives with that term with a signficance of 0
    for term in neg_CC_GOterms:

        if term not in pos_CC_GOterms:

            # write to outfile
            combined_CC_file.write(f"{term}\t0\tpositive\n")

    # Negative Correlated Genes #

    # for term in the list of top15 go terms for negatively correlated genes
    for row in zip(df_neg_CC['name'],df_neg_CC['p_fdr_bh']):

        GOterm = row[0]

        # log transform the pvalue, multiply by -1 so it becomes positive
        pvalue = np.log10(float(row[1])) * -1

        # write to out file
        combined_CC_file.write(f"{GOterm}\t{pvalue}\tnegative\n")

    # now go through each positive term, and if it is not in the top15 positive terms, add a row for the positives with that term with a signficance of 0
    for term in pos_CC_GOterms:

        if term not in neg_CC_GOterms:

            # write to outfile
            combined_CC_file.write(f"{term}\t0\tnegative\n")

    # now write the output to a tsv file which will be used
    df_pos_BP.to_csv(path_or_buf=outfile_pos_BP,sep='\t',index=False)
    df_pos_MF.to_csv(path_or_buf=outfile_pos_MF,sep='\t',index=False)
    df_pos_CC.to_csv(path_or_buf=outfile_pos_CC,sep='\t',index=False)
    df_neg_BP.to_csv(path_or_buf=outfile_neg_BP,sep='\t',index=False)
    df_neg_MF.to_csv(path_or_buf=outfile_neg_MF,sep='\t',index=False)
    df_neg_CC.to_csv(path_or_buf=outfile_neg_CC,sep='\t',index=False)


# now reopen the files and combine the negative and positive for a given process


# open file for R plotting
with open(f"{args.output_dir}/top_GOs_plotting.R",'w') as Rplotter_file:

    ############################
    # START OF R PLOTTING FILE #
    ############################

    # remotes::install_github("wilkelab/cowplot")

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

    # set a working directory
    Rplotter_file.write(f"setwd('{args.output_dir}')\n")

    # load in the different files
    Rplotter_file.write(f"df_BP <- read.table(file='{combined_filepath_BP}',header=TRUE,sep='\\t')\n")
    Rplotter_file.write(f"df_MF <- read.table(file='{combined_filepath_MF}',header=TRUE,sep='\\t')\n")
    Rplotter_file.write(f"df_CC <- read.table(file='{combined_filepath_CC}',header=TRUE,sep='\\t')\n")

    # start graphics driver
    # load graphics driver
#    Rplotter_file.write(f"pdf('{args.output_dir}/positive_top_GOterms.pdf', width = 12, height = 12)\n")
    Rplotter_file.write(f"pdf('{args.output_dir}/BP_combined_top_{args.number_of_terms}_GOterms.pdf', width = 30, height = 30)\n")

    Rplotter_file.write("ggplot(df_BP, aes(x = name, y = p_fdr_bh,fill=correlation_type)) +\n")
    Rplotter_file.write("geom_bar(stat = 'identity',position = 'identity') +\n")
    Rplotter_file.write("labs(title = 'Biological Processes',x = 'GO Term',y = '-log10(p-value)',fill = 'Correlation Type') +\n")

    # centre plot tile
#    Rplotter_file.write("theme(plot.title=element_text(hjust=0.5))+\n")

    # remove the grid from backgorund of plot
    # Rplotter_file.write("theme_minimal() +\n")
    Rplotter_file.write("theme_bw(base_size = 18)+\n")
    # Rplotter_file.write('theme(panel.background = element_rect(color="#FFFFFF", fill="white"),\n')
    Rplotter_file.write('theme(axis.line = element_line(colour = "black"),')
    Rplotter_file.write("panel.border = element_blank(),\n")
    Rplotter_file.write("panel.grid.minor = element_blank(),\n")
    Rplotter_file.write("panel.grid.major = element_blank(),\n")
    Rplotter_file.write("axis.text = element_text(size=40), \n")
    Rplotter_file.write("plot.title = element_text(size=50), \n")
    Rplotter_file.write("legend.text = element_text(size=50), \n")
    Rplotter_file.write("legend.title = element_text(size=50),\n")
    Rplotter_file.write("axis.title = element_text(size=65)) +\n")
    # flip plot and remove -values from x axis
    Rplotter_file.write("coord_flip() +\n")
    Rplotter_file.write("scale_y_continuous(labels=abs) +\n")
    Rplotter_file.write("theme(plot.title=element_text(hjust=0.5),plot.title.position = 'plot')+\n")
    # colour based on colour blind friendly colours
    Rplotter_file.write("scale_fill_OkabeIto()\n")
    #Rplotter_file.write("scale_fill_manual(values=cb_palette()) +\n")
    Rplotter_file.write("dev.off()\n")

    # Rplotter_file.write(f"pdf('{args.output_dir}/combined_top_{args.number_of_terms}_GOterms.pdf', width = 25, height = 25)\n")
    Rplotter_file.write(f"pdf('{args.output_dir}/CC_combined_top_{args.number_of_terms}_GOterms.pdf', width = 30, height = 30)\n")

    Rplotter_file.write("ggplot(df_CC, aes(x = name, y = p_fdr_bh,fill=correlation_type)) +\n")
    Rplotter_file.write("geom_bar(stat = 'identity',position = 'identity') +\n")
    # Rplotter_file.write("scale_fill_gradient(low = 'lightblue', high = 'darkblue') +\n")
    Rplotter_file.write("labs(title = 'Cellular Component',x = 'GO Term',y = '-log10(p-value)',fill = 'Correlation Type') +\n")
    Rplotter_file.write("theme_bw(base_size = 18)+\n")
    #Rplotter_file.write("theme_minimal() +\n")
    # Rplotter_file.write('theme(panel.background = element_rect(color="#FFFFFF", fill="white"))+\n')
    # Rplotter_file.write("theme(panel.grid.minor = element_blank())+\n")
    # Rplotter_file.write("theme(panel.grid.major = element_blank())+\n")
    Rplotter_file.write('theme(axis.line = element_line(colour = "black"),')
    Rplotter_file.write("panel.border = element_blank(),\n")
    Rplotter_file.write("panel.grid.minor = element_blank(),\n")
    Rplotter_file.write("panel.grid.major = element_blank(),\n")
    Rplotter_file.write("axis.text = element_text(size=40), \n")
    Rplotter_file.write("legend.text = element_text(size=50), \n")
    Rplotter_file.write("plot.title = element_text(size=50), \n")
    Rplotter_file.write("legend.title = element_text(size=50),\n")
    Rplotter_file.write("axis.title = element_text(size=65)) +\n")
    #Rplotter_file.write("theme(axis.text.x = element_text(angle = 45, hjust = 1)) +\n")
    Rplotter_file.write("coord_flip() +\n")
    Rplotter_file.write("scale_y_continuous(labels=abs) +\n")
    Rplotter_file.write("theme(plot.title=element_text(hjust=0.5),plot.title.position = 'plot')+\n")
    Rplotter_file.write("scale_fill_OkabeIto()\n")
    #Rplotter_file.write("scale_fill_manual(values=cb_palette()) +\n")

    Rplotter_file.write("dev.off()\n")

    Rplotter_file.write(f"pdf('{args.output_dir}/MF_combined_top_{args.number_of_terms}_GOterms.pdf', width = 30, height = 30)\n")
    Rplotter_file.write("ggplot(df_MF, aes(x = name, y = p_fdr_bh,fill=correlation_type)) +\n")
    Rplotter_file.write("geom_bar(stat = 'identity',position = 'identity') +\n")
    # Rplotter_file.write("scale_fill_gradient(low = 'lightblue', high = 'darkblue') +\n")
    Rplotter_file.write("labs(title = 'Molecular Function',x = 'GO Term',y = '-log10(p-value)',fill = 'Correlation Type') +\n")
    #Rplotter_file.write("theme_minimal() +\n")
    Rplotter_file.write("theme_bw(base_size = 18) +\n")
    Rplotter_file.write('theme(axis.line = element_line(colour = "black"),')
    Rplotter_file.write("panel.border = element_blank(),\n")
    Rplotter_file.write("panel.grid.minor = element_blank(),\n")
    Rplotter_file.write("panel.grid.major = element_blank(),\n")
#    Rplotter_file.write("axis.text.x = element_text(size=20)) +\n")
    Rplotter_file.write("axis.text = element_text(size=40), \n")
    Rplotter_file.write("legend.text = element_text(size=50), \n")
    Rplotter_file.write("plot.title = element_text(size=50), \n")
    Rplotter_file.write("legend.title = element_text(size=50),\n")
    Rplotter_file.write("axis.title = element_text(size=65)) +\n")
    Rplotter_file.write("coord_flip() +\n")
    Rplotter_file.write("scale_y_continuous(labels=abs) +\n")
    Rplotter_file.write("theme(plot.title=element_text(hjust=0.5),plot.title.position = 'plot')+\n")
    Rplotter_file.write("scale_fill_OkabeIto()\n")
    #Rplotter_file.write("scale_fill_manual(values=cb_palette()) +\n")


    # Rplotter_file.write("grid.arrange(plot1,plot2,plot3,ncol = 2, nrow = 2)\n")

    Rplotter_file.write("dev.off()\n")
