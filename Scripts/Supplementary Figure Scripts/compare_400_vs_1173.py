import argparse
import os
import logging
import pandas as pd

vars = argparse.ArgumentParser(description='Connecting scoary output with emapper output')

vars.add_argument('--output_directory',dest='output_dir', type=str, required=True, help='absolute path to output directory must be given')

vars.add_argument('--benjaminiH_1173genomes',dest='benjaminiH_1173genomes', type=str, required=True, help='absolute path to the directory with all the bejamini hochberg corrected values for 1173 genomes')

vars.add_argument('--benjaminiH_400genomes',dest='benjaminiH_400genomes', type=str, required=True, help='absolute path to the directory with all the bejamini hochberg corrected values for 400 genomes')

args = vars.parse_args()

# create log file directory if it doesnt exist
if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
    os.makedirs(f"{os.path.expanduser('~')}/log_files")

# create output directories if they dont exist
if not os.path.isdir(f"{args.output_dir}"):
    os.makedirs(f"{args.output_dir}")

# configure the logger we will be using
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/combining_pvalues.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
logger = logging.getLogger(__name__)

# open files for Rplotting
with open(f"{args.output_dir}/combined_pvalue_frequency_plot.R",'w') as Rplotter_file, open(f"{args.output_dir}/unique_genes_both_correlations.txt",'w') as out_file:

    # write headers for out file
    out_file.write('Gene\tPan_genome_Size\n')

    bejaminiH_1173genomes = f"{args.benjaminiH_1173genomes}/gene_families_associated_genes_and_annotations_all_csv_columns.tsv"
    bejaminiH_400genomes = f"{args.benjaminiH_400genomes}/gene_families_associated_genes_and_annotations_all_csv_columns.tsv"

    # load tsvs into pandas df
    # put the unique genes into the file
    benjH_genes_1173 = pd.read_csv(filepath_or_buffer=bejaminiH_1173genomes,sep='\t')['Gene'].unique().tolist()
    benjH_genes_400 = pd.read_csv(filepath_or_buffer=bejaminiH_400genomes,sep='\t')['Gene'].unique().tolist()
    # bh_pos_genes = pd.read_csv(filepath_or_buffer=benjaminiH_positive_1173genomes,sep='\t')['Gene'].unique().tolist()
    # bh_neg_genes = pd.read_csv(filepath_or_buffer=benjaminiH_negative_1173genomes,sep='\t')['Gene'].unique().tolist()

    for gene in benjH_genes_1173:
        out_file.write(f"{gene}\t1173\n")

    for gene in benjH_genes_400:
        out_file.write(f"{gene}\t400\n")

    # write R file code
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
}\n'''

    # write installation code to file
    Rplotter_file.write(installation_chunk)

    # library packages we will be using
    Rplotter_file.write('library(ggplot2)\n')
    Rplotter_file.write('library(gridExtra)\n')
    Rplotter_file.write('library(grid)\n')

    # set a working directory
    Rplotter_file.write(f"setwd('{args.output_dir}')\n")

    # load in the different files
    Rplotter_file.write(f"master_df <- read.table(file='unique_genes_both_correlations.txt',header=TRUE,sep='\\t')\n")

    # change numerical values to factor to be plotted
    Rplotter_file.write(f"master_df$Pan_genome_Size <- as.factor(master_df$Pan_genome_Size)\n")
    # start graphics driver
    # load graphics driver
    Rplotter_file.write(f"pdf('{args.output_dir}/400_vs_1173_genomes_comparison.pdf', width = 16, height = 18)\n")

    # bonferroni positive vs negative
    # plot positively correlated genes. scale fill manual changes colour, labs for labelling, theme axis and theme panel removes the gird background and adjusts things to our liking
    Rplotter_file.write("ggplot(master_df, aes(x = Pan_genome_Size, fill = Pan_genome_Size)) +\n")
    Rplotter_file.write("geom_bar() +\n")
    Rplotter_file.write("scale_fill_manual(values=c('darkorange','darkblue')) +\n")
    Rplotter_file.write("labs(x = 'Pan-genome Size',y = 'Frequency of Correlated Genes') +\n")
    Rplotter_file.write("theme_minimal() +\n")
    Rplotter_file.write('theme(panel.background = element_rect(color="#FFFFFF", fill="white")) +\n')
    Rplotter_file.write('theme(panel.grid.minor = element_blank()) + \n')
    Rplotter_file.write('theme(panel.grid.major = element_blank()) + \n')
    Rplotter_file.write('''theme(axis.title.y = element_text(),axis.text = element_text(size = 30),
        axis.line.y = element_line(),
        axis.title=element_text(size = 40),
        axis.ticks.y = element_line(),
        axis.ticks.x = element_line(),
        axis.line.x = element_line())+\n''')
#    Rplotter_file.write('theme_bw(base_size = 17)+\n')
    Rplotter_file.write('theme(legend.position = "none")+\n')
    Rplotter_file.write('scale_y_continuous(expand = c(0, 0), limits = c(0, 21000))\n')

   # Rplotter_file.write('plot1 <- plot1 + scale_fill_discrete(name = "Pan-genome Size")\n')
    #Rplotter_file.write('plot1\n')
#    Rplotter_file.write('scale_y_continuous(labels = c("0","5000","10000","15000","20000"), breaks = c(0,5000,10000,15000,20000))\n')

    Rplotter_file.write("dev.off()\n")

# run the plotting script
# os.system(f"Rscript {args.output_dir}/stats_and_plots/cog_frequency_plots.R")
