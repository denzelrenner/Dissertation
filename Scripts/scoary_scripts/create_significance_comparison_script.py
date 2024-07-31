import argparse
import os
import logging
import pandas as pd

vars = argparse.ArgumentParser(description='Connecting scoary output with emapper output')

vars.add_argument('--output_directory',dest='output_dir', type=str, required=True, help='absolute path to output directory must be given')

vars.add_argument('--benjaminiH_1173genomes',dest='benjaminiH_1173genomes', type=str, required=True, help='absolute path to the directory with all the bejamini hochberg corrected values for 1173 genomes')

vars.add_argument('--bonferroni_1173genomes',dest='bonferroni_1173genomes', type=str, required=True,help='absolute path to the directory with all the bonferroni corrected values for 1173 genomes')

#vars.add_argument('--benjaminiH_400genomes',dest='benjaminiH_400genomes', type=str, required=True, help='absolute path to the directory with all the bejamini hochberg corrected values for 400 genomes')

#vars.add_argument('--bonferroni_400genomes',dest='bonferroni_400genomes', type=str, required=True,help='absolute path to the directory with all the bonferroni corrected values for 400 genomes')

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
with open(f"{args.output_dir}/combined_pvalue_frequency_plot.R",'w') as Rplotter_file, open(f"{args.output_dir}/pos_unique_genes_both_correlations.txt",'w') as pos_out_file,open(f"{args.output_dir}/neg_unique_genes_both_correlations.txt",'w') as neg_out_file:

    # write headers for out file
    pos_out_file.write('Gene\tCorrelation\tCorrection_Method\n')
    neg_out_file.write('Gene\tCorrelation\tCorrection_Method\n')

    # set path for positive and negative correlated gene files
#    bonferroni_positive_400genomes = f"{args.bonferroni_400genomes}/master_positive_correlated_genes.tsv"
#    bonferroni_negative_400genomes = f"{args.bonferroni_400genomes}/master_negative_correlated_genes.tsv"
#    benjaminiH_positive_400genomes = f"{args.benjaminiH_400genomes}/master_positive_correlated_genes.tsv"
#    benjaminiH_negative_400genomes = f"{args.benjaminiH_400genomes}/master_negative_correlated_genes.tsv"

    bonferroni_positive_1173genomes = f"{args.bonferroni_1173genomes}/master_positive_correlated_genes.tsv"
    bonferroni_negative_1173genomes = f"{args.bonferroni_1173genomes}/master_negative_correlated_genes.tsv"
    benjaminiH_positive_1173genomes = f"{args.benjaminiH_1173genomes}/master_positive_correlated_genes.tsv"
    benjaminiH_negative_1173genomes = f"{args.benjaminiH_1173genomes}/master_negative_correlated_genes.tsv"

    # load tsvs into pandas df
    # put the unique genes into the file
    bon_pos_genes = pd.read_csv(filepath_or_buffer=bonferroni_positive_1173genomes,sep='\t')['Gene'].unique().tolist()
    bon_neg_genes = pd.read_csv(filepath_or_buffer=bonferroni_negative_1173genomes,sep='\t')['Gene'].unique().tolist()
    bh_pos_genes = pd.read_csv(filepath_or_buffer=benjaminiH_positive_1173genomes,sep='\t')['Gene'].unique().tolist()
    bh_neg_genes = pd.read_csv(filepath_or_buffer=benjaminiH_negative_1173genomes,sep='\t')['Gene'].unique().tolist()

    for gene in bon_pos_genes:
        pos_out_file.write(f"{gene}\tpositive\tBonferroni\n")

    for gene in bon_neg_genes:
        neg_out_file.write(f"{gene}\tnegative\tBonferroni\n")

    for gene in bh_pos_genes:
        pos_out_file.write(f"{gene}\tpositive\tBenjaminiH\n")

    for gene in bh_neg_genes:
        neg_out_file.write(f"{gene}\tnegative\tBenjaminiH\n")

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
   # Rplotter_file.write(f"n400genomes_bonferroni_positive <- read.table(file='{bonferroni_positive_400genomes}',header=TRUE,sep='\\t')\n")
   # Rplotter_file.write(f"n400genomes_bonferroni_negative <- read.table(file='{bonferroni_negative_400genomes}',header=TRUE,sep='\\t')\n")
   # Rplotter_file.write(f"n400genomes_benjaminiH_positive <- read.table(file='{benjaminiH_positive_400genomes}',header=TRUE,sep='\\t')\n")
   # Rplotter_file.write(f"n400genomes_benjaminiH_negative <- read.table(file='{benjaminiH_negative_400genomes}',header=TRUE,sep='\\t')\n")

   # Rplotter_file.write(f"n1173genomes_bonferroni_positive <- read.table(file='{bonferroni_positive_1173genomes}',header=TRUE,sep='\\t')\n")
   # Rplotter_file.write(f"n1173genomes_bonferroni_negative <- read.table(file='{bonferroni_negative_1173genomes}',header=TRUE,sep='\\t')\n")
   # Rplotter_file.write(f"n1173genomes_benjaminiH_positive <- read.table(file='{benjaminiH_positive_1173genomes}',header=TRUE,sep='\\t')\n")
   # Rplotter_file.write(f"n1173genomes_benjaminiH_negative <- read.table(file='{benjaminiH_negative_1173genomes}',header=TRUE,sep='\\t')\n")

    # depending on which genomes you are looking at, add a column for the genome size
   # Rplotter_file.write(f"n400genomes_bonferroni_positive$Genome_Size <- '400'\n")
   # Rplotter_file.write(f"n400genomes_bonferroni_negative$Genome_Size <- '400'\n")
   # Rplotter_file.write(f"n400genomes_benjaminiH_positive$Genome_Size <- '400'\n")
   # Rplotter_file.write(f"n400genomes_benjaminiH_negative$Genome_Size <- '400'\n")

   # Rplotter_file.write(f"n1173genomes_bonferroni_positive$Genome_Size <- '1173'\n")
   # Rplotter_file.write(f"n1173genomes_bonferroni_negative$Genome_Size <- '1173'\n")
   # Rplotter_file.write(f"n1173genomes_benjaminiH_positive$Genome_Size <- '1173'\n")
   # Rplotter_file.write(f"n1173genomes_benjaminiH_negative$Genome_Size <- '1173'\n")

    # combine the positive and negaitve dataset between the two
   # Rplotter_file.write(f"combined_bonferroni_positive <- rbind(n400genomes_bonferroni_positive,n1173genomes_bonferroni_positive)\n")
   # Rplotter_file.write(f"combined_bonferroni_negative <- rbind(n400genomes_bonferroni_negative,n1173genomes_bonferroni_negative)\n")
   # Rplotter_file.write(f"combined_benjaminiH_positive <- rbind(n400genomes_benjaminiH_positive,n1173genomes_benjaminiH_positive)\n")
   # Rplotter_file.write(f"combined_benjaminiH_negative <- rbind(n400genomes_benjaminiH_negative,n1173genomes_benjaminiH_negative)\n")


    # start graphics driver
    # load graphics driver
    #Rplotter_file.write(f"pdf('{args.output_dir}/pvalue_distribution.pdf', width = 12, height = 12)\n")

    # bonferroni positive 400 vs 1173
    # plot positively correlated genes. scale fill gradient changes colour based on frequency, labs for labelling, theme axis so axis marks and text look better
    #Rplotter_file.write("plot1 <- ggplot(combined_bonferroni_positive, aes(x = Genome_Size, fill = Genome_Size)) +\n")
    #Rplotter_file.write("geom_bar() +\n")
    # Rplotter_file.write("scale_fill_gradient(low = 'lightblue', high = 'darkblue') +\n")
   # Rplotter_file.write("labs(title = 'Pvalue comparison Bonferroni Positive',x = 'Genome Size',y = 'Positively Correlated Genes') +\n")
    #Rplotter_file.write("theme_minimal() \n")
    # Rplotter_file.write("theme(axis.text.x = element_text(angle = 45, hjust = 1))\n")

    # bonferroni negative 400 vs 1173
    #Rplotter_file.write("plot2 <- ggplot(combined_bonferroni_negative, aes(x = Genome_Size, fill = Genome_Size)) +\n")
    #Rplotter_file.write("geom_bar() +\n")
    #Rplotter_file.write("labs(title = 'Pvalue comparison Bonferroni Negative',x = 'Genome Size',y = 'Negatively Correlated Genes') +\n")
    #Rplotter_file.write("theme_minimal() \n")

    # benjaminiH positive 400 vs 1173
    #Rplotter_file.write("plot3 <- ggplot(combined_benjaminiH_positive, aes(x = Genome_Size, fill = Genome_Size)) +\n")
    #Rplotter_file.write("geom_bar() +\n")
    #Rplotter_file.write("labs(title = 'Pvalue comparison BenjaminiH Positive',x = 'Genome Size',y = 'Positively Correlated Genes') +\n")
    #Rplotter_file.write("theme_minimal() \n")

    # benjaminiH negative 400 vs 1173
    #Rplotter_file.write("plot4 <- ggplot(combined_benjaminiH_negative, aes(x = Genome_Size, fill = Genome_Size)) +\n")
    #Rplotter_file.write("geom_bar() +\n")
    #Rplotter_file.write("labs(title = 'Pvalue comparison BenjaminiH Negative',x = 'Genome Size',y = 'Negatively Correlated Genes') +\n")
    #Rplotter_file.write("theme_minimal() \n")

    # arrange the plots on the same page
    #Rplotter_file.write("grid.arrange(plot1,plot2,plot3,plot4,ncol=2,nrow=2)\n")

    # save plot
    #Rplotter_file.write("dev.off()\n")

    Rplotter_file.write(f"master_df_pos <- read.table(file='pos_unique_genes_both_correlations.txt',header=TRUE,sep='\\t')\n")
    Rplotter_file.write(f"master_df_neg <- read.table(file='neg_unique_genes_both_correlations.txt',header=TRUE,sep='\\t')\n")

    # start graphics driver
    # load graphics driver
    Rplotter_file.write(f"pdf('{args.output_dir}/number_of_genes_comparison.pdf', width = 16, height = 18)\n")

    # bonferroni positive vs negative
    # plot positively correlated genes. scale fill manual changes colour, labs for labelling, theme axis and theme panel removes the gird background and adjusts things to our liking
    Rplotter_file.write("plot1 <- ggplot(master_df_pos, aes(x = Correction_Method, fill = Correction_Method)) +\n")
    Rplotter_file.write("geom_bar() +\n")
    Rplotter_file.write("scale_fill_manual(values=c('darkorange','darkblue')) +\n")
    Rplotter_file.write("labs(x = 'Correction Method',y = 'Number of Positively Correlated Genes') +\n")
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
#    Rplotter_file.write('scale_y_continuous(labels = c("0","5000","10000","15000","20000"), breaks = c(0,5000,10000,15000,20000))\n')

    Rplotter_file.write("plot2 <- ggplot(master_df_neg, aes(x = Correction_Method, fill = Correction_Method)) +\n")
    Rplotter_file.write("geom_bar() +\n")
    Rplotter_file.write("scale_fill_manual(values=c('darkorange','darkblue')) +\n")
    Rplotter_file.write("labs(x = 'Correction Method',y = 'Number of Negatively Correlated Genes') +\n")
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
#    Rplotter_file.write('scale_y_continuous(labels = c("0","5000","10000","15000","20000"), breaks = c(0,5000,10000,15000,20000))\n')
    Rplotter_file.write('grid.arrange(plot1,plot2,ncol=2)\n')

    Rplotter_file.write("dev.off()\n")

# run the plotting script
# os.system(f"Rscript {args.output_dir}/stats_and_plots/cog_frequency_plots.R")

