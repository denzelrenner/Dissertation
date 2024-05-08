# install packages you need
######

# library packages we will be using
library(ggplot2)

# get the command line argument passed which should just be an output directory 
# followed by the name of the file with all the data inside of it
args <- commandArgs(trailingOnly = TRUE)

# parse the input flags for the R script
working_dir <- args[1]

metrics_file <- args[2]

# set a working directory 
setwd(working_dir)

# load in the file with all the metrics
metrics_data <- read.table(metrics_file,header = T)

# load graphics driver
pdf("assembly_length_distribution.pdf", width = 6, height = 4)

# plot the assembly lengths for all the different species
ggplot(data = metrics_data,aes(x=Species_Name,y=Assembly_Length)) +
  geom_point()

dev.off()

# load graphics driver
pdf("gc_content_distribution.pdf", width = 6, height = 4)

# plot the assembly gc content for all the different species
ggplot(data = metrics_data,aes(x=Species_Name,y=GC_Percentage)) +
  geom_point()

dev.off()
