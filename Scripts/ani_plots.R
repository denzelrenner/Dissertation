# install required packages 

# library packages we will be using
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)

# get the command line argument passed which should just be an output directory 
# followed by the name of the file with all the data inside of it
#args <- commandArgs(trailingOnly = TRUE)


# parse the input flags for the R script
#working_dir <- args[1]

#metrics_file <- args[2]

# set a working directory 
#setwd(working_dir)

# load in the file with all the metrics
ani_values <- read.table('ANI_values.txt',header = T)

# load graphics driver
pdf("ANI_distribution.pdf", width = 6, height = 4)

#ggplot(data = metrics_data,aes(x=Complete)) +
#  geom_histogram(color='black',fill='red') +
#  theme_classic()

# Reshape the data
# data_long <- pivot_longer(metrics_data, cols = c(Complete,Fragmented,Missing,Single_Copy), names_to = "metric", values_to = "value")

#title = "Complete (C)",
# create plot for distribution of complete buscos
ggplot(ani_values, aes(x = ANI_Values)) +
  geom_histogram(bins = 30, fill = "blue",color='black', alpha = 0.7) +
  #facet_wrap(~ metric, scales = "free_x") +
  labs( x = "%ANI", y = "Number of Pairwise Comparisons") +
  theme_minimal() +
  theme(panel.background = element_rect(color="#FFFFFF", fill="white")) +
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.title.y = element_text(),
        axis.text.y = element_text(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(),
        axis.line.x = element_line()) +
  #scale_y_continuous(expand = c(0, 0), breaks = seq(10, max(hist(metrics_data$Complete)$counts), by = 10))
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(labels = c("0.7","0.75","0.8","0.85","0.9","0.95","1"), breaks = c(0.7,0.75,0.8,0.85,0.9,0.95,1)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "orange")
#scale_x_continuous(expand = c(0, 0))

dev.off()


