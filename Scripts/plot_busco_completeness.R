# install needed packages

# some of the code here was edited and adapted from the busco plot

# library packages we will be using
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)

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
pdf("busco_completeness.pdf", width = 6, height = 4)

#ggplot(data = metrics_data,aes(x=Complete)) +
#  geom_histogram(color='black',fill='red') +
#  theme_classic()

# Reshape the data
data_long <- pivot_longer(metrics_data, cols = c(Complete,Fragmented,Missing,Single_Copy), names_to = "metric", values_to = "value")

# Create the faceted plot
#ggplot(data_long, aes(x = value)) +
#  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
#  facet_wrap(~ metric, scales = "free_x") +
#  labs(title = "Faceted Histograms of Metrics", x = "Value", y = "Count") +
#  theme_minimal()

#title = "Complete (C)",
# create plot for distribution of complete buscos
plot1 <- ggplot(metrics_data, aes(x = Complete)) +
  geom_histogram(bins = 30, fill = "blue",color='black', alpha = 0.7) +
  #facet_wrap(~ metric, scales = "free_x") +
  labs( x = "%Complete", y = "Number of Genomes") +
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
  scale_x_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100))
  #scale_x_continuous(expand = c(0, 0))

#plot1+geom_vline(xintercept = 90, linetype = "dashed", color = "orange")

plot2 <- ggplot(metrics_data, aes(x = Missing)) +
  geom_histogram(bins = 30, fill = "purple",color='black', alpha = 0.7) +
  #facet_wrap(~ metric, scales = "free_x") +
  labs(x = "%Missing", y = "Number of Genomes") +
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
  scale_x_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100))
#scale_x_continuous(expand = c(0, 0))

plot3 <- ggplot(metrics_data, aes(x = Fragmented)) +
  geom_histogram(bins = 30, fill = "purple",color='black', alpha = 0.7) +
  #facet_wrap(~ metric, scales = "free_x") +
  labs(x = "%Fragmented", y = "Number of Genomes") +
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
  scale_x_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100))
#scale_x_continuous(expand = c(0, 0))

plot4 <- ggplot(metrics_data, aes(x = Total_Length)) +
  geom_histogram(bins = 30, fill = "purple",color='black', alpha = 0.7) +
  #facet_wrap(~ metric, scales = "free_x") +
  labs(x = "Assembly Length (in millions)", y = "Number of Genomes") +
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
  scale_y_continuous(expand = c(0, 0))
  #scale_x_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100))
#scale_x_continuous(expand = c(0, 0))

grid.arrange(plot1+geom_vline(xintercept = 90, linetype = "dashed", color = "orange"),plot2,plot3,plot4,ncol=2,nrow=2,top = textGrob("BUSCO Assessment Results", gp = gpar(fontsize = 16)))

dev.off()

#axis.title.y = element_text(),  # Ensure y-axis title is not blank
#axis.text.y = element_text(),   # Ensure y-axis labels are not blank
#axis.line.y = element_line()    # Ensure y-axis line is visible
#)
