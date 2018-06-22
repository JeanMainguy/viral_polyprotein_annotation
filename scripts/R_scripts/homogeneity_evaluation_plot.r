#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  else if (length(args)==1) {
  # default output file base
  args[2] = "plot"
}
csv_file = args[1]
basename = args[2]
output_dir = basename# paste("results/figures", basename, sep = "/")

data = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$unclassified_cluster == "False",] 
# Density plot

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(group=inflation, colour=inflation), alpha=0) + theme_minimal()

output_file = paste(output_dir, 'homogeneity_density_inflation.png', sep ='')

png(filename=output_file,  width = 900, height = 500)
p
dev.off()


#COLOR BY DOMAIN ANNOTATIONS
#We don't display the name of first node that appears less than  time
#They are labeled as Other
first_nodes = table(data$first_node)
first_nodes = data.frame(first_nodes, stringsAsFactors = FALSE)
colnames(first_nodes) = c('First_node', 'Freq')
other_nodes = first_nodes$First_node[first_nodes$Freq < 50]
data$first_node[data$first_node %in% other_nodes] = "Other"

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=first_node, group=first_node), alpha=0) + theme_minimal()

output_file = paste(output_dir, 'homogeneity_density_first_node.png', sep ='')

png(filename=output_file,  width = 1800, height = 600)
p
dev.off()


data$polyprotein = ifelse(data$nb_polyprotein == 0, FALSE, TRUE)
p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=polyprotein, group=polyprotein), alpha=0) + theme_minimal()

output_file = paste(output_dir, 'homogeneity_density_polyprot.png', sep ='')

png(filename=output_file,  width = 1800, height = 600)
p
dev.off()
