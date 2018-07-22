#!/usr/bin/env Rscript
library(ggplot2)
output_dir = './'

file="data/clustering_result/Viruses/homogeneity_evalue_coverage_inflation/Viruses_evalue_1e-20coverage40_I2_tax_homogeneity.csv"
# file="data/clustering_result/Viruses/inflation_test_no_ceil_200/Viruses_evalue_1e-20coverage40_I2_tax_homogeneity.csv"

data = read.csv(file = file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$unclassified_cluster == "False",]

#COLOR BY FISR NODE
#We don't display the name of first node that appears less than n time
#They are labeled as Other
first_nodes = table(data$first_node)
first_nodes = data.frame(first_nodes, stringsAsFactors = FALSE)
colnames(first_nodes) = c('First_node', 'Freq')
other_nodes = first_nodes$First_node[first_nodes$Freq < 20]
data$first_node[data$first_node %in% other_nodes] = "Other"

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=first_node, group=first_node, fill=first_node), alpha=0.4) + theme_minimal()+
theme(axis.text.x = element_text(size=13),
      axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
      legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold'))

output_file = paste(output_dir, 'homogeneity_density_first_node.png', sep ='')

  png(filename=output_file,  width = 1368, height = 768)
p
dev.off()

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(), alpha=1) + theme_minimal()+
theme(axis.text.x = element_text(size=13),
      axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
      legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
labs( y="Density", x="Evalue")


output_file = paste(output_dir, 'homogeneity_density_inflation.png', sep ='')


  png(filename=output_file,  width = 1368, height = 768)
p
dev.off()

data$polyprotein = ifelse(data$nb_polyprotein == 0, FALSE, TRUE)
p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=polyprotein, group=polyprotein, fill=polyprotein), alpha=0.4) + theme_minimal()+
theme(axis.text.x = element_text(size=13),
      axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
      legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold'))

output_file = paste(output_dir, 'density_polyprot.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()
