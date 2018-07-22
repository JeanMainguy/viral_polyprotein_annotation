#!/usr/bin/env Rscript
library(ggplot2)
output_dir = './results/figures/'


file='data/clustering_result/Viruses/clustering_parameter_variation_homogeneity_evaluation/merge_Viruses_evalue_1e-20_I2_homogeneity.csv'
data = read.csv(file = file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$unclassified_cluster == "False",]

data$coverage =  paste(data$coverage, "%")

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=coverage, group=coverage, fill=coverage), alpha=0.4) + theme_minimal()+
theme(axis.text.x = element_text(size=13),
      axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
      legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold'))

output_file = paste(output_dir, 'homogeneity_density_coverage.png', sep ='')

  png(filename=output_file,  width = 1368, height = 768)
p
dev.off()
