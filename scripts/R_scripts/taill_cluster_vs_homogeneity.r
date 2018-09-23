#!/usr/bin/env Rscript
library(ggplot2)

file='results/clustering_evaluation/homogeneity_results_merge/merge_Viruses_evalue1e-30_I2_coverage_all_homogeneityonly_poly.csv'
file='results/clustering_evaluation/homogeneity_results_merge/merge_Viruses_evalue1e-30_I_all_coverage20_homogeneityonly_poly.csv'


file='data/clustering_result/Viruses/clustering_parameter_variation_homogeneity_evaluation/Viruses_evalue_1e-60coverage20_I1_4_homogeneity.csv'


data = read.csv(file = file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$unclassified_cluster == "False",]

coverage = 20
inflation="1.4"
evalue=1e-60

output_dir = "results/clustering_evaluation/"

data$Cluster_with_polyprotein = ifelse(data$nb_polyprotein == 0, 'False', 'True')


data$inflation = ifelse(grepl('_', data$inflation), gsub("_", ".", data$inflation), paste(data$inflation, '.0', sep=''))
data_poly = data[data$nb_polyprotein != 0,]
data_poly = data_poly[data_poly$inflation == inflation,]
data_poly = data_poly[data_poly$coverage == coverage,]
data_poly = data_poly[data_poly$evalue == evalue,]

caption_text = paste('\nClustering Parameters:\nInflation:',inflation, '   Evalue max:', evalue, '   coverage',coverage)

p = ggplot(data_poly, aes(y=nb_of_protein, x=homogeneity, size= nb_polyprotein), alpha=0.4) + 
  geom_point(alpha=0.5) +  
  labs( y="size of the cluster", size='number of\nannotated polyproteins', x="homogeneity", caption_text=caption_text)  + theme_minimal() 
#  geom_smooth(method=lm)
p



