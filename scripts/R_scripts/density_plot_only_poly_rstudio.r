#!/usr/bin/env Rscript
library(ggplot2)
library(grid)

csv_fileI='results/clustering_evaluation/homogeneity_results_merge/merge_Viruses_evalue1e-30_I_all_coverage20_homogeneityonly_poly.csv'
csv_fileC='results/clustering_evaluation/homogeneity_results_merge/merge_Viruses_evalue1e-30_I2_coverage_all_homogeneityonly_poly.csv'
csv_fileE="results/clustering_evaluation/homogeneity_results_merge/merge_Viruses_evalue_all_I2_coverage20_homogeneityonly_poly.csv"
  
evalue = 1e-30
inflation = 2
coverage = 20


output_dir = 'results/clustering_evaluation/homogeneity_results_merge/'

#EVALUE variation
capation_text = '\nClustering parameters:\n'

data = read.csv(file = csv_fileE, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$unclassified_cluster == "False",]

data = data[data$inflation == inflation,]
output_dir = paste(output_dir, '_I', inflation, sep="")
capation_text = paste(capation_text, "Inflation=", inflation, sep='')

data = data[data$coverage == coverage,]
output_dir = paste(output_dir, '_coverage_all', sep="")
capation_text = paste(capation_text, " Coverage=", coverage, '%', sep='')

#Filter some of Evalue to display less data
data = data[data$evalue != 1e-100,]
data = data[data$evalue != 1e-120,]
data = data[data$evalue != 1e-70,]
data = data[data$evalue != 1e-40,]
data = data[data$evalue != 1e-20,]
data = data[data$evalue != 1e-60,]
data = data[data$evalue != 1e-10,]

data$evalue_string = paste('Evalue=',data$evalue)
data = data[data$nb_polyprotein != 0,] # take only polyprotein

p = ggplot(data, aes(x=homogeneity)) + 
  geom_density(aes(colour=reorder(evalue, evalue), group=reorder(evalue, evalue), fill=reorder(evalue, evalue)), alpha=0.1) +
  theme_minimal()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_text(size=22,face="bold"),
        plot.caption = element_text(size=22, hjust=0))+
  labs(caption = capation_text) +
  scale_fill_discrete(name="Evalue") +
  scale_color_discrete(name="Evalue")
p
output_file = paste(output_dir, '_polyprotein_only.png', sep ='')

png(filename=output_file,  width = 1500, height = 750)
p
dev.off()

#COVERAGE variation
capation_text = '\nClustering parameters:\n'

data = read.csv(file = csv_fileC, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$unclassified_cluster == "False",]

data = data[data$inflation == inflation,]
output_dir = paste(output_dir, '_I', inflation, sep="")
capation_text = paste(capation_text, "Inflation=", inflation, sep='')

data = data[data$evalue == evalue,]
output_dir = paste(output_dir, '_evalue', evalue, sep="")
capation_text = paste(capation_text, " Evalue =", evalue, sep='')

#Filter some of  to display less data
data = data[data$coverage != 10,]

data$coverage = paste(data$coverage, '%')
data = data[data$nb_polyprotein != 0,] # take only polyprotein

p = ggplot(data, aes(x=homogeneity)) + 
  geom_density(aes(colour=coverage, group=coverage, fill=coverage), alpha=0.1) +
  theme_minimal()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_text(size=22,face="bold"),
        plot.caption = element_text(hjust=0, size=22))+
  labs(caption = capation_text) +
  scale_fill_discrete(name="Coverage\nthreshold") +
  scale_color_discrete(name="Coverage\nthreshold")
p
output_file = paste(output_dir, '_polyprotein_only.png', sep ='')

png(filename=output_file,  width = 1500, height = 750)
p
dev.off()


#INFLATION variation
capation_text = '\nClustering parameters:\n'

data = read.csv(file = csv_fileI, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$unclassified_cluster == "False",]

data = data[data$coverage == coverage,]
output_dir = paste(output_dir, '_coverage_all', sep="")
capation_text = paste(capation_text, " Coverage=", coverage, '%', sep='')

data = data[data$evalue == evalue,]
output_dir = paste(output_dir, '_evalue', evalue, sep="")
capation_text = paste(capation_text, " Evalue =", evalue, sep='')

data$inflation = gsub("_", ".", data$inflation)

#Filter some of  to display less data
data = data[data$inflation != '1.8',]

data$inflation = paste('I=',data$inflation)
data = data[data$nb_polyprotein != 0,] # take only polyprotein

p = ggplot(data, aes(x=homogeneity)) + 
  geom_density(aes(colour=inflation, group=inflation, fill=inflation), alpha=0.1) +
  theme_minimal()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_text(size=22,face="bold"),
        plot.caption = element_text(hjust=0, size=22))+
  labs(caption = capation_text) +
  scale_fill_discrete(name="Inflation") +
  scale_color_discrete(name="Inflation")
p
output_file = paste(output_dir, '_polyprotein_only.png', sep ='')

print(output_file)
png(filename=output_file,  width = 1500, height = 750)
p
dev.off()