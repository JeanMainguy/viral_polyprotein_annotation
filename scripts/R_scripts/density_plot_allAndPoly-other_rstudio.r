#!/usr/bin/env Rscript
library(ggplot2)
library(grid)

csv_file='results/clustering_evaluation/homogeneity_results_merge/merge_Viruses_evalue1e-30_I2_coverage20_homogeneity.csv'
output_dir = 'results/clustering_evaluation/homogeneity_results_merge/'
evalue = 1e-30
inflation = 2
coverage = 20


capation_text = '\nClustering parameters:\n'
legend=''

data = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

data = data[data$unclassified_cluster == "False",]


  data = data[data$evalue == evalue,]
  output_dir = paste(output_dir, '_evalue', evalue, sep="")
  capation_text = paste(capation_text, "Evalue=", evalue, sep='')


  data = data[data$inflation == inflation,]
  output_dir = paste(output_dir, '_I', inflation, sep="")
  capation_text = paste(capation_text, " Inflation=", inflation,sep='')



  data = data[data$coverage == coverage,]
  output_dir = paste(output_dir, '_coverage_all', sep="")
  capation_text = paste(capation_text, " Coverage=", coverage, '%', sep='')
  fix_parameter = fix_parameter +1

data$coverage =  paste(data$coverage, "%")
data$inflation =  paste('I=',data$inflation)



  print(fix_parameter)
  print(data$inflation[1:6])
  print(data$Variable_Parameter)
  p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=inflation, group=inflation, fill=inflation), alpha=0.4) +
  theme_minimal() +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_blank(), legend.title = element_blank(),
        plot.caption = element_text(hjust=0, size=22))+
        labs(caption = capation_text)

  output_file = paste(output_dir, '_all.png', sep ='')
  p
  
  png(filename=output_file,  width = 1083, height = 750)
  p
  dev.off()

  # Comparaison between Polyprotein clusters and non-polyprotein clusters
  data$Cluster_with_polyprotein = ifelse(data$nb_polyprotein == 0, 'False', 'True')
  
  p = ggplot(data, aes(x=homogeneity)) +
        geom_density(aes(colour=Cluster_with_polyprotein,
                          group=Cluster_with_polyprotein,
                          fill=Cluster_with_polyprotein),
                          alpha=0.4) +
                        theme_minimal()+
    
                        theme(axis.text.x = element_text(size=20),
                              axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
                              legend.text = element_text( size=20), legend.title = element_text(size=22, face='bold'),
                              plot.caption = element_text(hjust=0, size=22), legend.position="right")+
                              labs(caption = capation_text)   +
                              scale_fill_discrete(name="Cluster with \npolyprotein") +
                             scale_color_discrete(name="Cluster with \npolyprotein")

  p
  output_file = paste(output_dir, '_poly_vs_non-poly.png', sep ='')
  print(output_file)

  png(filename=output_file,  width = 1083, height = 750)
  p
  dev.off()
  p

