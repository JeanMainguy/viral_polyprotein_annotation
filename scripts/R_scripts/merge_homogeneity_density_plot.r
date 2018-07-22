#!/usr/bin/env Rscript
library(ggplot2)
library(grid)
args = commandArgs(trailingOnly=TRUE)
# output_dir = './results/figures/'
# file='data/clustering_result/Viruses/clustering_parameter_variation_homogeneity_evaluation/merge_Viruses_evalue_1e-20_I2_homogeneity.csv'

# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
  print(length(args))
  print(args[1])
  stop("5 arguments must be supplied input file", call.=FALSE)
}

csv_file = args[1]
output_dir = args[2]
evalue= args[3]
inflation=args[4]
coverage=args[5]

csv_file='results/clustering_evaluation/homogeneity_results_merge/merge_Viruses_evalue1e-30_I2_coverage20_homogeneity.csv'
output_dir = 'results/clustering_evaluation/homogeneity_results_merge/'
evalue = 1e-30
inflation = 2
coverage = 20


capation_text = 'Clustering parameters:'
legend=''
fix_parameter = 0

data = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
summary(data)
data = data[data$unclassified_cluster == "False",]

data$Variable_Parameter = ''

if ( evalue != '*'){
  data = data[data$evalue == evalue,]
  output_dir = paste(output_dir, '_evalue', evalue, sep="")
  capation_text = paste(capation_text, "Evalue =", evalue)
  fix_parameter = fix_parameter +1
} else {
  data$Variable_Parameter =  paste(data$Variable_Parameter,"evalue", data$evalue)
  legend = paste(legend, 'Evalue')
}

if ( inflation != '*'){
  data = data[data$inflation == inflation,]
  output_dir = paste(output_dir, '_I', inflation, sep="")
  capation_text = paste(capation_text, "Inflation =", inflation)
    fix_parameter = fix_parameter +1
} else {
  data$Variable_Parameter =  paste(data$Variable_Parameter, 'inflation',data$inflation)
  legend = paste(legend, 'Inflation')
}

if ( coverage != '*'){
  data = data[data$coverage == coverage,]
  output_dir = paste(output_dir, '_coverage_all', sep="")
  capation_text = paste(capation_text, "Coverage =", coverage, '%')
  fix_parameter = fix_parameter +1
  print(data$coverage[5])
} else {
    data$Variable_Parameter =  paste(data$Variable_Parameter, "coverage",data$coverage, '%')
    legend = paste(legend, 'Coverage')
}

data$coverage =  paste(data$coverage, "%")
data$inflation =  paste('I=',data$inflation)
# data$evalue =  paste(data$evalue)



  #All

  print(fix_parameter)
  print(data$inflation[1:6])
  print(data$Variable_Parameter)
  p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=inflation, group=inflation, fill=inflation), alpha=0.4) +
  theme_minimal() +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_blank(),
        legend.text = element_text( size=20), legend.title = element_blank(),
        plot.caption = element_text(size=22))+
        labs(caption = capation_text)

  output_file = paste(output_dir, '_all.png', sep ='')

  png(filename=output_file,  width = 1500, height = 750)
  p
  dev.off()

  # Comparaison between Polyprotein clusters and non-polyprotein clusters
  data$Cluster_with_polyprotein = ifelse(data$nb_polyprotein == 0, 'True', 'False')
  data$Variable_Parameter_pol_vs_non =  paste(data$Variable_Parameter, data$Cluster_with_polyprotein)
  p = ggplot(data, aes(x=homogeneity)) +
        geom_density(aes(colour=Cluster_with_polyprotein,
                          group=Variable_Parameter_pol_vs_non,
                          fill=Cluster_with_polyprotein),
                          alpha=0.4) +
                        theme_minimal()+
                        theme(axis.text.x = element_text(size=20),
                              axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
                              legend.text = element_text( size=20), legend.title = element_text(size=22, face='bold'),
                              plot.caption = element_text(size=22))+
                              labs(caption = capation_text)

  output_file = paste(output_dir, '_poly_vs_non-poly.png', sep ='')


  png(filename=output_file,  width = 1500, height = 750)
  p
  dev.off()

  summary(data)
  print(fix_parameter)
  #Only Polyprotein cluster
  data = data[data$nb_polyprotein != 0,] # take only polyprotein

  p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=Variable_Parameter, group=Variable_Parameter, fill=Variable_Parameter), alpha=0.4) +
  theme_minimal()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_blank(),
        plot.caption = element_text(size=22))+
        labs(caption = capation_text)

  output_file = paste(output_dir, '_polyprotein_only.png', sep ='')

  png(filename=output_file,  width = 1500, height = 750)
  p
  dev.off()
