#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# output_dir = './results/figures/'
# file='data/clustering_result/Viruses/clustering_parameter_variation_homogeneity_evaluation/merge_Viruses_evalue_1e-20_I2_homogeneity.csv'

# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
  print(length(args))
  print(args[1])
  stop("5 arguments must be supplied input file", call.=FALSE)
}
#  else if (length(args)==1) {
# # default arguments
#   args[2] = "tmp/$USER/homogeneity_density"
#   args[1] = 1e-20
#   args[2] = 2
#   args[3] = "*"
#
# }
csv_file = args[1]
output_dir = args[2]
evalue= args[3]
inflation=args[4]
coverage=args[5]


data = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
summary(data)
data = data[data$unclassified_cluster == "False",]

data_comparaison = data

if ( evalue != '*'){
  data = data[data$evalue == evalue,]
  output_dir = paste(output_dir, '_evalue', evalue, sep="")
} else {
  data_comparaison = data_comparaison[data_comparaison$evalue == 1e-20,]
}

if ( inflation != '*'){
  data = data[data$inflation == inflation,]
  output_dir = paste(output_dir, '_I', inflation, sep="")
} else {
  data_comparaison = data_comparaison[data_comparaison$inflation == 2,]
}

if ( coverage != '*'){
  data = data[data$coverage == coverage,]
  output_dir = paste(output_dir, '_coverage_all', sep="")
} else {
  data_comparaison = data_comparaison[data_comparaison$coverage == 20,]
}

data$coverage =  paste(data$coverage, "%")
data$inflation =  paste('I=',data$inflation)
# data$evalue =  paste(data$evalue)

#All
p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=coverage, group=coverage, fill=coverage), alpha=0.4) + theme_minimal()+
theme(axis.text.x = element_text(size=13),
      axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
      legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold'))

output_file = paste(output_dir, '_all.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()

#Only Polyprotein cluster
data = data[data$nb_polyprotein != 0,] # take only polyprotein

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=coverage, group=coverage, fill=coverage), alpha=0.4) + theme_minimal()+
theme(axis.text.x = element_text(size=13),
      axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
      legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold'))

output_file = paste(output_dir, '_polyprotein_only.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()


# Comparaison between Polyprotein clusters and non-polyprotein clusters
data_comparaison$polyprotein_cluster = ifelse(data_comparaison$nb_polyprotein == 0, 'False', 'True')
p = ggplot(data_comparaison, aes(x=homogeneity)) +
  geom_density(aes(colour=polyprotein_cluster,
                    group=polyprotein_cluster,
                    fill=polyprotein_cluster),
                    alpha=0.4,
                    size=2) +
theme_minimal()+
theme(axis.text.x = element_text(size=13),
      axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
      legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold'))

output_file = paste(output_dir, '_poly_vs_non-poly.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()
