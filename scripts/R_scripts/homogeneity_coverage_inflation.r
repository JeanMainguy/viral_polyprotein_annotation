#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  else if (length(args)==1) {
  # default output file base
  args[2] = "test/test_plot"
}
csv_file = args[1]
output_dir = args[2]


data = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# data = data[data$first_node == "Retro-transcribing viruses",]
data = data[data$unclassified_cluster == "False",]
# data = data[data$evalue == "1e-20",]
# data = data[ data$coverage < 90,]
# data = data[ data$coverage > 0,]
# data = data[ data$inflation != 6,]
# data = data[ data$coverage != 60,]
# data = data[ data$coverage != 10,]
# data = data[ data$coverage != 20,]
# data = data[ data$coverage != 30,]
# data = data[ data$coverage != 70,]
# data = data[ data$coverage != 80,]
# data = data[ data$coverage == 40,]
# data = data[ data$inflation != 3,]
# data = data[ data$inflation != 4,]
data = data[ data$inflation != "1",]
data = data[ data$inflation != "7",]
# data = data[ data$inflation != "5",]
data$inflation_coverage <- paste(data$coverage,data$inflation)
data$coverage <- paste(data$coverage,"%")

# print(data$inflation_coverage)

# p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(group=inflation_coverage, colour=inflation)) + theme_minimal()
#
# output_file = paste(output_dir, 'homogeneity_density_inflation.png', sep ='')
#
# png(filename=output_file,  width = 1368, height = 768)
# p
# dev.off()
#
# p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(group=inflation_coverage, colour=inflation_coverage)) + theme_minimal()
#
# output_file = paste(output_dir, 'homogeneity_density_inflation_coverage.png', sep ='')
#
# png(filename=output_file,  width = 1368, height = 768)
# p
# dev.off()
#
# p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(group=inflation_coverage, colour=coverage)) + theme_minimal()
#
# output_file = paste(output_dir, 'homogeneity_density_coverage.png', sep ='')
#
# png(filename=output_file,  width = 1368, height = 768)
# p
# dev.off()
#We don't display the name of first node that appears less than n time
#They are labeled as Other
# n = 100
# first_nodes = table(data$first_node)
# first_nodes = data.frame(first_nodes, stringsAsFactors = FALSE)
# colnames(first_nodes) = c('First_node', 'Freq')
# other_nodes = first_nodes$First_node[first_nodes$Freq < n]
# data$first_node[data$first_node %in% other_nodes] = "Other"

#ALL
data$coverage_evalue =  paste(data$coverage,data$evalue)
p <- ggplot(data, aes(x=coverage_evalue, y=homogeneity, fill=inflation)) +
  geom_boxplot()  + theme_minimal() + scale_fill_brewer(palette="Set3")#+  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")


  output_file = paste(output_dir, 'box_plot_homogeneity_inflation_coverage.png', sep ='')
  png(filename=output_file,  width = 1368, height = 768)
  p
  dev.off()

#POLYPROT AND NON POLYPROT CLUSTER
data$polyprotein = data$polyprotein = ifelse(data$nb_polyprotein == 0, '', 'Polyprotein')
data$coverage_polyprot =  paste(data$coverage,data$polyprotein)

p <- ggplot(data, aes(x=coverage_polyprot, y=homogeneity, fill=inflation)) +
  geom_boxplot()  + theme_minimal() + scale_fill_brewer(palette="Set3")#+  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")


  output_file = paste(output_dir, 'box_plot_homogeneity_inflation_coverage_polyprot.png', sep ='')
  png(filename=output_file,  width = 1368, height = 768)
  p
  dev.off()

excluded_node='dsDNA viruses'
excluded_node="Caulimoviridae"
with=paste(excluded_node)
without=paste('without', excluded_node)

  data$includedsDNAViruses = ifelse(grepl(excluded_node, data$valid_branches, fixed=TRUE), with, without)
  data$includedsDNAViruses = ifelse(grepl(excluded_node, data$shared_taxonomy, fixed=TRUE), with, without)

  data_no_dsDNA = data[data$includedsDNAViruses == FALSE,]

#WITH AND WITHOUT EXCLUDED NODE
data$coverage_excluded_node =  paste(data$coverage,data$includedsDNAViruses)

  p <- ggplot(data, aes(x=coverage_excluded_node, y=homogeneity, fill=inflation)) +
    geom_boxplot()  + theme_minimal() + scale_fill_brewer(palette="Set3")+
  theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number of domain annotations", fill="Interpro Database:", x="Overlap distance (amino acid)") #+  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")


    output_file = paste(output_dir, 'box_plot_homogeneity_inflation_coverage_excludedNode.png', sep ='')
    png(filename=output_file,  width = 1368, height = 768)
    p
    dev.off()

  # data$includedsDNAViruses = ifelse(grepl('dsDNA viruses', data$valid_branches, fixed=TRUE), TRUE, FALSE)
  # data$includedsDNAViruses = ifelse(grepl('dsDNA viruses', data$shared_taxonomy, fixed=TRUE), TRUE, FALSE)
  # data$include_dsDNAViruses_and_Inflation <- paste(data$includedsDNAViruses,data$inflation)
  # data_no_dsDNA = data[data$includedsDNAViruses == FALSE,]


  # p <- ggplot(data, aes(x=inflation_coverage, y=homogeneity, fill=coverage)) +
  #   geom_boxplot(position=position_dodge(1))  + theme_minimal()+ stat_summary(fun.y=mean, geom="point", shape=19, size=4)
  #
  #
  #   output_file = paste(output_dir, 'box_plot_homogeneity_coverage.png', sep ='')
  #   png(filename=output_file,  width = 1368, height = 768)
  #   p
  #   dev.off()
  #
  # p <- ggplot(data, aes(x=inflation_coverage, y=homogeneity, fill=inflation)) +
  #   geom_boxplot(position=position_dodge(1))  + theme_minimal() + stat_summary(fun.y=mean, geom="point", shape=19, size=4)
  #
  #
  #   output_file = paste(output_dir, 'box_plot_homogeneity_inflation.png', sep ='')
  #   png(filename=output_file,  width = 1368, height = 768)
  #   p
  #   dev.off()
