#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# function for number of observations
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

# function for mean labels
mean.n <- function(x){
  return(c(y = median(x)*0.97, label = round(mean(x),2)))
  # experiment with the multiplier to find the perfect position
}

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  else if (length(args)==1) {
  # default output file base
  args[2] = "test/test_plot"
}
csv_file = args[1]
output_dir = args[2]
# output_dir = basename# paste("results/figures", basename, sep = "/")

data = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$unclassified_cluster == "False",]
#Density plot

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(group=inflation, colour=inflation), alpha=0) + theme_minimal()

output_file = paste(output_dir, 'homogeneity_density_inflation.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()


#COLOR BY FISR NODE
#We don't display the name of first node that appears less than n time
#They are labeled as Other
first_nodes = table(data$first_node)
first_nodes = data.frame(first_nodes, stringsAsFactors = FALSE)
colnames(first_nodes) = c('First_node', 'Freq')
other_nodes = first_nodes$First_node[first_nodes$Freq < 10]
data$first_node[data$first_node %in% other_nodes] = "Other"

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=first_node, group=first_node, fill=first_node), alpha=0.4) + theme_minimal()

output_file = paste(output_dir, 'homogeneity_density_first_node.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()


##BOX PLOT
#Box plots basiques
p <- ggplot(data, aes(x=first_node, y=homogeneity, fill=first_node)) +
  geom_boxplot()  + theme_minimal()
  output_file = paste(output_dir, 'box_plot_homogeneity_firstnode.png', sep ='')

  png(filename=output_file,  width = 1368, height = 768)
  p
  dev.off()

data$includedsDNAViruses = ifelse(grepl('dsDNA viruses', data$valid_branches, fixed=TRUE), TRUE, FALSE)
data$includedsDNAViruses = ifelse(grepl('dsDNA viruses', data$shared_taxonomy, fixed=TRUE), TRUE, FALSE)
data$include_dsDNAViruses_and_Inflation <- paste(data$includedsDNAViruses,data$inflation)
data_no_dsDNA = data[data$includedsDNAViruses == FALSE,]

# data$includedsDNAViruses = ifelse(grepl('dsDNA viruses, no RNA stage', data$first_node, fixed=TRUE), TRUE, FALSE)

p = ggplot(data_no_dsDNA, aes(x=homogeneity)) + geom_density(aes(color=inflation,  group=inflation), alpha=0.4) + theme_minimal()

output_file = paste(output_dir, 'includedsDNAViruses_vs_all_density.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()

data$polyprotein = ifelse(data$nb_polyprotein == 0, FALSE, TRUE)
p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=polyprotein, group=polyprotein, fill=polyprotein), alpha=0.4) + theme_minimal()

output_file = paste(output_dir, 'homogeneity_density_polyprot.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()
##BOX PLOT
#Box plots basiques
p <- ggplot(data, aes(x=polyprotein, y=homogeneity, fill=polyprotein)) +
  geom_boxplot(alpha=0.7)  + theme_minimal()

  output_file = paste(output_dir, 'box_plot_homogeneity_polyprotein.png', sep ='')

  png(filename=output_file,  width = 1368, height = 768)
  p
  dev.off()



  p <- ggplot(data, aes(x=first_node, y=homogeneity, fill=polyprotein)) +
    geom_boxplot(position=position_dodge(1))  + theme_minimal() +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median)


    output_file = paste(output_dir, 'box_plot_homogeneity_firstnode_poly.png', sep ='')
    png(filename=output_file,  width = 1368, height = 768)
    p
    dev.off()

  p <- ggplot(data, aes(x=inflation, y=homogeneity, fill=polyprotein)) +
    geom_boxplot(position=position_dodge(1))  + theme_minimal() +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median)


    output_file = paste(output_dir, 'box_plot_homogeneity_inflation_poly.png', sep ='')
    png(filename=output_file,  width = 1368, height = 768)
    p
    dev.off()

  # p <- ggplot(data, aes(x=inflation, y=homogeneity, fill=includedsDNAViruses)) +
  #   geom_boxplot(position=position_dodge(1))  + theme_minimal() +
  #   stat_summary(fun.data = give.n, geom = "text", fun.y = median)
  #
  #
  #   output_file = paste(output_dir, 'box_plot_homogeneity_inflation_dsDNA.png', sep ='')
  #   png(filename=output_file,  width = 1368, height = 768)
  #   p
  #   dev.off()

  p <- ggplot(data, aes(x=inflation, y=homogeneity, fill=coverage)) +
    geom_boxplot(position=position_dodge(1))  + theme_minimal() +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median)


    output_file = paste(output_dir, 'box_plot_inflation_coverage_all.png', sep ='')
    png(filename=output_file,  width = 1368, height = 768)
    p
    dev.off()

  p <- ggplot(data, aes(x=coverage, y=homogeneity, fill=inflation)) +
    geom_boxplot(position=position_dodge(1))  + theme_minimal() +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median)


    output_file = paste(output_dir, 'box_plot_coverage__inflationall.png', sep ='')
    png(filename=output_file,  width = 1368, height = 768)
    p
    dev.off()

# #Plot le nombre de sequence/genome par cluster
# data = data[data$nb_polyprotein != 0,]
# data_cnt = data.frame( table(data$number_of_genome))
# colnames(data_cnt) = c("number_of_genome", "Number_of_Cluster")
#
# p<-ggplot(data=data_cnt, aes(x=number_of_genome, y=Number_of_Cluster)) +
#   geom_bar(stat="identity") +  theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, size=15, vjust=0.5, hjust=1),
#         axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
#         legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
#   labs( y="Number of genome per clusters", x="number_of_genome")
#
#
# output_file = paste(output_dir, 'Cluster_with_Polyprotein.png', sep ='')
# png(filename=output_file,  width = 800, height = 500)
# p
# dev.off()
#
# p = ggplot(data_cnt, aes(x=number_of_genome)) + geom_density() + theme_minimal()
#
# output_file = paste(output_dir, 'density_nb_of_genome.png', sep ='')
#
# png(filename=output_file,  width = 1368, height = 768)
# p
# dev.off()
