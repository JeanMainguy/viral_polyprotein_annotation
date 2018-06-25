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


#COLOR BY FISR NODE
#We don't display the name of first node that appears less than n time
#They are labeled as Other
first_nodes = table(data$first_node)
first_nodes = data.frame(first_nodes, stringsAsFactors = FALSE)
colnames(first_nodes) = c('First_node', 'Freq')
other_nodes = first_nodes$First_node[first_nodes$Freq < 30]
data$first_node[data$first_node %in% other_nodes] = "Other"

p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=first_node, group=first_node, fill=first_node), alpha=0.4) + theme_minimal()

output_file = paste(output_dir, 'homogeneity_density_first_node.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()


data$polyprotein = ifelse(data$nb_polyprotein == 0, FALSE, TRUE)
p = ggplot(data, aes(x=homogeneity)) + geom_density(aes(colour=polyprotein, group=polyprotein), alpha=0) + theme_minimal()

output_file = paste(output_dir, 'homogeneity_density_polyprot.png', sep ='')

png(filename=output_file,  width = 1368, height = 768)
p
dev.off()

plot(data$number_of_genome, data$homogeneity)
#Plot le nombre de sequence/genome par cluster
data = data[data$nb_polyprotein != 0,]
data_cnt = data.frame( table(data$number_of_genome))
colnames(data_cnt) = c("number_of_genome", "Number_of_Cluster")

p<-ggplot(data=data_cnt, aes(x=number_of_genome, y=Number_of_Cluster)) +
  geom_bar(stat="identity") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number of genome per clusters", x="number_of_genome")


output_file = paste(output_dir, 'Cluster_with_Polyprotein.png', sep ='')
png(filename=output_file,  width = 800, height = 500)
p
dev.off()
