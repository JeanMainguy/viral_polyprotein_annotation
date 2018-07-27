#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
if (FALSE){
  csv_file='data/alignment/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-30coverage20_I2/stat_of_all_cluster.csv'
  basename = "test/"
  window = 20
}
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  else if (length(args)==1) {
  # default output file
  args[2] = "plot"
  args[3] = 10
}
csv_file = args[1]
basename = args[2]
window = args[3]
output_dir = basename # paste("results/figures", basename, sep = "/")
output_dir = paste(output_dir, window, sep='')
print(csv_file)

data = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = TRUE)
data = data[data$window == window,] 
data_cnt = data.frame( table(data$nb_of_cleavage_site, data$completness_category, data$round_std))
colnames(data_cnt) = c("nb_of_cleavage_site_per_groupe", "Completness", "standard_deviation","Nb_of_Cleavage_Site_Group")

p<-ggplot(data=data_cnt, aes(x=standard_deviation, y=Nb_of_Cleavage_Site_Group, fill=Completness)) +
  geom_bar(stat="identity") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number of cleavage site groups", fill="Completness:", x="standard_deviation")


output_file = paste(output_dir, 'standard_deviation_distribution.png', sep ='')
png(filename=output_file,  width = 900, height = 500)
p
dev.off()

p<-ggplot(data=data_cnt, aes(x=nb_of_cleavage_site_per_groupe, y=Nb_of_Cleavage_Site_Group, fill=reorder(Completness, Nb_of_Cleavage_Site_Group))) +
  geom_bar(stat="identity") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size=12, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number Cleavage Site groups", fill="Completness:", x="number of cleavage site per group")


  output_file = paste(output_dir, 'cleavage_site_groupe_distribution_completeness.png', sep ='')
  png(filename=output_file,  width = 1000, height = 500)
  p
  dev.off()


p<-ggplot(data=data_cnt, aes(x=nb_of_cleavage_site_per_groupe, y=Nb_of_Cleavage_Site_Group, fill=standard_deviation)) +
  geom_bar(stat="identity") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size=12, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number Cleavage Site groups", fill="Standard deviation:", x="number of cleavage site per group")

  output_file = paste(output_dir, 'cleavage_site_groupe_distribution_standard_deviation.png', sep ='')

  png(filename=output_file,  width = 1000, height = 500)
  p
  dev.off()
