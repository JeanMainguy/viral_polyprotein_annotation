library(ggplot2)

csv_file = "data/alignment/Viruses/RefSeq_download_date_2018-08-13/Viruses_evalue_1e-20coverage60_I1_4/stat/stat_alignments.csv"
csv_file = "data/alignment/Viruses/RefSeq_download_date_2018-08-13/Viruses_evalue_1e-20coverage60_I1_8/stat/stat_alignments.csv"

#output_dir = "data/alignment/Viruses/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-40coverage70_I2/stat/"
args = commandArgs(trailingOnly=TRUE)
output_dir = 'test/'

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  else if (length(args)==1) {
  # default output file base
  args[2] = "test/test_plot"
}
csv_file = args[1]
output_dir = args[2]



data = read.csv(file = csv_file, sep = '\t', header = T, stringsAsFactors = T)

I= data$inflation[1]
E=data$evalue[1]
C=data$coverage[1]
window=data$window[1]
score_threshold=data$confidence_score_threshold[1]

score_threshold=data$confidence_score_threshold[1]

caption_text = paste('\nParameters:\nWindow:', window, "    Confidence value treshold:", score_threshold,
                     '\ninflation:',I,
                     '  Evalue max',E, '  Min Coverage:',C)

#Make distinction between very good cluster that have all cleavages site group have 100% completeness.

valid_cluster= "Majority of cleavage site groups are valid"
invalid_cluster="None of the cleavage site groups are valid"
single_cluster = "Cluster with only one annotated polyprotein"
okish_cluster = "At least one cleavage site group is valid"
perfect_cluster = "All cleavage site groups are valid"

valid_cluster_global = 'Valid cluster'
invalid_cluster_global = 'Invalid cluster'
single = 'cluster with only one\nannotated polyprotein'


data$aln_validity_display[data$aln_validity == "Number of valid cleavage_site groups >= minimum number of cleavage site/cds annotated"] = valid_cluster
data$aln_validity_display[data$aln_validity == "None of the cleavage site groups are valid"] = invalid_cluster
data$aln_validity_display[data$aln_validity == "single annotated polyprotein"] = single_cluster
data$aln_validity_display[data$aln_validity == "All cleavage sites groups have a good score"] = perfect_cluster
data$aln_validity_display[data$aln_validity == "At least one cleavage site group is valid"] = okish_cluster

data$global_category[data$aln_validity == "Number of valid cleavage_site groups >= minimum number of cleavage site/cds annotated"] = valid_cluster_global
data$global_category[data$aln_validity == "All cleavage sites groups have a good score"] = valid_cluster_global
data$global_category[data$aln_validity == "At least one cleavage site group is valid"] = valid_cluster_global

data$global_category[data$aln_validity == "single annotated polyprotein"] = single
data$global_category[data$aln_validity == "None of the cleavage site groups are valid"] = invalid_cluster_global



table_status = data.frame( table(data$global_category,data$aln_validity_display ))
colnames(table_status) = c("cluster_validity_global","cluster_validity", 'number_of_cluster')
order_category = c(perfect_cluster, valid_cluster ,okish_cluster, single_cluster, invalid_cluster)
table_status$cluster_validity = factor(table_status$cluster_validity, levels = order_category)

p<-ggplot(data=table_status, aes(x=cluster_validity_global, y=number_of_cluster, fill=cluster_validity)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(x =NULL, y="Number of clusters", fill="Cluster categories", caption = caption_text) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=15), legend.title = element_text(size=18,face="bold"),
        plot.caption = element_text(hjust=0, size=15) )+
  scale_x_discrete(limits=c(valid_cluster_global,  single, invalid_cluster_global))
p = p+scale_fill_manual(values=c('forestgreen',   "dodgerblue3", "darkviolet", "dimgrey",  "firebrick3" )) #green-yellow yellow grey green
#c(perfect_cluster, valid_cluster ,okish_cluster, single_cluster, invalid_cluster))
p

name=paste("distribution_cluster_category_treshold_score_", score_threshold, sep='')
output_file = paste(output_dir, name,".png", sep ='')
png(filename=output_file,  width = 1000, height = 450)
p
dev.off()

valid_cluster_global = 'Unannotated cds\nfrom valid cluster'
invalid_cluster_global = 'Unannotated cds\nfrom invalid cluster'
single = 'Unannotated cds\nfrom cluster with\nonly one\nannotated polyprotein'



data$global_category_prot = valid_cluster
data$global_category_prot[data$aln_validity == "Number of valid cleavage_site groups >= minimum number of cleavage site/cds annotated"] = valid_cluster_global
data$global_category_prot[data$aln_validity == "All cleavage sites groups have a good score"] = valid_cluster_global
data$global_category_prot[data$aln_validity == "At least one cleavage site group is valid"] = valid_cluster_global

data$global_category_prot[data$aln_validity == "single annotated polyprotein"] = single
data$global_category_prot[data$aln_validity == "None of the cleavage site groups are valid"] = invalid_cluster_global


table_status = data.frame( table(data$global_category_prot,data$aln_validity_display))
colnames(table_status) = c("cluster_validity_global","cluster_validity", 'number_of_cluster')
order_category = c(perfect_cluster, valid_cluster ,okish_cluster, single_cluster, invalid_cluster)
table_status$cluster_validity = factor(table_status$cluster_validity, levels = order_category)

table_status$nb_unannotated_seq = 0
table_status$nb_unannotated_seq[table_status$cluster_validity ==  valid_cluster & table_status$number_of_cluster >  0] = sum(data$nb_unannotated_seq[data$aln_validity_display == valid_cluster] )
table_status$nb_unannotated_seq[table_status$cluster_validity ==  invalid_cluster & table_status$number_of_cluster >  0] = sum(data$nb_unannotated_seq[data$aln_validity_display == invalid_cluster] )
table_status$nb_unannotated_seq[table_status$cluster_validity ==  okish_cluster & table_status$number_of_cluster >  0] = sum(data$nb_unannotated_seq[data$aln_validity_display == okish_cluster] )
table_status$nb_unannotated_seq[table_status$cluster_validity ==  single_cluster & table_status$number_of_cluster >  0] = sum(data$nb_unannotated_seq[data$aln_validity_display == single_cluster] )
table_status$nb_unannotated_seq[table_status$cluster_validity ==  perfect_cluster & table_status$number_of_cluster >  0] = sum(data$nb_unannotated_seq[data$aln_validity_display == perfect_cluster] )

p<-ggplot(data=table_status, aes(x=cluster_validity_global, y=nb_unannotated_seq, fill=cluster_validity)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(x =NULL, y="Number of unnannotated cds", fill="Cluster categories", caption = caption_text) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=15), legend.title = element_text(size=18,face="bold"),
        plot.caption = element_text(hjust=0, size=15) )+
  scale_x_discrete(limits=c(valid_cluster_global,  single, invalid_cluster_global))
p = p+scale_fill_manual(values=c('forestgreen',   "dodgerblue3", "darkviolet", "dimgrey",  "firebrick3" )) #green-yellow yellow grey green
#c(perfect_cluster, valid_cluster ,okish_cluster, single_cluster, invalid_cluster))
p
name=paste("unnannotated_CDS_distribution_cluster_category_treshold_score_", score_threshold, sep='')
output_file = paste(output_dir, name,".png", sep ='')
png(filename=output_file,  width = 1000, height = 450)
p
dev.off()
