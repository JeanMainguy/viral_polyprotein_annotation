library(ggplot2)

csv_file = "data/alignment/Viruses/RefSeq_download_date_2018-08-13/Viruses_evalue_1e-60coverage50_I2/stat/cross_validation_stat_alignments.csv"
output_dir = ""
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

data = read.csv(file = csv_file, sep = '\t', header = T, stringsAsFactors = T)
csv_file = "data/alignment/Viruses/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-50coverage60_I1_4/stat/stat_alignments.csv"
data_aln = read.csv(file = csv_file, sep = '\t', header = T, stringsAsFactors = T)



I= data$inflation[1]
E=data$evalue[1]
C=data$coverage[1]
window=data$window[1]
score_threshold=data$confidence_score_threshold[1]

caption_text = paste('\n',length(data[data$nb_seq_with_annotation==1,]),
                     'clusters with only one annotated polyprotein',
                     '\nParameters:\nWindow:', window, "    confidence score threshold:", score_threshold,
                     '\ninflation:',I,
                     '  Evalue max',E, '  Min Coverage:',C)

caption_text = paste('\nParameters:\nWindow:', window, "    confidence score threshold:", score_threshold,
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

#Remove cluster with single annotated cds
data = data[data$nb_seq_with_annotation>1,] 

data$nb_unannotated_seq
p = ggplot(data, aes(x=precision, y=recall, color=aln_validity_display, alpha=0.4)) + 
  geom_point(aes(size=nb_unannotated_seq)) +
  theme_minimal() +
  labs(size="Number of sequences", color="Cluster categories", caption = caption_text, alpha=NULL) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=15), legend.title = element_text(size=18,face="bold"),
        plot.caption = element_text(hjust=0, size=17) )+ 
  scale_alpha_continuous(guide = FALSE)  + 
  scale_color_manual(values=c('forestgreen',   "darkviolet","dodgerblue3",  "firebrick3" ))
p

name=paste("cross_validation_precision_vs_recall_", score_threshold,"_window_", window, sep='')
output_file = paste(output_dir, name,".png", sep ='')
png(filename=output_file,  width = 1000, height = 450)
p
dev.off()

