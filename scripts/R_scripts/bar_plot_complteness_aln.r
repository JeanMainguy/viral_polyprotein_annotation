#!/usr/bin/env Rscript
library(ggplot2)

### 

csv_file="data/alignment/Viruses/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-50coverage40_I1_4/stat_of_all_cluster_summary.csv"
csv_file='data/alignment/Viruses/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-60coverage40_I2/stat_of_all_cluster_summary.csv'

I= '2'
E='1e-60'
C='40'

caption_text = paste('\nParamters:\ninflation:',I, '  Evalue max',E, '  Coverage:',C, '\nWindow:', 30)

data = read.csv(file = csv_file, sep = '\t', header = T, stringsAsFactors = T)

#Make distinction between very good cluster that have all cleavages site group have 100% completeness.  

unannotated_valid="unannotated sequence\nin valid cluster"
unannotated_invalid="unannotated sequence\nin invalid cluster"
unannotated_single = "unannotated sequence in\ncluster with only one\nannotated polyprotein"
unannotated_perfect = 'unannotated sequence\nin perfect cluster'

data$aln_validity_prot[data$aln_validity == "majority of cleavage site groups agree"] = unannotated_valid
data$aln_validity_prot[data$aln_validity == "cleavage site groups are divergent"] = unannotated_invalid
data$aln_validity_prot[data$aln_validity == "single annotated polyprotein"] = unannotated_single
data$aln_validity_prot[data$aln_validity == "All cleavage sites groups are complete"] = unannotated_perfect




p<-ggplot(data=data, aes(x=aln_validity_prot, y=nb_unannotated_seq, fill=aln_validity_prot)) + 
  geom_bar(stat="identity") + 
  theme_minimal() +
  guides(fill=FALSE)+
  labs(x =NULL, y="Number of unnannotated proteins", fill=NULL,caption = caption_text) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=18),  
        plot.caption = element_text(hjust=0, size=22) )+
  scale_x_discrete(limits=c(unannotated_invalid, unannotated_perfect, unannotated_valid ,unannotated_single ))
p = p+scale_fill_manual(values=c( "grey60",  "coral", 'darkolivegreen3',  "cadetblue1")) #green-yellow yellow grey green 
p

unannotated_valid= "valid cluster\nmajority of\ncleavage site groups agree"
unannotated_invalid="invalid cluster\ncleavage site\ngroups are divergent"
unannotated_single = "cluster with only one\nannotated polyprotein"
unannotated_perfect = "perfect cluster\nAll cleavage sites\ngroups agree"

data$aln_validity_display[data$aln_validity == "majority of cleavage site groups agree"] = unannotated_valid
data$aln_validity_display[data$aln_validity == "cleavage site groups are divergent"] =unannotated_invalid
data$aln_validity_display[data$aln_validity == "single annotated polyprotein"] = unannotated_single 
data$aln_validity_display[data$aln_validity == "All cleavage sites groups are complete"] = unannotated_perfect
table_status = data.frame( table(data$aln_validity_display))
colnames(table_status) = c("cluster_validity", 'number_of_cluster')

p<-ggplot(data=table_status, aes(x=cluster_validity, y=number_of_cluster, fill=cluster_validity)) + 
  geom_bar(stat="identity") + 
  theme_minimal() +
  guides(fill=FALSE)+
  labs(x =NULL, y="Number of clusters", fill=NULL, caption = caption_text) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=18),
        plot.caption = element_text(hjust=0, size=22) )+
  scale_x_discrete(limits=c(unannotated_invalid, unannotated_perfect, unannotated_valid ,unannotated_single ))
p = p+scale_fill_manual(values=c( "dimgrey",  "firebrick3", 'forestgreen',  "dodgerblue3")) #green-yellow yellow grey green 
p





