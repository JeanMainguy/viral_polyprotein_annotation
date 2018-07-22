#!/usr/bin/env Rscript
library(ggplot2)

csv_file="results/stat_viral_protein/stat_peptides_domains_Viruses.csv"
output_dir = "results/figures/"
name= 'barplot_peptide_includs_domains.png'

data_o = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = TRUE)
data = data_o[data_o$polyprotein_outline == 'True', ]
#data = data_o

#SIMPLE HISTOGRAMME WITH OVERLAPPING AND NON OVERLAPPING

table_status = data.frame( table(data$cover_by_domain, data$overlapped_by_domain))

colnames(table_status) = c("includs_full_domain", "includs_partially_domain", 'Nb_of_Peptides')


table_status$full_partial = paste(table_status$includs_full_domain,table_status$includs_partially_domain, sep='_' )

table_status$status = ifelse(table_status$full_partial == "True_False", 
                             "only full domains included ", 
                             "NA")

table_status$status = ifelse(table_status$full_partial == "True_True", 
                             "full and partial domains included", 
                             table_status$status)

table_status$status = ifelse(table_status$full_partial == "False_False", 
                             "includs no domains", 
                             table_status$status)

table_status$status = ifelse(table_status$full_partial == "False_True", 
                             "only partial domains included", 
                             table_status$status)



table_status$includs_domains = ifelse((table_status$includs_full_domain == "True"), 
                                      "Includs domains",
                                      'Includs no domains')
table_status$includs_domains = ifelse((table_status$includs_partially_domain == "True"), 
                                      "Includs domains",
                                      table_status$includs_domains)


p<-ggplot(data=table_status, aes(x=includs_domains, y=Nb_of_Peptides, fill=reorder(status, Nb_of_Peptides ))) + 
  geom_bar(stat="identity") + 
  theme_minimal() +
  labs(x =NULL, y="Number of domain annotations", fill=NULL) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=18) )

p = p+scale_fill_manual(values=c("goldenrod","darkorange3", "azure4", "firebrick4"))
p


output = paste(output_dir, name, sep='')

png(filename=output,  width = 2000, height = 900, res=300)
p
dev.off() 


data$coverage = (data$length_covered_by_all_domains/(data$len)) * 100 
summary(data$coverage)
