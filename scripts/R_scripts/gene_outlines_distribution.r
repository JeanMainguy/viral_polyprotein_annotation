#!/usr/bin/env Rscript
library(ggplot2)

csv_file="results/stat_viral_protein/stat_proteins_Viruses.csv"
output_dir = "results/figures/"

data_original = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

data = data_original[data_original$has_peptide == "True" ,]

data$status = ifelse(data$polyprotein_outline == "True", 'Polyprotein', data$non_polyprotein_explanation)

data$status = ifelse(data$non_polyprotein_explanation == "Intein outline: extein includes intein", 
                     'Intein-extein\nannotation', 
                     data$status)
data$status = ifelse(data$non_polyprotein_explanation == "Intein outline: extein surounds intein", 
                     'Intein-extein\nannotation', 
                     data$status)

data$status = ifelse(data$non_polyprotein_explanation == "single peptide annotaion covering the whole CDS", 
                     'single peptide annotaion covering almost the whole CDS', 
                     data$status)

data$status = ifelse(data$status == "single peptide annotaion covering almost the whole CDS", 
                     'Single peptide\nannotation covering\nalmost the whole CDS', 
                     data$status)

#SIMPLE HISTOGRAMME WITH outline type
table = data.frame( table(data$status))
colnames(table) = c("outline", "Nb_protein")

p <- ggplot(data=table, aes(x=reorder(outline, Nb_protein), y=Nb_protein, fill=outline)) +geom_bar(stat="identity") + theme_minimal() +
  labs(x =NULL, y="Number of CDS", fill=NULL) + 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"))  + guides(fill=FALSE)

#p = p+scale_fill_manual(values=c("#E69F00", "#56B4E9"), guide=FALSE)
p = p+  geom_text(aes(label=Nb_protein), vjust=-0.1, color="black",
                  position = position_dodge(0.95), size=7)
p

name="annotated_protein_outlines_distribution"
output_file = paste(output_dir, name,".png", sep ='')
png(filename=output_file,  width = 1000, height = 600)
p
dev.off()

