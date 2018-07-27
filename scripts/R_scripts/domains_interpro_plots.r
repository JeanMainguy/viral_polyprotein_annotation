library(ggplot2)

csv_file =  "results/stat_viral_protein/domain_stat_Viruses.csv"

data = read.csv(csv_file, header = TRUE, sep = "\t")
data = data[data$polyprotein_outline == "True",]
data = data[data$protein_fully_covered == "True",]
#SIMPLE HISTOGRAMME WITH OVERLAPPING AND NON OVERLAPPING
table = data.frame( table(data$overlapping))
colnames(table) = c("Fully_included", "Nb_domains")
table$Overlap = ifelse(table$Fully_included == 'False', 'Non overlapping', 'Overlapping')

p<-ggplot(data=table, aes(x=Overlap, y=Nb_domains, fill=Overlap)) +geom_bar(stat="identity") + theme_minimal() +
  labs(x =NULL, y="Number of domain annotations", fill=NULL) + theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=18) )

p = p+scale_fill_manual(values=c( "#56B4E9", "firebrick"), guide=FALSE)
p = p+  geom_text(aes(label=Nb_domains), vjust=-0.1, color="black",
                  position = position_dodge(0.95), size=7)

p

png(filename="results/figures/domain_annotations/Overlaping_and_NonOverlaping_Domains.png",  width = 600, height = 600)
p
dev.off()



