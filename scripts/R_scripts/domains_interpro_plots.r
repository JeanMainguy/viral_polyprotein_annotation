library(ggplot2)

csv_file =  "results/stat_viral_protein/stat_domains_Viruses.csv"

data = read.csv(csv_file, header = TRUE, sep = "\t")

#SIMPLE HISTOGRAMME WITH OVERLAPPING AND NON OVERLAPPING
table = data.frame( table(data$fully.included))
colnames(table) = c("Fully_included", "Nb_domains")
table$Overlap = ifelse(table$Fully_included == 'True', 'Non overlapping', 'Overlapping')

p<-ggplot(data=table, aes(x=Overlap, y=Nb_domains, fill=Overlap)) +geom_bar(stat="identity") + theme_minimal() +
  labs(x =NULL, y="Number of domain annotations", fill=NULL) + theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=18) )

p = p+scale_fill_manual(values=c("#E69F00", "#56B4E9"), guide=FALSE)
p = p+  geom_text(aes(label=Nb_domains), vjust=-0.1, color="black",
                  position = position_dodge(0.95), size=7)

png(filename="results/figures/Overlaping_and_NonOverlaping_Domains.png",  width = 400, height = 500)
p
dev.off()



