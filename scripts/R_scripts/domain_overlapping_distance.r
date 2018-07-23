library(ggplot2)
library(scales)

csv_file =  "results/stat_viral_protein/stat_domains_Viruses.csv"

data = read.csv(csv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data$len = data$end_in_prot - data$start_in_prot +1 

data = data[data$overlapping == "True",]
data = data[data$protein_fully_covered == "True",]

data$overlapping_len_pc = data$OverlappingDistance / data$len  

overlap_data = data[data$overlapping == 'True', ]
#overlap_data = overlap_data[overlap_data$Expected_peptides_number == 'True', ]


#Group overlap distance
overlap_data$groupedby5 <-cut(overlap_data$OverlappingDistance , seq(0, 150, 10))

#DISTANCE OF OVERLAPING
overlap_database = data.frame( table(overlap_data$groupedby5, overlap_data$method))
colnames(overlap_database) = c("overlap", "method", "Nb_domains")

p<-ggplot(data=overlap_database, aes(x=overlap, y=Nb_domains)) +
  geom_bar(stat="identity",fill = "#f28759") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold")) +
  labs( y="Number of domain annotations", x="Overlap distance (amino acid)")+
  scale_fill_manual(values=c("red"))

png(filename="results/figures/domain_annotations/Overlap_distance_of_domain_annotation",  width = 700, height = 500)
p
dev.off()

#COLOR BY DATABASE
p<-ggplot(data=overlap_database, aes(x=overlap, y=Nb_domains, fill=reorder(method, Nb_domains))) +
  geom_bar(stat="identity") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number of domain annotations", fill="Interpro Database:", x="Overlap distance (amino acid)")

png(filename="results/figures/domain_annotations/overlap_distance_by_database.png",  width = 700, height = 500)
p
dev.off()

#COLOR BY DOMAIN ANNOTATIONS
#We don't display the name of domain that appears less than 5 time
#They are label as Other
domains = table(overlap_data$name)
domains = data.frame(domains, stringsAsFactors = FALSE)
colnames(domains) = c('Domain', 'Freq')
other_domain = domains$Domain[domains$Freq < 10]
overlap_data$name[overlap_data$name %in% other_domain] = "Other"

nb_domains = length(domains$Domain) - length(other_domain)

domain_annotations_data = data.frame( table(overlap_data$groupedby5, overlap_data$name), stringsAsFactors = FALSE)
colnames(domain_annotations_data) = c("overlap", "Domain", "Nb_domains")
domain_annotations_data = rev(domain_annotations_data)

nb_of_annotations = length(domains$Domain) - length(other_domain) -1
domain_annotations_data$sortOther[domain_annotations_data$Domain == "Other"] = 0
domain_annotations_data$sortOther[domain_annotations_data$Domain != "Other"] = 1

p<-ggplot(data=domain_annotations_data, aes(x=overlap, y=Nb_domains, fill=reorder(Domain, sortOther))) +
  geom_bar(stat="identity") +  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number of domain annotations", fill="Domain Name:", x="Overlap distance (amino acid)") 
p = p +  scale_fill_manual(values=c("#999999", hue_pal()(nb_domains)))
p

png(filename="results/figures/domain_annotations/Overlap_distance_by_domain_name",  width = 700, height = 500)
p
dev.off()

#PERCENTAGE OF THE LENGTH THAT OVERLAP


overlap_data$sortOther[overlap_data$name == "Other"] = 0
overlap_data$sortOther[overlap_data$name != "Other"] = 1

p = ggplot(overlap_data, aes(x=overlapping_len_pc)) +
  geom_histogram(aes(y=..density.., fill=reorder(name, sortOther), colour=reorder(name, sortOther)),      # Histogram with density instead of count on y-axis
                 binwidth=0.05) + 
  scale_fill_manual(values=c("#999999", hue_pal()(nb_domains))) + 
  scale_colour_manual(values=c("#999999", hue_pal()(nb_domains)))
  p
  p = ggplot(overlap_data, aes(x=overlapping_len_pc)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=0.01) 
  p
  geom_vline(aes(xintercept=sp_length_threshold),   # Ignore NA values for mean
             color="black", size=.4,  linetype="dashed")+ 
  annotate("text", x=sp_length_threshold+280, y = 0.02, label = signal_threshold_txt, size = 5) +   theme_minimal()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_text(size=22,face="bold"))+
  labs( x="length first peptide annotation (aa)", fill="First peptide\nannotation type", colour="First peptide\nannotation type") +
  scale_x_continuous(breaks=c(sp_length_threshold, seq(200,max(data$potential_signal_P_first_pep_len_aa)+1, 100)))
  p2 = p + geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot + 
  p
  p2

png(filename="results/figures/domain_annotations/Overlap_distance_of_domain_annotation",  width = 700, height = 500)
p
dev.off()



