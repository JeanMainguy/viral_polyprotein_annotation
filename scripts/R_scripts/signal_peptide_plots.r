#!/usr/bin/env Rscript
library(ggplot2)

csv_file="results/stat_viral_protein/stat_proteins_Viruses.csv"
output_dir = "results/figures/"
sp_length_threshold = 90
signal_threshold_txt = paste("Signal peptide\nlength threshold =", sp_length_threshold)
len_txt = nchar(signal_threshold_txt)

data = read.csv(file = csv_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data = data[data$potential_signal_P_first_pep_len != 'None', ]
data$potential_signal_P_first_pep_len = sapply(data$potential_signal_P_first_pep_len, as.numeric) 

data$potential_signal_P_first_pep_len_aa = data$potential_signal_P_first_pep_len / 3 

p = ggplot(data, aes(x=potential_signal_P_first_pep_len_aa)) +
  geom_histogram(aes(y=..density.., fill=potential_signal_P_first_pep_type, colour=potential_signal_P_first_pep_type),      # Histogram with density instead of count on y-axis
                 binwidth=10) +
  
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
  
  name="first_peptide_distribution"
  output_file = paste(output_dir, name,".png", sep ='')
  png(filename=output_file,  width = 1400, height = 600)
  p
  dev.off()
  
