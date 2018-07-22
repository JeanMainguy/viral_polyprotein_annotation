#!/usr/bin/env Rscript
library(ggplot2)

csv_file="results/clustering_evaluation/homogeneity_summary_stat.csv"
coverage = 20
inflation="2.0"
evalue=1e-30

output_dir = "results/clustering_evaluation/"

data = read.csv(file = csv_file, sep = ',', header = TRUE, stringsAsFactors = TRUE)
data = data[data$type == 'poly',]
data$inflation = ifelse(grepl('_', data$inflation), gsub("_", ".", data$inflation), paste(data$inflation, '.0', sep=''))

#data$inflation = gsub("_", ".", data$inflation)

# COVERAGE FIXED = 20

dataC = data[data$coverage == coverage,]
caption_text = paste('\nCoverage threshold:',coverage, '%' )

#Filter some of  to display less dataC
dataC = dataC[dataC$inflation != '1.8',]
#dataC = dataC[dataC$inflation == '2',]
#dataC = dataC[dataC$coverage == 20,]
dataC = dataC[dataC$evalue != 1e-100,]
dataC = dataC[dataC$evalue != 1e-120,]
dataC = dataC[dataC$evalue != 1e-70,]
dataC = dataC[dataC$evalue != 1e-40,]
dataC = dataC[dataC$evalue != 1e-20,]
dataC = dataC[dataC$evalue != 1e-60,]
#dataC = dataC[dataC$evalue != 1e-10,]

#dataC$inflation =  paste('I=',dataC$inflation, sep='')
dataC$coverage =  paste(dataC$coverage, "%")
# dataC$coverage =  paste(dataC$coverage, dataC$evalue)
# dataC$inflation_evalue_coverage =  paste(dataC$inflation_evalue, dataC$coverage)
# dataC$type =  paste(dataC$type, dataC$evalue)

dataC$evalue_s = paste('Evalue', dataC$evalue)
# scale_fill_gradient(low = "red", high = "blue", limit=c(0,1))+

p = ggplot(dataC, aes(x=reorder(evalue, evalue), y=inflation, z=median)) +
  geom_tile(aes(fill = median))+ 
  scale_fill_gradient2(low = "yellow", high = "blue",mid='green',  midpoint = 0.6, limit=c(min(data$median),max(data$median)))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_text(size=22,face="bold"),
        plot.caption = element_text(hjust=0, size=22))+
  labs( y="Inflation", fill="Homogeneity\nmedian", x="Evalue", caption = caption_text)

p
name=paste("homogeneity_median_heatmap_coverage_of_", coverage, sep='')
output_file = paste(output_dir, name,".png", sep ='')
png(filename=output_file,  width = 950, height = 950)
p
dev.off()


# Inflation FIXED = 2

dataC = data[data$inflation == inflation,]
caption_text = paste('\nInflation:',inflation)

#Filter some of  to display less dataC
dataC = dataC[dataC$inflation != '1.8',]
#dataC = dataC[dataC$inflation == '2',]
#dataC = dataC[dataC$coverage == 20,]
dataC = dataC[dataC$evalue != 1e-100,]
dataC = dataC[dataC$evalue != 1e-120,]
dataC = dataC[dataC$evalue != 1e-70,]
dataC = dataC[dataC$evalue != 1e-40,]
dataC = dataC[dataC$evalue != 1e-20,]
dataC = dataC[dataC$evalue != 1e-60,]
#dataC = dataC[dataC$evalue != 1e-10,]

dataC$inflation =  paste('I=',dataC$inflation, sep='')
dataC$coverage =  paste(dataC$coverage, "%")
# dataC$coverage =  paste(dataC$coverage, dataC$evalue)
# dataC$inflation_evalue_coverage =  paste(dataC$inflation_evalue, dataC$coverage)
# dataC$type =  paste(dataC$type, dataC$evalue)

dataC$evalue_s = paste('Evalue', dataC$evalue)
# scale_fill_gradient(low = "red", high = "blue", limit=c(0,1))+

p = ggplot(dataC, aes(x=reorder(evalue, evalue), y=coverage, z=median)) +
  geom_tile(aes(fill = median))+ 
  scale_fill_gradient2(low = "yellow", high = "blue",mid='green',  midpoint = 0.6, limit=c(min(data$median),max(data$median)))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_text(size=22,face="bold"),
        plot.caption = element_text(hjust=0, size=22))+
  labs( y="Coverage", fill="Homogeneity\nmedian", x="Evalue", caption = caption_text)

p
name=paste("homogeneity_median_heatmap_inflation_of_", inflation, sep='')
output_file = paste(output_dir, name,".png", sep ='')
png(filename=output_file,  width = 950, height = 950)
p
dev.off()

# Evalue FIXED

dataC = data[data$evalue == evalue,]
caption_text = paste('\nEvalue:',evalue)

#Filter some of  to display less dataC
dataC = dataC[dataC$inflation != '1.8',]
#dataC = dataC[dataC$inflation == '2',]
#dataC = dataC[dataC$coverage == 20,]
dataC = dataC[dataC$evalue != 1e-100,]
dataC = dataC[dataC$evalue != 1e-120,]
dataC = dataC[dataC$evalue != 1e-70,]
dataC = dataC[dataC$evalue != 1e-40,]
dataC = dataC[dataC$evalue != 1e-20,]
dataC = dataC[dataC$evalue != 1e-60,]
#dataC = dataC[dataC$evalue != 1e-10,]

dataC$inflation =  paste('I=',dataC$inflation, sep='')
dataC$coverage =  paste(dataC$coverage, "%")
# dataC$coverage =  paste(dataC$coverage, dataC$evalue)
# dataC$inflation_evalue_coverage =  paste(dataC$inflation_evalue, dataC$coverage)
# dataC$type =  paste(dataC$type, dataC$evalue)

dataC$evalue_s = paste('Evalue', dataC$evalue)
# scale_fill_gradient(low = "red", high = "blue", limit=c(0,1))+

p = ggplot(dataC, aes(x=inflation, y=coverage, z=median)) +
  geom_tile(aes(fill = median))+ 
  scale_fill_gradient2(low = "yellow", high = "blue",mid='green',  midpoint = 0.6, limit=c(min(data$median),max(data$median)))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_text(size=22,face="bold"),
        plot.caption = element_text(hjust=0, size=22))+
  labs( y="Coverage", fill="Homogeneity\nmedian", x="Inflation", caption = caption_text)

p
name=paste("homogeneity_median_heatmap_evalue_of_", evalue, sep='')
output_file = paste(output_dir, name,".png", sep ='')
png(filename=output_file,  width = 950, height = 950)
p
dev.off()

#Heatmap global 

data_all = data
#Filter some of  to display less data
data_all = data_all[data_all$inflation != '1.8',]
#data_all = data_all[data_all$inflation == '2',]
#data_all = data_all[data_all$coverage == 20,]
data_all = data_all[data_all$evalue != 1e-100,]
data_all = data_all[data_all$evalue != 1e-120,]
data_all = data_all[data_all$evalue != 1e-70,]
data_all = data_all[data_all$evalue != 1e-40,]
data_all = data_all[data_all$evalue != 1e-20,]
data_all = data_all[data_all$evalue != 1e-60,]
data_all = data_all[data_all$evalue != 1e-10,]
data_all = data_all[data_all$coverage != 10,]
#data_all = data_all[data_all$evalue != 1e-10,]
data_all$coverage =paste(data_all$coverage, '%', sep='')

data_all$coverage_inflation = paste(data_all$coverage, 'I=',data_all$inflation)

p = ggplot(data_all, aes(x=reorder(evalue, evalue), y=coverage_inflation, z=median)) +
  geom_tile(aes(fill = median))+ 
  scale_fill_gradient2(low = "yellow", high = "blue",mid='green',  midpoint = 0.6, limit=c(min(data$median),max(data$median)))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18), axis.title=element_text(size=22,face="bold"),
        legend.text = element_text( size=20), legend.title =element_text(size=22,face="bold"),
        plot.caption = element_text(hjust=0, size=22))+
  labs( y="Coverage and Inflation", fill="Homogeneity\nmedian", x="Evalue")

p

output_file = paste(output_dir, "homogeneity_median_heatmap_ALL.png", sep ='')
png(filename=output_file,  width = 900, height = 1000)
p
dev.off()
