#!/usr/bin/env Rscript
library(ggplot2)
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

print(csv_file)

data = read.csv(file = csv_file, sep = ',', header = TRUE, stringsAsFactors = TRUE)


# inflation='2'
# # data = data[data$type == 'all',]
# print(data$median)
# data = data[data$type != 'dsDNA',]
# data = data[data$type == 'no dsDNA',]

data = data[data$type == 'poly',]
# data = data[data$coverage == 50,]
# data = data[data$evalue >= 1e-40,]
# data = data[data$evalue <= 1e-20,]
# data = data[data$type != 'poly',]
# print(data[data$mean == max(data$mean),])
# print(data[data$median == max(data$median),])
# data = data[data$inflation == inflation,]
# data = data[data$coverage == '40',]
# data = data[data$coverage == '50',]
data$inflation =  paste('I=',data$inflation, sep='')
data$coverage =  paste(data$coverage, "%")
data$coverage_inflation =  paste(data$coverage, data$inflation)
# data$coverage =  paste(data$coverage, data$evalue)
# data$inflation_evalue_coverage =  paste(data$inflation_evalue, data$coverage)
# data$type =  paste(data$type, data$evalue)
data$evalue_s = paste('Evalue', data$evalue)
# p = ggplot(data, aes(x=coverage , fill=inflation,  lower = Q1, upper = Q3,
# middle = median,
# ymin = min,
# ymax = max)) +
# geom_boxplot( stat = "identity") + theme(axis.text.x = element_text(angle = 90, size=5, vjust=0.5, hjust=1),
#       axis.text.y = element_text(size=5), axis.title=element_text(size=15,face="bold"),
#       legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
# labs( y="Number of domain annotations", fill="Interpro Database:", x="Overlap distance (amino acid)")
#
# output_file = paste(output_dir, 'from_summary.png', sep ='')
# png(filename=output_file,  width = 1368, height = 768)
# p
# dev.off()





#
# p = ggplot(data, aes(x=reorder(evalue, evalue), y=coverage_inflation, z=mean)) +
# geom_tile(aes(fill = mean)) + theme_bw() + scale_fill_gradient(low = "red", high = "blue", limit=c(0,1))+
# theme(axis.text.x = element_text(size=13),
#       axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
#       legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
# labs( y="Coverage and inflation", fill="mean", x="Evalue")

#
# output_file = paste(output_dir, 'mean_hedatmap.png', sep ='')
# png(filename=output_file,  width = 600, height = 800)
# p
# dev.off()

p = ggplot(data, aes(x=reorder(evalue, evalue), y=coverage_inflation, z=median)) +
geom_tile(aes(fill = median)) + theme_bw() + scale_fill_gradient(low = "red", high = "blue", limit=c(0,1))+
theme(axis.text.x = element_text(size=13),
      axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"),
      legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
labs( y="Coverage and inflation", fill="median", x="Evalue")

output_file = paste(output_dir, 'median_hedatmap.png', sep ='')
png(filename=output_file,  width = 800, height = 800)
p
dev.off()
