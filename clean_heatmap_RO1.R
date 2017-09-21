#clean heatmap for RO1
library(pheatmap)
library("RColorBrewer") 


data1 = read.csv(file.choose())
bulk_single= data1[data1[,3] %in% c("KRT18","KRT7","CD36","CD68","CD163"),c(3,10:20,29:48)]
rownames(bulk_single) = bulk_single[,1]
bulk_single = bulk_single[,-1]
bulk = bulk_single[,1:11]
single = data.frame(mean_bulk_DCIS = rowSums(bulk[,1:5])/5, 
                    mean_bulk_MO = rowSums(bulk[,6:11])/6, bulk_single[,12:31])

single2 = single[,c(1:6,10,12,9,8,13,14,11,7,15:22)]
single3 = single2[,c(1:8,13,14,15:22,12)]


myBreaks = unique(c(seq(0, 5.98, length=49), 6, seq(6.01,max(log2(bulk+1)), length=50)))
myColor <- colorRampPalette(c("darkblue", "white", "brown1"))(100)
png("heatmap_bulk.png", width = 1200, height = 700, res=300)
print(pheatmap(log2(bulk+1),breaks = myBreaks, color = myColor, 
               cluster_rows=FALSE, 
               cluster_cols=FALSE,
               gaps_col =  c(5),show_colnames = FALSE))
dev.off()

myBreaks = unique(c(seq(0, 7.98, length=49), 7, seq(7.01,max(log2(single3+1)), length=50)))
myColor <- colorRampPalette(c("darkblue", "white", "brown1"))(100)
png("heatmap_sc.png", width = 2400, height = 800, res=300)
print(pheatmap(log2(single3+1),
               cluster_rows=FALSE, breaks = myBreaks, color = myColor, 
               cluster_cols=FALSE,
               gaps_col =  c(1,2,8,10),show_colnames = FALSE))
dev.off()
