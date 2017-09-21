#############
#heatmaps
#############
library(heatmaply)
library(pheatmap)
which(change1$gene.names == "VGLL1")
EN_IDC_genes = change1$gene.names[c(1:100,1135:1234)]
to_heatmap = change[change$gene.names %in% EN_IDC_genes,c(1,4,6)]
to_heatmap = to_heatmap[order(to_heatmap$IDC_Mono, decreasing = TRUE),]
rownames(to_heatmap) = to_heatmap$gene.names
to_heatmap = to_heatmap[,-1]
heatmaply(log2(to_heatmap+1), dendrogram = "none",showticklabels = c(T,F))
Im_genes_to_hm =unlist(ImGenes)
Im_genes_to_hm = Im_genes_to_hm[!Im_genes_to_hm %in% c("TLR3","MUC4","TGFB2")]
to_heatmap_IMM = change[change$gene.names %in% Im_genes_to_hm, c(1,4,6)]
rownames(to_heatmap_IMM) = to_heatmap_IMM[,1]
to_heatmap_IMM = to_heatmap_IMM[,-1]

heatmaply(log2(to_heatmap[,c(1,2)]+1), dendrogram = "row",showticklabels = c(T,F))
heatmaply(log2(to_heatmap_IMM[,c(1,2)]+1), dendrogram = "row")

print(pheatmap(to_heatmap_IMM,
               cluster_rows=TRUE, 
               cluster_cols=FALSE))

all = cor(al3_mono,use="pairwise.complete.obs",method="pearson")
all = data.frame(all)

library(heatmaply)
heatmaply(all)

apply(Mono_ex,2,function(x){sum(x != 0, na.rm = TRUE)})
complete= Mono_ex[rowSums(!is.na(Mono_ex)) != 0,]
genes = data.frame(rownames(complete))

#remove LN and normal
trimmed = Mono_ex[,-5]
cor_trimmed = cor(trimmed,use="pairwise.complete.obs",method="pearson")
heatmaply(cor_trimmed)
#removeLN from comparison and try comparing DGE between DCIS and IDC cause they have most number of cases
fold = trimmed[,c(1,2,3,4)]
fold[fold <=1] = NA
fold = fold[complete.cases(fold),]
fold = cbind(names = rownames(fold), fold)
fold = fold[,c(1,3,5,4,2)]
library(dplyr)
viral_res = viral_res[,-1]
heatmaply(log2(viral_res[,1:3]), dendrogram = "none", 
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "darkblue", high = "red", midpoint = 8.5),
          margins = c(90,80,NA,0))

trimmed1 = trimmed[complete.cases(trimmed[,2]),]
heatmaply(log2(trimmed1+1))