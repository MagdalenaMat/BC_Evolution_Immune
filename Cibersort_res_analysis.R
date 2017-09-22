mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}


#############
path = "/home/magda/Desktop/cibersortX/run2/CV_filtered"
setwd(path)
files <- list.files(path=path)

#extract names
namesGEP = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_bulk_","_unames.txt_GEPs_CVFilteredAll.txt"),c("",""),i)
  namesGEP = c(namesGEP, nam)
}

#put all files in the list
filteredGEP = list()
for (i in 1: length(files)){
  r = read.table(files[i], sep = "\t", header = TRUE)
  filteredGEP[[i]] = r
}
names(filteredGEP) = namesGEP

gene_count = data.frame()
for(i in namesGEP){
  gene_count = rbind(gene_count, apply(filteredGEP[[i]], 2, function(x){sum(!is.na(x))}))
}

rownames(gene_count) = namesGEP
colnames(gene_count) =colnames(filteredGEP[[1]])

#extract all gene names from the files and put them in the list to find intersect 
#of genes acros all the histological subtypes
results = list()
for (i in namesGEP){
  results[[i]] = as.character(filteredGEP[[i]][,1])
}

universe = Reduce(intersect, results)

#extract genes in monocyte signature common for all samples 
Mono_ex = data.frame(matrix(NA, nrow = 9748, ncol = 1))
for (i in namesGEP){
  Mono = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,6)] #select gene_names and Monocytes
  Mono_ex = cbind(Mono_ex,Mono)
}

all(as.character(Mono_ex[,10])==as.character(Mono_ex[,4]))
rownames(Mono_ex) = as.character(Mono_ex[,2])
Mono_ex = Mono_ex[, seq(3, length(Mono_ex), 2)]
colnames(Mono_ex) = c("AVL","DCIS","EN","IDC","LN","MetECE","normal")

Mono_ex[Mono_ex<=1] <- NA
predicted_genes_Monocytes = apply(Mono_ex, 2, function(x){sum(!is.na(x))})

#write.table(predicted_genes_Monocytes, "predicted_genes_Monocytes_in_stages.txt", col.names = FALSE, quote = FALSE)

#Mono_ex_DCIS = Mono_ex[!is.na(Mono_ex$DCIS_Mono),] 
#cor(Mono_ex,use="pairwise.complete.obs",method="pearson") # can't use spearman here
al2_mono = Mono_ex[rowSums(!is.na(Mono_ex)) >= 2,]
al3_mono = Mono_ex[rowSums(!is.na(Mono_ex)) >= 3,]
# library(heatmaply)
# heatmaply_na(al3_mono, showticklabels = c(T,F))

##########
#plot immune related gene values in different stages 
###########
library(plotly)
library(reshape2)
library(dplyr)
library(ggplot2)
to_reshape = data.frame(gene.names = rownames(al2_mono), al2_mono)
to_reshape = to_reshape[,c(1,8,6,4,3,5,2,7)]

#compute log2 Fold Changes
change = mutate(to_reshape, EN_DCIS = log2(to_reshape$DCIS/to_reshape$EN), 
                DCIS_IDC = log2(to_reshape$IDC/to_reshape$DCIS),
                EN_IDC = log2(to_reshape$IDC/to_reshape$EN))

ImGenes = list()
ImGenes[["HLA"]] = c("HLA.DMA","HLA.DPB1","HLA.DRB1","HLA.DQB1")
ImGenes[["TNF"]] = c("TNFSF12","TNFRSF1A","TRADD")
ImGenes[["NFkB"]] = c("NFKBIA","MYD88")
#ImGenes[["B2M"]] = c("B2M")
ImGenes[["CD"]] = c("CD14","CD163")
ImGenes[["IFN"]] = c("IFNGR1","IFNGR2","NFKBIA")
ImGenes[["selected"]] = c("CD14", "CD163","NFKBIA")
genes = c("CD","HLA","IFN","TNF","NFkB", "selected")  
genes = c("CD", "HLA","IFN", "selected")


for(i in genes){
  Immuno = change[change$gene.names %in% ImGenes[[i]],]
  Immuno = Immuno[,c(1,2,4,5,6)]
  Immuno$gene.names = as.factor(as.character(Immuno$gene.names))
  print(Immuno)
  long_Imm = melt(Immuno, id.vars = "gene.names")
  colnames(long_Imm) = c("gene.names","stage","gene.count")
  print(ggplot(long_Imm,aes(gene.names,gene.count, group = stage, fill =stage)) + 
          geom_col(position = "dodge")+ theme_classic()+
          theme(text = element_text(size=20))+
          scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red")))
}

png('selectedv2.png',width = 1200, height = 1200, res = 120)
print(ggplot(long_Imm,aes(gene.names,gene.count, group = stage, fill =stage)) + 
        geom_col(position = "dodge")+ theme_classic()+
        theme(text = element_text(size=25))+
        scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red")))
dev.off()

