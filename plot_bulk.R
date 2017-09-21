normal = read.table("/home/magda/Desktop/cibersortX/bulk_normal_unames.txt", sep = "\t", header = TRUE)
EN = read.table("/home/magda/Desktop/cibersortX/bulk_EN_unames.txt", sep = "\t", header = TRUE)
DCIS = read.table("/home/magda/Desktop/cibersortX/bulk_DCIS_unames.txt", sep = "\t", header = TRUE)
IDC = read.table("/home/magda/Desktop/cibersortX/bulk_IDC_unames.txt", sep = "\t", header = TRUE)

Im_genes = c("CTSB","MMP9","CD14", "CD163","HLA-DMA","HLA-DPB1","HLA-DRB1","HLA-DQB1","NFKBIA",
             "TNFSF12","TNFRSF1A","TRADD","CD80","APOE", "CSF1R","CD47","TGFB2","MYD88","TLR3","IFNAR1","IFNAR2","IFNGR1","IFNGR2","STAT1")

Im_genes = c("CD14", "CD163","NFKBIA")
stages = list(normal,EN,DCIS,IDC)
names(stages) = c("normal","EN","DCIS","IDC")

for(i in 1:length(stages)){
  rownames(stages[[i]])= stages[[i]][,1]
  stages[[i]]= stages[[i]][rownames(stages[[i]]) %in% Im_genes,-1]
}

gene_names = rownames(stages[[1]])

# res = data.frame(gene=c(),stage=c(),value=c())
# for(i in 1:length(gene_names)){
#   for(j in 1:length(stages)){
#     new.row = data.frame(gene=gene_names[i],
#                          stage=names(stages)[j],
#                          value=t(stages[[j]])[,i])
#     res = rbind(res, new.row)
#   }
# }

res2 = data.frame()
for(i in 1:length(stages)){
  nstag = data.frame(t(stages[[i]]), stage = names(stages)[i])
  res2 = rbind(res2,nstag)
}


library(ggplot2)
library(ggpubr)
genes = colnames(res2)[1:3]
for(gene in genes){
  pvalues = t.test(res2[[gene]] ~ res2$stage, res2)$p.value
  formatedp = format.pval(pvalues,digits = 2)
  lbl = paste0("Welch's t-test\n p = ",formatedp,sep=" ")
  print(ggboxplot(res2, x = "stage", y= gene, color = "stage", palette = "jco",
                  add = "jitter", outlier.colour = NA)+
          theme(panel.grid.major = element_blank(), legend.position = "none", axis.text=element_text(size=12),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                plot.title = element_text(hjust = 0.5))+
          annotate("text", x = 1.5, y = max(res2[[gene]]), label = lbl, parse = FALSE, size = 6) +
          theme(axis.text=element_text(size=18),
                axis.title=element_text(size=16,face="bold")))
}

# for(gene in genes){
#   print(ggboxplot(res2, x = "stage", y = gene,
#                  color = "stage",
#                  add = "jitter"))
# }

for(gene in genes){
  p = ggboxplot(res2, x = "stage", y = gene,
                color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"),
                add = "jitter", size = 2) +
    theme(panel.grid.major = element_blank(), legend.position = "none", 
          axis.text=element_text(size=24), axis.ticks = element_line(size = 2, colour = "black"),
          axis.line = element_line(size = 2, colour = "black"), 
          plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=22,face="bold"))
  my_comparisons <- list( c("normal", "EN"), c("EN", "DCIS"), c("DCIS","IDC"),c("normal", "IDC") )
  print(p + stat_compare_means(comparisons = my_comparisons, size = 2))
}
