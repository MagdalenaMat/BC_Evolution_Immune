normal = read.table("/home/magda/Desktop/cibersortX/bulk_normal_unames.txt", sep = "\t", header = TRUE)
EN = read.table("/home/magda/Desktop/cibersortX/bulk_EN_unames.txt", sep = "\t", header = TRUE)
DCIS = read.table("/home/magda/Desktop/cibersortX/bulk_DCIS_unames.txt", sep = "\t", header = TRUE)
IDC = read.table("/home/magda/Desktop/cibersortX/bulk_IDC_unames.txt", sep = "\t", header = TRUE)

#Im_genes = c("CTSB","MMP9","CD14", "CD163","HLA-DMA","HLA-DPB1","HLA-DRB1","HLA-DQB1","NFKBIA",
#             "TNFSF12","TNFRSF1A","TRADD","CD80","APOE", "CSF1R","CD47","TGFB2","MYD88","TLR3","IFNAR1","IFNAR2","IFNGR1","IFNGR2","STAT1")

Im_genes = c("CD14", "CD163","NFKBIA")
stages = list(EN,DCIS,IDC)
names(stages) = c("EN","DCIS","IDC")

#select rows wth desired genes and set rownames
for(i in 1:length(stages)){
  rownames(stages[[i]])= stages[[i]][,1]
  stages[[i]]= stages[[i]][rownames(stages[[i]]) %in% Im_genes,-1]
}

gene_names = rownames(stages[[1]])

#stack dataframes together to have it in good format for plotting 
res2 = data.frame()
for(i in 1:length(stages)){
  nstag = data.frame(t(stages[[i]]), stage = names(stages)[i])
  res2 = rbind(res2,nstag)
}


library(ggplot2)
library(ggpubr)
genes = colnames(res2)[1:3]

###t-test for two genes comparisone
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


#multiple genes
for(gene in genes){
  p = ggboxplot(res2, x = "stage", y = gene,outlier.colour = NA,
                color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"), 
                size = 2)+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 2) +
    theme(panel.grid.major = element_blank(), legend.position = "none", 
          axis.text=element_text(size=24), axis.ticks = element_line(size = 2, colour = "black"),
          axis.line = element_line(size = 2, colour = "black"), 
          plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=22,face="bold"))
  my_comparisons <- list( c("EN", "DCIS"), c("DCIS","IDC"),c("EN", "IDC") )
  print(p + stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 2))
}


p= ggboxplot(res2, x = "stage", y = gene,outlier.colour = NA,
              color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"), size =2, add = "jitter") +
  theme(panel.grid.major = element_blank(), legend.position = "none", 
        axis.text=element_text(size=24), axis.ticks = element_line(size = 2, colour = "black"),
        axis.line = element_line(size = 2, colour = "black"), 
        plot.title = element_text(hjust = 0.5), 
        axis.title=element_text(size=22,face="bold"))
my_comparisons <- list( c("EN", "DCIS"), c("DCIS","IDC"),c("EN", "IDC") )
print(p + stat_compare_means(comparisons = my_comparisons, label = "p.signif", size = 2))
