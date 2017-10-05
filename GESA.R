DE_EN_DCIS = c(sum(change$EN_DCIS >= 0.5, na.rm = TRUE), sum(change$EN_DCIS <= -0.5, na.rm = TRUE))
DE_DCIS_IDC = c(sum(change$DCIS_IDC >= 0.5, na.rm = TRUE), sum(change$DCIS_IDC <= -0.5, na.rm = TRUE))
DE_EN_IDC = c(sum(change$EN_IDC >= 0.5, na.rm = TRUE), sum(change$EN_IDC <= -0.5, na.rm = TRUE))

diff = data.frame(DE_EN_DCIS,DE_DCIS_IDC, DE_EN_IDC)
rownames(diff) = c("up","down")

genes =c("CTSB","MMP9","CD14", "CD163","B2M","HLA.DMA","HLA.DPB1","HLA.DRB1","HLA.DQB1","NFKBIA","NFKBIE",
             "TNFSF12","TNFRSF1A","TRADD","APOE", "CSF1R","CD47","TGFB2","MYD88","TLR3","TLR4",
         "TLR1","TLR6","TRAF3","TGFA","TGFB1","TGFB3","IFNAR1","IFNAR2","IFNGR1","IFNGR2","STAT1",
         "IL10RA", "IL12A","FOXO3", "DDX18","DDX23","DDX39A","DDX41","DDX56","DHX9"," IFI16","IFI44",
         "ATF6", "USP7", "CASP7", "VEGEFA","GSDMD", "CYLD", "CASP1","PYCARD","RHOA","RIPK1")

genes = c("TLR3","TLR4","TLR1","TLR6","IL10RA","IL12A", "CYLD")
genes = c("CASP1","PYCARD","GSDMD", "TLR3", "TLR4", "DHX9")
length(genes)

change2 = change[change$gene.names %in% genes,c(1,4,5,6)]
change2$gene.names =as.factor(change2$gene.names)
long_Mono = melt(change2, id.vars = "gene.names")
p <- plot_ly(long_Mono, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
    add_lines() %>%  layout(showlegend = FALSE)
p
DE_EN_DCIS_up = change[which(change$EN_DCIS >= 0.5), c(1,9)]
DE_EN_DCIS_down = change[which(change$EN_DCIS <= -0.5), c(1,9)]
DE_DCIS_IDC_up = change[which(change$DCIS_IDC >= 0.5), c(1,10)]
DE_DCIS_IDC_down = change[which(change$DCIS_IDC <= -0.5), c(1,10)]
DE_EN_IDC_up = change[which(change$EN_IDC >= 0.5), c(1,11)]
DE_EN_IDC_down = change[which(change$EN_IDC <= -0.5), c(1,11)]

DE = list(DE_EN_DCIS_up, DE_EN_DCIS_down, DE_DCIS_IDC_up, DE_DCIS_IDC_down, DE_EN_IDC_up, DE_EN_IDC_down)
names(DE) = c("DE_EN_DCIS_up", "DE_EN_DCIS_down", "DE_DCIS_IDC_up","DE_DCIS_IDC_down", "DE_EN_IDC_up", "DE_EN_IDC_down")

DE_EI = list()
for(element in names(DE)){
  DE_EI[[element]] = bitr(DE[[element]]$gene.names, 
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)
}

for(element in names(DE)){
  DE[[element]] = merge(DE[[element]], DE_EI[[element]], by.x = "gene.names", by.y = "SYMBOL")
}
# for(i in names(DE)){
#   write.table(DE[[i]], file = paste(i,".txt",sep = ""), quote = FALSE, row.names = FALSE,col.names = FALSE)
# }



######## prepare for clusterProfiler
DE_EN_DCIS = change[which(change$EN_DCIS >= 0.5 | change$EN_DCIS <= -0.5), c(1,9)]
DE_DCIS_IDC = change[which(change$DCIS_IDC >= 0.5 | change$DCIS_IDC <= -0.5), c(1,10)]
DE_EN_IDC = change[which(change$EN_IDC >= 0.5 | change$EN_IDC <= -0.5), c(1,11)]

DEfull = list(DE_EN_DCIS, DE_DCIS_IDC, DE_EN_IDC)
names(DEfull) = c("DE_EN_DCIS", "DE_DCIS_IDC", "DE_EN_IDC")

lapply(DEfull, head)


#stack the elements of the list
gnames.change = data.frame()
for(name in names(DEfull)){
  element = data.frame(DEfull[[name]][,1:2], comparison = sub("DE_","",name))
  colnames(element) = c("gene.name","fold.change","comparison")
  gnames.change = rbind(gnames.change, element)
}
gnames.change = mutate(gnames.change, diff = "upregulated")
gnames.change$diff[gnames.change$fold.change <= -0.5] = "downregulated"


gene.df <- bitr(gnames.change$gene.name, fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
data_cC = merge(gnames.change, gene.df, by.x = "gene.name", by.y = "SYMBOL")


formula_res <- compareCluster(ENTREZID~comparison+diff, data= data_cC, 
                              fun="enrichGO", 
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05)

head(as.data.frame(formula_res))
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)
?compareCluster
