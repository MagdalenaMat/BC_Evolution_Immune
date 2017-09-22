DE_EN_DCIS = c(sum(change$EN_DCIS >= 0.5, na.rm = TRUE), sum(change$EN_DCIS <= -0.5, na.rm = TRUE))
DE_DCIS_IDC = c(sum(change$DCIS_IDC >= 0.5, na.rm = TRUE), sum(change$DCIS_IDC <= -0.5, na.rm = TRUE))
DE_EN_IDC = c(sum(change$EN_IDC >= 0.5, na.rm = TRUE), sum(change$EN_IDC <= -0.5, na.rm = TRUE))

diff = data.frame(DE_EN_DCIS,DE_DCIS_IDC, DE_EN_IDC)
rownames(diff) = c("up","down")

genes =c("CTSB","MMP9","CD14", "CD163","HLA.DMA","HLA.DPB1","HLA.DRB1","HLA.DQB1","NFKBIA","NFKBIE",
             "TNFSF12","TNFRSF1A","TRADD","B2M","APOE", "CSF1R","CD47","TGFB2","MYD88","TLR3","TLR4",
         "TLR1","TLR6","TRAF3","TGFA","TGFB1","TGFB3","IFNAR1","IFNAR2","IFNGR1","IFNGR2","STAT1",
         "IL10RA", "IL12A","FOXO3", "DDX18","DDX23","DDX39A","DDX41","DDX56","DHX9"," IFI16","IFI44",
         "ATF6", "USP7", "CASP7", "VEGEFA","GSDMD", "CYLD", "CASP1","PYCARD","RHOA","RIPK1")
length(genes)

change2 = change[change$gene.names %in% genes,]

DE_EN_DCIS_up = change[which(change$EN_DCIS >= 0.5), 1]
DE_EN_DCIS_down = change[which(change$EN_DCIS <= -0.5), 1]
DE_DCIS_IDC_up = change[which(change$DCIS_IDC >= 0.5), 1]
DE_DCIS_IDC_down = change[which(change$DCIS_IDC <= -0.5), 1]
DE_EN_IDC_up = change[which(change$EN_IDC >= 0.5), 1]
DE_EN_IDC_down = change[which(change$EN_IDC <= -0.5), 1]

DE = list(DE_DCIS_IDC_up,DE_DCIS_IDC_down, DE_EN_DCIS_up, DE_EN_DCIS_down, DE_EN_IDC_up, DE_EN_IDC_down)
names(DE) = c("DE_DCIS_IDC_up","DE_DCIS_IDC_down", "DE_EN_DCIS_up", "DE_EN_DCIS_down", "DE_EN_IDC_up", "DE_EN_IDC_down")

for(i in names(DE)){
  write.table(DE[[i]], file = paste(i,".txt",sep = ""), quote = FALSE, row.names = FALSE,col.names = FALSE)
}
