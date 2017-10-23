
library("Rgraphviz")
library("DOSE")
library("clusterProfiler")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

gene = read.table(file.choose(), sep = "\t")
gene = as.character(gene$V1)
gene1 = gene[!gene %in% c("HLA.B","HLA.C", "HLA.DPA1","NKX2.8")] 

gene.df <- bitr(gene1, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

# ggo <- groupGO(gene     = gene.df$ENTREZID,
#                OrgDb    = org.Hs.eg.db,
#                ont      = "BP",
#                level    = 3,
#                readable = TRUE)
# 
# head(ggo)

ego2 <- enrichGO(gene         = gene.df$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable = TRUE)
head(ego2)
ego2@result$Description


barplot(ego2, showCategory=22)
dotplot(ego2, showCategory=22)
enrichMap(ego2, n= 22)
cnetplot(ego2, categorySize="pvalue") #foldChange= named vector
plotGOgraph(ego2)

path = "/home/magda/Desktop/breast_cancer_progression/DE_gene_lists"
setwd(path)
files <- list.files(path=path)

DE = list()
for(file in files){
  gene.list = read.table(file, sep ="\t")
  gene.list = as.character(gene.list$V1)
  gene.df <- bitr(gene.list, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  DE[[file]] = gene.df
}
lapply(DE, head)

entrezIDs = list()
for(name in names(DE)){
  entrezIDs[[name]] = DE[[name]]$ENTREZID
}
lapply(entrezIDs, head)

##########################
gmtfile <- file.choose()
c7 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c7)
head(egmt)

egmtc7 = list()
for(element in names(entrezIDs)){
  egmtc7[[element]]  = enricher(entrezIDs[[element]], TERM2GENE=c7)
}

for(element in names(egmtc7)){
  print(dotplot(egmtc7[[element]]))
}
#############
data(gcSample)#list with genes chr vectors
lapply(gcSample, head)

ck = compareCluster(geneCluster = entrezIDs, 
                    fun = "enrichGO", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

head(as.data.frame(ck))


data(geneList, package="DOSE")
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))
dotplot(ck)
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)
