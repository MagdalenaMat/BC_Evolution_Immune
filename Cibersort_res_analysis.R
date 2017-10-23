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
#path = "/home/magda/Desktop/cibersortX/boot.files/boot.run/CVfiltered"
setwd(path)
files <- list.files(path=path)

#extract names
namesGEP = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_","_unames.txt_GEPs_CVFilteredAll.txt"),c("",""),i)
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
Mono_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  Mono = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,6)] #select gene_names and Monocytes
  Mono_ex = cbind(Mono_ex,Mono)
}

all(as.character(Mono_ex[,2])==as.character(Mono_ex[,4]))
rownames(Mono_ex) = as.character(Mono_ex[,2])
Mono_ex = Mono_ex[, seq(3, length(Mono_ex), 2)]
colnames(Mono_ex) = c("AVL","DCIS","EN","IDC","LN","MetECE","normal")
#colnames(Mono_ex) = namesGEP

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
colnames(to_reshape)
to_reshape = to_reshape[,c(1,8,6,4,3,5,2,7)]

#compute log2 Fold Changes
change = mutate(to_reshape, EN_DCIS = log2(to_reshape$DCIS/to_reshape$EN), 
                DCIS_IDC = log2(to_reshape$IDC/to_reshape$DCIS),
                EN_IDC = log2(to_reshape$IDC/to_reshape$EN))
change$gene.names = sub("HLA.","HLA-", change$gene.names) 
change$gene.names = sub("NKX2.","NKX2-", change$gene.names) 
change$gene.names = sub("WT1.","WT1-", change$gene.names)

ImGenes = list()
ImGenes[["TLR_NFkB"]] = c("TLR3","TLR4","TLR1","TLR6","IL10RA","IL12A", "CYLD", "NFKBIA","NFKBIE", "MYD88")
ImGenes[["NLRs"]] = c("CASP1","PYCARD","GSDMD", "TLR3", "TLR4", "DHX9")
ImGenes[["HLA"]] = c("HLA.DMA","HLA.DPB1","HLA.DRB1","HLA.DQB1")
ImGenes[["TNF"]] = c("TNFSF12","TNFRSF1A","TRADD")
ImGenes[["NFkB"]] = c("NFKBIA","MYD88")
ImGenes[["B2M"]] = c("B2M")
ImGenes[["CD"]] = c("CD14","CD163")
ImGenes[["IFN"]] = c("IFNGR1","IFNGR2","NFKBIA")
ImGenes[["selected"]] = c("CD14", "CD163","NFKBIA")
genes = c("TLR_NFkB","NLRs","CD","HLA","IFN","TNF","NFkB", "B2M","selected")  
#genes = c("CD", "HLA","IFN", "selected")

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


########for bootstraped

gnames = unlist(ImGenes)
to_reshape = to_reshape[to_reshape[,1] %in% gnames,]
to_plot = list()
for(i in c("EN","DCIS","IDC")){
  to_plot[[i]] = as.matrix(to_reshape[,grepl(i,colnames(to_reshape))])
}

#stack dataframes together to have it in good format for plotting 
res.boot = data.frame()
for(i in 1:length(to_plot)){
  nstag = data.frame(t(to_plot[[i]]), stage = names(to_plot)[i])
  res.boot = rbind(res.boot,nstag)
}

for(gene in gnames){
  p = ggboxplot(res.boot, x = "stage", y = gene,outlier.colour = NA,
                color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"))+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
    theme(panel.grid.major = element_blank(), legend.position = "none")
  my_comparisons <- list( c("EN", "DCIS"), c("DCIS","IDC"),c("EN", "IDC"))
  print(p + stat_compare_means(comparisons = my_comparisons, size = 1))
}

for(gene in gnames){
  p = ggboxplot(res.boot, x = "stage", y = gene,outlier.colour = NA,
                color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"))+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
    theme(panel.grid.major = element_blank(), legend.position = "none")
  my_comparisons <- list(c("normal", "EN"), c("EN", "DCIS"), c("DCIS","IDC"),c("normal", "DCIS"),c("EN","IDC"),c("normal", "IDC"))
  print(p + stat_compare_means(comparisons = my_comparisons, size = 1))
}


#####for data from 100x bootstrapp
library(reshape2)
library(dplyr)
library(ggplot2)

boot.immuno = read.table("/home/magda/Desktop/bootstrap100x_plots/long_res_boot_100_ImGenes.txt", sep = "\t", header = TRUE, row.names = 1)
boot.immuno = boot.immuno[boot.immuno$stage %in% c("EN","DCIS","IDC"),]
boot.immuno = boot.immuno[,colnames(boot.immuno) %in% c("CD14","CD163","NFKBIA","PYCARD","IL10RA","stage")]
long_Bimmuno = melt(boot.immuno, id.vars = "stage")
names(long_Bimmuno) = c("stage","gene.name","expression")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
df1 <- data_summary(long_Bimmuno, varname="expression", 
                    groupnames=c("stage", "gene.name"))
df1$gene.name=as.factor(df1$gene.name)
df1$stage = factor(df1$stage, levels = c("EN","DCIS","IDC"))
p <- ggplot(df1, aes(x=gene.name, y=expression, fill=stage)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=expression-sd, ymax=expression+sd), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + 
  scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red"))+
  theme_classic() +theme(axis.text=element_text(size=24, colour = "black"), axis.ticks = element_line(size = 2, colour = "black"),
                         axis.line = element_line(size = 2, colour = "black"), 
                         plot.title = element_text(hjust = 0.5), 
                         axis.title=element_text(size=22,face="bold"), 
                         legend.text = element_text(size=20), legend.title = element_text(size = 20)) 


#####################boxplots
library(ggpubr)
boot.immuno$stage = factor(boot.immuno$stage, levels = c("EN","DCIS","IDC"))
gnames = c("CD14","CD163","NFKBIA","PYCARD","IL10RA")
for(gene in gnames){
  p = ggboxplot(boot.immuno, x = "stage", y = gene,outlier.colour = NA,
                color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"))+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
    theme(panel.grid.major = element_blank(), legend.position = "none")
  my_comparisons <- list( c("EN", "DCIS"), c("DCIS","IDC"),c("EN", "IDC"))
  print(p + stat_compare_means(comparisons = my_comparisons, size = 1))
}

