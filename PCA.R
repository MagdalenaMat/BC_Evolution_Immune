#prepare samples accross diferent experiments for PCA

path = "/home/magda/Desktop/breast_cancer_progression/data"
setwd(path)
files <- list.files(path=path)

design.df = data.frame(name = NA, batch = NA)
for(i in 1:length(files)){
  dat = read.csv(files[i])
  dat = dat[,-c(1,2)]
  batch.info = data.frame(name = colnames(dat), batch = files[i])
  design.df = rbind(design.df, batch.info)
}
design.df = design.df[-1,]

#check if gene names are same and in the same order in all files  
names = list()
for(i in 1:length(files)){
  print(i)
  namb = as.character(read.csv(files[i])[,2])
  names[[i]] = namb
}
for(i in 1:(length(names)-1) ){
  print(all(names[[i]][1] == names[[c(i+1)]][1]))
}


gene.names = as.character(read.csv(files[1])[,1])

names <- data.frame(matrix(gene.names,58037,1))
dim(names)

results = list()
desine = data.frame(name = NA, batch = NA, subtype = NA)
for(j in c("DCIS", "EN", "IDC", "AVL", "LN", "ECE", "normal")){
  for (i in 1: length(files)){
    r = read.csv(files[i])
    cols = as.matrix(r[,grepl(j,colnames(r)), drop = FALSE])
    if(ncol(cols) != 0){
      results[[j]] = cbind(results[[j]],cols)
      des = data.frame(name = colnames(cols), batch = files[i], subtype = j)
      desine = rbind(desine, des)
    }
    print(paste(j,files[i],ncol(cols), colnames(cols)))
  }
  results[[j]] = data.frame(GeneNames = gene.names, results[[j]])
  print(ncol(results[[j]]))
}

desine = desine[-1,]
dim(desine)
sum(lengths(results))

#disolve EC
ECE = results[["ECE"]]
results[["ECE"]] = ECE[,!grepl("et",colnames(ECE))]
results[["met_ECE"]] = ECE[,c(1,grep("et",colnames(ECE)))]

#delete scDCIS
dim(results[["DCIS"]])
results[["DCIS"]] = results[["DCIS"]][,!grepl("DCIS_single|DCIS_ablation", colnames(results[["DCIS"]]))]

desine1 = desine[!grepl("DCIS_single|DCIS_ablation",desine$name),]
desine1$subtype[grepl("et", desine1$name)] = "met_ECE"
dim(desine1)

#annotate
library(AnnotationDbi)
library(org.Hs.eg.db)
anot.counts = list()
for(i in names(results)){
  mixture = results[[i]]
  genes = as.character(mixture$GeneNames)
  genes.list <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  list1 = data.frame(lapply(genes.list, type.convert), stringsAsFactors=FALSE)
  list2 = t(list1)
  list2 = data.frame(list2)
  list2$ens_id = rownames(list2)
  ex1 = merge(list2, mixture, by.x="ens_id", by.y="GeneNames")
  ex1$ens_id = NULL
  ex2 = ex1[complete.cases(ex1),]
  colnames(ex2)[1] = "genes"
  ex3 = ex2[!duplicated(ex2[,1]),]
  anot.counts[[i]] = ex3
}

##############do PCA on results

for(i in 1:(length(anot.counts)-1) ){
  print(all(anot.counts[[i]][1] == anot.counts[[c(i+1)]][1]))
}

#nam.fac = as.character()
breast.data = data.frame(gene.names = anot.counts[[1]][,1])
for(name in names(anot.counts)){
  bdat = anot.counts[[name]]
  bdat = bdat[,2:length(colnames(bdat))]
  breast.data = cbind(breast.data, bdat)
  #nam.fac = c(nam.fac, rep(name, length(colnames(bdat))))
}


rownames(breast.data) = breast.data[,1]
breast.data = breast.data[,-1]

breast.data = breast.data[rowSums(breast.data) > 300,]
breast.data = as.matrix(breast.data)

colnames(breast.data)[is.na(pmatch(colnames(breast.data), desine1$name))]

desine2 = desine1[match(colnames(breast.data), desine1$name),]
desine2$name == colnames(breast.data)
desine2$batch = factor(desine2$batch)
desine2$batch = sub(".raw_reads.csv","", desine2$batch)
desine2$subtype = factor(desine2$subtype, levels = c("normal","EN","DCIS","IDC","ECE","met_ECE", "LN", "AVL"))

setwd("~/Desktop/breast_cancer_progression")
write.table(desine2, "desine.BC.txt", sep = "\t")
write.table(breast.data, "breast.data.combined.txt", sep = "\t", row.names = TRUE)
# nam.fac = factor(nam.fac, levels = c("normal","EN","DCIS","IDC","ECE","met_ECE", "LN", "AVL"))
# design1 = data.frame(samples = colnames(breast.data),subtype = nam.fac)
# design2 = merge(design1, design.df, by.x = "samples", by.y = "name" )
# design1$samples[is.na(pmatch(design1$samples,design.df$name))]
# design$Id[is.na(pmatch(design$Id,colnames(data.ordered)))]

library(DESeq2)
library(plotly)
dds <- DESeqDataSetFromMatrix(countData = breast.data,
                              colData = desine2,
                              design= ~ subtype)


vsd <- vst(dds, blind=FALSE)

PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
p <- plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = "Set2", marker = list(size = 20)) %>%
  layout(title = "DESeq2::vst transformed original data")
p


PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p <- plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = "Set2", marker = list(size = 20)) %>%
  layout(title = "DESeq2::vst transformed original data")
p
