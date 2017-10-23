####samples

path = "/home/magda/Desktop/breast_cancer_progression/data"
setwd(path)
files <- list.files(path=path)

#check if gene names are sme and in the same order in all files  
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
for(j in c("DCIS", "EN", "IDC", "AVL", "LN", "ECE", "normal")){
  for (i in 1: length(files)){
    r = read.csv(files[i])
    cols = as.matrix(r[,grepl(j,colnames(r)), drop = FALSE])
    if(ncol(cols) != 0){
      results[[j]] = cbind(results[[j]],cols)
    }
    print(paste(j,files[i],ncol(cols), colnames(cols)))
  }
  results[[j]] = data.frame(GeneNames = gene.names, results[[j]])
  print(ncol(results[[j]]))
}

#disolve EC
ECE = results[["ECE"]]
results[["ECE"]] = ECE[,!grepl("et",colnames(ECE))]
results[["met_ECE"]] = ECE[,c(1,grep("et",colnames(ECE)))]

#delete scDCIS
results[["DCIS"]] = results[["DCIS"]][,c(1:120,124:129)]

design = as.character()
for(name in names(results)){
  design = c(design, rep(name,ncol(results[[name]])-1))
}


# calculate TPM
TPM = list()
for(i in c("DCIS", "EN", "IDC", "AVL", "LN", "ECE", "normal", "met_ECE")){
  colSums = apply(results[[i]][,-1],2,sum)
  TPM[[i]] = sweep(results[[i]][,-1]*1000000,2,colSums,"/")
  TPM[[i]] = cbind(gene = results[[i]][,1], TPM[[i]])
}
#annotate
library(AnnotationDbi)
library(org.Hs.eg.db)
for(i in c("DCIS", "EN", "IDC", "AVL", "LN", "ECE", "normal", "met_ECE")){
  mixture =TPM[[i]]
  genes = as.character(mixture$gene)
  genes.list <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  list1 = data.frame(lapply(genes.list, type.convert), stringsAsFactors=FALSE)
  list2 = t(list1)
  list2 = data.frame(list2)
  list2$ens_id = rownames(list2)
  ex1 = merge(list2, mixture, by.x="ens_id", by.y="gene")
  ex1$ens_id = NULL
  ex2 = ex1[complete.cases(ex1),]
  colnames(ex2)[1] = "genes"
  
  setwd("~/Desktop/cibersortX")
  write.table(ex2, file = paste("bulk_", i,".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE) 
  
  #DCIS = read.table(file.choose(), sep = "\t", header = TRUE)
  ex3 = ex2[!duplicated(ex2[,1]),]
  write.table(ex3, file = paste("bulk_",i,"_unames.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

}
lengths(results)


