#############
path = "/home/magda/Desktop/cibersortX/run2/CV_filtered"
setwd(path)
files <- list.files(path=path)

#extract all gene names from the files and put them in the list to find intersect 
#of genes acros all the histological subtypes
results = list()
for (i in 1: length(files)){
  r = read.table(files[i], sep = "\t", header = TRUE)
  results[[i]] = as.character(r[,1])
}

universe = Reduce(intersect, results)

#extract genes in monocyte signature common for all samples 
Mono_ex = data.frame(matrix(NA, nrow = 9748, ncol = 1))
for (i in 1: length(files)){
  r = read.table(files[i], sep = "\t", header = TRUE)
  Mono = r[r[,1] %in% universe,c(1,6)] #select gene_names and Monocytes
  Mono_ex = cbind(Mono_ex,Mono)
}

all(as.character(Mono_ex[,2])==as.character(Mono_ex[,4]))
rownames(Mono_ex) = as.character(Mono_ex[,2])
Mono_ex = Mono_ex[, seq(3, length(Mono_ex), 2)]
colnames(Mono_ex) = c("AVL_Mono","DCIS_Mono","EN_Mono","IDC_Mono","LN_Mono","MetECE_Mono","normal_Mono")

Mono_ex[Mono_ex<=1] <- NA
apply(Mono_ex, 2, function(x){sum(!is.na(x))})

Mono_ex_DCIS = Mono_ex[!is.na(Mono_ex$DCIS_Mono),] 
cor(Mono_ex,use="pairwise.complete.obs",method="pearson") # can't use spearman here
al3_mono = Mono_ex[rowSums(!is.na(Mono_ex)) >= 3,]

###########
#plot gene values in different stages 
###########

library(reshape2)
library(dplyr)
to_reshape = data.frame(gene.names = rownames(al3_mono), al3_mono)
to_reshape = to_reshape[,c(1,8,6,4,3,5,2,7)]
change = mutate(to_reshape, EN_DCIS = ((to_reshape$EN_Mono/to_reshape$DCIS_Mono)-1)*100, 
                DCIS_IDC = ((to_reshape$IDC_Mono/to_reshape$DCIS_Mono)-1)*100)

long_Mono = melt(to_reshape, id.vars = "gene.names")
p <- plot_ly(long_Mono, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>%  layout(showlegend = FALSE)


EN_subset = to_reshape[!is.na(to_reshape$EN_Mono),c(1,4,5,6)]
long_EN = melt(EN_subset, id.vars = "gene.names")

p <- plot_ly(long_EN, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>%  layout(showlegend = FALSE)  



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
fold = mutate(fold, DCIS_EN = log2((EN_Mono+1)/(DCIS_Mono +1)), IDC_EN = log2((EN_Mono+1)/(IDC_Mono +1)))
fold = fold[order(fold[,5]),]
filter = c("IFNA","IFNG") #"DHX","TRIM"
viral_res = fold[grep(paste(filter,collapse="|"),fold[,1]),]
rownames(viral_res) = viral_res[,1]
viral_res = viral_res[,-1]
heatmaply(log2(viral_res[,1:3]), dendrogram = "none", 
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "darkblue", high = "red", midpoint = 8.5),
          margins = c(90,80,NA,0))

trimmed1 = trimmed[complete.cases(trimmed[,2]),]
heatmaply(log2(trimmed1+1))

#######
path = "/home/magda/Desktop/cibersortX/run2/GEP"
setwd(path)
files <- list.files(path=path)

#extract genes in monocyte signature common for all samples 
LM22genes_ex = data.frame(matrix(NA, nrow = 2696, ncol = 1))
for (i in 1: length(files)){
  r = read.csv(files[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
  stage = r[r[,1] %in% rownames(al3_mono),] 
  colnames(stage) = paste(colnames(stage),rep(files[i],10))
  LM22genes_ex = cbind(LM22genes_ex,stage)
}

all(as.character(LM22genes_ex[,2])==as.character(LM22genes_ex[,12]))
rownames(LM22genes_ex) = as.character(LM22genes_ex[,2])
LM22genes_ex1 = LM22genes_ex[, -c(1,seq(2, length(LM22genes_ex), 10))]

LM22genes_ex1[LM22genes_ex1 <=0] = NA # strong assumption that 0's are not detected TODO

library(softImpute)

soft_imputed = soft_imputed_ex$u %*% diag(soft_imputed_ex$d) %*% t(soft_imputed_ex$v)

soft_imputed_ex = softImpute(as.matrix(LM22genes_ex1),rank.max=2,maxit=1000,thresh = 1e-5)
plot(soft_imputed_ex$u)
plot(soft_imputed_ex$v)
soft_imputed = soft_imputed_ex$u %*% diag(soft_imputed_ex$d) %*% t(soft_imputed_ex$v)

############
sub_to_SI = list()
for (i in 1: length(files)){
  r = read.csv(files[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
  stage = r[r[,1] %in% rownames(al3_mono),] 
  colnames(stage) = paste(colnames(stage),rep(files[i],10))
  sub_to_SI[[i]] = stage
}

for(i in 1: length(files)){
  sub_to_SI[[i]][sub_to_SI[[i]] <=0] = NA # strong assumption that 0's are not detected TODO
  rownames(sub_to_SI[[i]]) = sub_to_SI[[i]][,1]
  sub_to_SI[[i]] = sub_to_SI[[i]][,-1]
}

sI_expr = list()
for(i in 1: length(files)){
  sI_expr[[i]] = softImpute(as.matrix(LM22genes_ex1),rank.max=2,maxit=1000,thresh = 1e-5)
}

plot(sI_expr[[4]]$u)
points(sI_expr[[5]]$u, col = 2)

library(plotly)
DCIS = data.frame(sI_expr[[4]]$u)
colnames(DCIS) = c("dim1","dim2")
rownames(DCIS) = rownames(sub_to_SI[[4]])
IDC = data.frame(sI_expr[[5]]$u)
colnames(IDC) = c("dim1","dim2")
rownames(IDC) = rownames(sub_to_SI[[5]])


plot_ly() %>% 
  add_trace(data=DCIS, x = ~dim1, y = ~dim2,  text = rownames(DCIS), type="scatter", mode="markers") %>% 
  add_trace(data=IDC, x = ~dim1, y = ~dim2,  text = rownames(IDC), type="scatter", mode = "markers")
