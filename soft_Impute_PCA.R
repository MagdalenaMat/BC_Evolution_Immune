###########################################################################################################################################
path = "/home/magda/Desktop/cibersortX/run2/GEP"
setwd(path)
files <- list.files(path=path)

#extract genes in monocyte signature common for all samples 
# LM22genes_ex = data.frame(matrix(NA, nrow = 5372, ncol = 1))
# for (i in 1: length(files)){
#   r = read.csv(files[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
#   stage = r[r[,1] %in% rownames(al2_mono),] 
#   colnames(stage) = paste(colnames(stage),rep(files[i],10))
#   LM22genes_ex = cbind(LM22genes_ex,stage)
# }
# 
# all(as.character(LM22genes_ex[,2])==as.character(LM22genes_ex[,12]))
# rownames(LM22genes_ex) = as.character(LM22genes_ex[,2])
# LM22genes_ex1 = LM22genes_ex[, -c(1,seq(2, length(LM22genes_ex), 10))]
# 
# LM22genes_ex1[LM22genes_ex1 <=0] = NA # strong assumption that 0's are not detected TODO

library(softImpute)

#soft_imputed = soft_imputed_ex$u %*% diag(soft_imputed_ex$d) %*% t(soft_imputed_ex$v)

# soft_imputed_ex = softImpute(as.matrix(LM22genes_ex1),rank.max=2,maxit=1000,thresh = 1e-5)
# plot(soft_imputed_ex$u)
# plot(soft_imputed_ex$v)
# soft_imputed = soft_imputed_ex$u %*% diag(soft_imputed_ex$d) %*% t(soft_imputed_ex$v)

############
sub_to_SI = list()
for (i in 1: length(files)){
  r = read.csv(files[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
  stage = r[r[,1] %in% rownames(al2_mono),] 
  colnames(stage) = paste(colnames(stage),rep(files[i],10))
  sub_to_SI[[i]] = stage
}

for(i in 1: length(files)){
  sub_to_SI[[i]][sub_to_SI[[i]] <=0] = NA # strong assumption that 0's are not detected TODO
  rownames(sub_to_SI[[i]]) = sub_to_SI[[i]][,1]
  sub_to_SI[[i]] = sub_to_SI[[i]][,-1]
}

# for(i in 1: length(files)){
#   sub_to_SI[[i]] = t(t(sub_to_SI[[i]]) - colMeans(sub_to_SI[[i]],na.rm = TRUE))
# }


#log2
for(i in 1: length(files)){
  sub_to_SI[[i]] = log2(sub_to_SI[[i]]+1) 
}

#gene center
for(i in 1: length(files)){
  sub_to_SI[[i]] = sub_to_SI[[i]] - rowMeans(sub_to_SI[[i]],na.rm = TRUE)
}

# #/sd
# for(i in 1: length(files)){
#   sub_to_SI[[i]] = t(t(sub_to_SI[[i]]) / apply(sub_to_SI[[i]],2,sd,na.rm = TRUE))
# }

sI_expr = list()
for(i in 1: length(files)){
  sI_expr[[i]] = softImpute(as.matrix(sub_to_SI[[i]]),rank.max=2,maxit=1000,thresh = 1e-5)
}


v_list = list()
for(i in 1: length(sI_expr)){
  v = data.frame(sI_expr[[i]]$v)
  colnames(v) = c("dim1","dim2")
  rownames(v) = colnames(sub_to_SI[[i]])
  v_list[[i]] = v
}


sI_Mono = data.frame()
for(i in 1: length(v_list)){
  mon = v_list[[i]][5,]
  sI_Mono = rbind(sI_Mono,mon)
}
rownames(sI_Mono) = c("LN","normal","EN","DCIS","IDC","AVL","metECE")

plot_ly(data = sI_Mono,  x = ~dim1, y = ~dim2,  
        text = rownames(sI_Mono), type="scatter", mode="markers") %>%
  add_annotations(rownames(sI_Mono)) %>%
  layout(title = "log2 gene centered")

names_stages = c("LN","normal","EN","DCIS","IDC","AVL","metECE")
p = plot_ly()
for(i in 1: length(v_list)){
  p = p %>% add_trace(data= v_list[[i]], x = ~dim1, y = ~dim2, name = names_stages[i],
                      text = rownames(v_list[[i]]), type="scatter", mode="markers") %>%
    layout(title = "log2 gene centered")
}
p



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

DCIS_cells = data.frame(sI_expr[[4]]$v)
colnames(DCIS_cells) = c("dim1","dim2")
rownames(DCIS_cells) = colnames(sub_to_SI[[4]])
IDC_cells = data.frame(sI_expr[[5]]$v)
colnames(IDC_cells) = c("dim1","dim2")
rownames(IDC_cells) = colnames(sub_to_SI[[5]])
plot_ly() %>% 
  add_trace(data=DCIS_cells, x = ~dim1, y = ~dim2,  text = rownames(DCIS_cells), type="scatter", mode="markers") %>% 
  add_trace(data=IDC_cells, x = ~dim1, y = ~dim2,  text = rownames(IDC_cells), type="scatter", mode = "markers")