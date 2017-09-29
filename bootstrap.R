for(i in c("bulk_normal_unames.txt","bulk_EN_unames.txt","bulk_DCIS_unames.txt","bulk_IDC_unames.txt")){
  for(j in 1:10){
    setwd("~/Desktop/cibersortX")
    data = read.table(i,header=T, sep="\t", row.names=1)
    n = ncol(data)
    boot.mixture = data[,sample(x = 1:n, size = n, replace = TRUE)]
    boot.mixture = cbind(genes = rownames(boot.mixture),boot.mixture)
    setwd("~/Desktop/cibersortX/boot.files")
    write.table(boot.mixture, paste("boot_",j, i, sep = ""), sep = "\t",row.names = FALSE)
    }
}

for(i in Sys.glob("*unames.txt")[1]){ 
  test = read.table(i,header=T, sep="\t", row.names=1)
}
str(test)
