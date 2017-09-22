for(i in c("bulk_EN_unames.txt","bulk_DCISunames.txt","bulk_IDC_unames.txt")){
  i = read.table(i,header=T, sep="\t", row.names=1)
  for(j in 10){
    n = ncol(i)
    boot.mixture = i[sample(x = 1:n, size = n, replace = TRUE),]
    write.table(boot.mixture, paste("boot_",j, i, sep = ""), row.names = TRUE)
    mixture = 
    sigmatrix = "LM22.txt"
    classes = "GEPs_LM9.txt"
    Monocytes = "Monocytes.txt"
    
    CIBERSORTxGEP(mixture, sigmatrix, classes, cibresults = NA, label= i, groundtruth= Monocytes, 
                  maxsamples=NA , degclasses = "", doComBat = TRUE, referenceprofiles = "LM22_source_GEPs.txt", 
                  dobg = FALSE, redo=TRUE, threads="",plots=FALSE)}
}

for(i in c("bulk_EN_unames.txt","bulk_DCIS_unames.txt","bulk_IDC_unames.txt")){
  for(j in 1:10){
    setwd("~/Desktop/cibersortX")
    data = read.table(i,header=T, sep="\t", row.names=1)
    n = ncol(data)
    boot.mixture = data[,sample(x = 1:n, size = n, replace = TRUE)]
    setwd("~/Desktop/cibersortX/boot_files")
    write.table(boot.mixture, paste("boot_",j, i, sep = ""), sep = "\t",row.names = TRUE)
    }
}


