#get design info

path = "/home/magda/Desktop/breast_cancer_progression/Megan/data"
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
