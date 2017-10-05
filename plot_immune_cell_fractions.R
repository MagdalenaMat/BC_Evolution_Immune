##fractions 

EN = read.table("/home/magda/Desktop/cibersortX/run2/CIBERSORTx_bulk_EN_unames.txt_Fractions_PostComBat.txt",
                quote = "", row.names = 1, header = T, sep = "\t")


DCIS = read.table("/home/magda/Desktop/cibersortX/run2/CIBERSORTx_bulk_DCIS_unames.txt_Fractions_PostComBat.txt",
                  quote = "", row.names = 1, header = T, sep = "\t")


IDC = read.table("/home/magda/Desktop/cibersortX/run2/CIBERSORTx_bulk_IDC_unames.txt_Fractions_PostComBat.txt",
                        quote = "", row.names = 1, header = T, sep = "\t")

colnames(DCIS)[13:16]

stages = list(EN,DCIS,IDC)
names(stages) = c("EN","DCIS","IDC")

#stack dataframes together to have it in good format for plotting 
res2 = data.frame()
for(i in c("EN","DCIS", "IDC")){
  nstag = data.frame(stages[[i]][,13:16], stage = names(stages[i]))
  res2 = rbind(res2,nstag)
}

cells = colnames(res2)[1:4]
for(cell in cells){
  p = ggboxplot(res2, x = "stage", y = cell,outlier.colour = NA,
                color = "stage")+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
    theme(panel.grid.major = element_blank(), legend.position = "none")
  my_comparisons <- list(c("EN","DCIS"),c("DCIS", "IDC"),c("EN","IDC"))
  print(p + stat_compare_means(comparisons = my_comparisons, size = 1))
}

