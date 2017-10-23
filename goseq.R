###goseq
to_goseq = list()
for(nm in tail(names(change),11)){
  to_goseq[[paste0("DE_",nm,"_up")]] = as.character(change[change[[nm]] >= 0.5, 1])
  to_goseq[[paste0("DE_",nm,"_down")]] = as.character(change[change[[nm]] <= -0.5, 1])
}

to_goseq_char = list()
for(i in names(DE)){
  to_goseq_char[[i]] = as.integer(to_goseq[[i]] %in% DE[[i]])
  names(to_goseq_char[[i]]) = to_goseq[[i]]
}

lengths(to_goseq)
lengths(to_goseq_char)
supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]

# library(goseq)
# go_res = list()
# for(genes in names(to_goseq_char)){
#   pwf=nullp(to_goseq_char[[genes]],"hg19","geneSymbol")
#   go_res[[genes]] = goseq(pwf, "hg19","geneSymbol", test.cats=c("GO:MF"), method = "Hypergeometric")
# }
# 
# lengths(go_res)
# View(go_res$DE_EN_DCIS_up)