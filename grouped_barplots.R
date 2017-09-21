#base r barplot
rownames(Immuno) = Immuno$gene.names
Immuno = Immuno[,-1]
Immuno = t(Immuno)

barplot(Immuno,beside = TRUE, col = c("dodgerblue3","orange","chartreuse4","red"),width=0.8)

legend("topleft", inset=.02, title="BC progression stages",
       c("normal","EN","DCIS","IDC"), fill=c("dodgerblue3","orange","chartreuse4","red"), horiz=FALSE, cex=0.8)

#plotly
Immuno = change[change$gene.names %in% ImGenes$selected,]
Immuno = Immuno[,c(1,2,4,5,6)]
Immuno$gene.names = as.factor(as.character(Immuno$gene.names))
plot_ly(Immuno, x = ~gene.names, y = ~normal, type = 'bar', name = 'normal') %>%
  add_trace(y = ~EN, name = 'EN') %>%
  add_trace(y = ~DCIS, name = 'DCIS') %>%
  add_trace(y = ~IDC, name = 'IDC') %>%
  layout(yaxis = list(title = 'normalized expression', titlefont = list(size = 20, color = "black"), tickfont = list(size=15, color = "black")),
         xaxis = list(title = '', tickfont = list(size = 28, color = "black")), 
         legend = list(font = list(size = 28, color = "black")),
         barmode = 'group',margin = list(b = 50))

#ggplot
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
