#select genes with biggest FC difference between EN and IDC
change1 = change[order(change$EN_IDC, decreasing = TRUE),c(1,11)]
up_EN_IDC = change1$gene.names[which(change1$EN_IDC >= 0.5)]
down_EN_IDC = change1$gene.names[which(change1$EN_IDC <= -0.5)]

scaled = change[change$gene.names %in% up_EN_IDC,c(1,4,6)]
scaled$gene.names = as.factor(as.character(scaled$gene.names))
scaled = data.frame(gene.names = scaled[,1], scaled[,c(2,3)]/scaled[,2])
long_scaled = melt(scaled,id.vars = "gene.names")
p <- plot_ly(long_scaled, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>% layout(showlegend = FALSE, margin = list(500,50,50,50))  
p

scaled = change[change$gene.names %in% down_EN_IDC,c(1,4,6)]
scaled$gene.names = as.factor(as.character(scaled$gene.names))
scaled = data.frame(gene.names = scaled[,1], scaled[,c(2,3)]/scaled[,3])
long_scaled = melt(scaled,id.vars = "gene.names")
p <- plot_ly(long_scaled, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>% layout(showlegend = FALSE, margin = list(500,50,50,50))  
p

add_annotations(long_EN_DCIS_up$gene.names) %>%
  
  not_scaled = change[change$gene.names %in% down_EN_IDC,c(1,4,5,6)]
not_scaled$gene.names = as.factor(as.character(not_scaled$gene.names))
long_not_scaled = melt(not_scaled, id.vars = "gene.names")
p <- plot_ly(long_not_scaled, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>% layout(showlegend = FALSE, margin = list(500,50,50,50))  
p


#select genes with biggest FC difference between EN and DCIS
change1 = change[order(change$EN_DCIS, decreasing = TRUE),c(1,9,10)]
To_GOEA_up_EN_DCIS = change1$gene.names[which(change1$EN_DCIS >= 0.5)]
To_GOEA_down_EN_DCIS = change1$gene.names[which(change1$EN_DCIS <= -1)]

##select genes with biggest FC difference between DCIS and IDC
change2 = change[order(change$DCIS_IDC, decreasing = TRUE),c(1,9,10)]
To_GOEA_up_DCIS_IDC = change2$gene.names[which(change2$DCIS_IDC >= 1)]
To_GOEA_down_DCIS_IDC = change2$gene.names[which(change2$DCIS_IDC <= -0.8)]




long_Mono = melt(to_reshape, id.vars = "gene.names")
p <- plot_ly(long_Mono, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>%  layout(showlegend = FALSE)


EN_subset = to_reshape[!is.na(to_reshape$EN_Mono),c(1,4,5,6)]
long_EN = melt(EN_subset, id.vars = "gene.names")

p <- plot_ly(long_EN, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>%  layout(showlegend = FALSE)  

up_EN_DCIS_to_reshape = to_reshape[to_reshape$gene.names %in% To_GOEA,c(1,4,5,6)]
up_EN_DCIS_to_reshape1 = mutate(up_EN_DCIS_to_reshape, difference = up_EN_DCIS_to_reshape$IDC_Mono - up_EN_DCIS_to_reshape$EN_Mono)
up_EN_DCIS_to_reshape1 = up_EN_DCIS_to_reshape1[order(up_EN_DCIS_to_reshape1$difference, decreasing = TRUE),]
up_EN_DCIS_to_reshape1 = up_EN_DCIS_to_reshape1[1:5,c(1:4)]
long_EN_DCIS_up = melt(up_EN_DCIS_to_reshape1, id.vars = "gene.names")
p <- plot_ly(long_EN_DCIS_up, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>% add_annotations(long_EN_DCIS_up$gene.names) %>% layout(showlegend = FALSE)  


down_EN_DCIS_to_reshape = to_reshape[to_reshape$gene.names %in% To_GOEA_down,c(1,4,5,6)]
down_EN_DCIS_to_reshape1 = mutate(down_EN_DCIS_to_reshape, difference = down_EN_DCIS_to_reshape$IDC_Mono - down_EN_DCIS_to_reshape$EN_Mono)
down_EN_DCIS_to_reshape1 = down_EN_DCIS_to_reshape1[order(down_EN_DCIS_to_reshape1$difference, decreasing = TRUE),]
down_EN_DCIS_to_reshape1 = down_EN_DCIS_to_reshape1[1:20,c(1:4)]
long_EN_DCIS_down = melt(down_EN_DCIS_to_reshape1, id.vars = "gene.names")
p <- plot_ly(long_EN_DCIS_down, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>% add_annotations(long_EN_DCIS_down$gene.names) %>% layout(showlegend = FALSE)  


up_DCIS_IDC_to_reshape = to_reshape[to_reshape$gene.names %in% To_GOEA_up,c(1,4,5,6)]
up_DCIS_IDC_to_reshape1 = mutate(up_DCIS_IDC_to_reshape, difference = up_DCIS_IDC_to_reshape$IDC_Mono - up_DCIS_IDC_to_reshape$EN_Mono)
up_DCIS_IDC_to_reshape1 = up_DCIS_IDC_to_reshape1[order(up_DCIS_IDC_to_reshape1$difference, decreasing = TRUE),]
up_DCIS_IDC_to_reshape1 = up_DCIS_IDC_to_reshape1[1:5,c(1:4)]
long_DCIS_IDC_up = melt(up_DCIS_IDC_to_reshape1, id.vars = "gene.names")
p <- plot_ly(long_DCIS_IDC_up, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>%  add_annotations(long_DCIS_IDC_up$gene.names) %>% layout(showlegend = FALSE)  


down_DCIS_IDC_to_reshape = to_reshape[to_reshape$gene.names %in% To_GOEA_down,c(1,4,5,6)]
down_DCIS_IDC_to_reshape1 = mutate(down_DCIS_IDC_to_reshape, difference = down_DCIS_IDC_to_reshape$IDC_Mono - down_DCIS_IDC_to_reshape$EN_Mono)
down_DCIS_IDC_to_reshape1 = down_DCIS_IDC_to_reshape1[order(down_DCIS_IDC_to_reshape1$difference, decreasing = TRUE),]
down_DCIS_IDC_to_reshape1 = down_DCIS_IDC_to_reshape1[1:5,c(1:4)]
long_DCIS_IDC_down = melt(down_DCIS_IDC_to_reshape1, id.vars = "gene.names")
p <- plot_ly(long_DCIS_IDC_down, x = ~variable, y = ~value, color = ~gene.names, text = ~gene.names) %>%
  add_lines() %>% add_annotations(long_DCIS_IDC_down$gene.names) %>% layout(showlegend = FALSE)  



