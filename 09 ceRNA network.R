library(multiMiR)
db.ver = multimir_dbInfoVersions()
db.ver[,1:3]
db.count

dePC <- read.table("dePC.txt",sep = "\t",header = T,row.names = 1)
dePC <- as.matrix(dePC)
example4 <- get_multimir(org     = 'hsa',
                         target  = dePC,
                         table   = 'validated',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)
table(gene2mir@data[["support_type"]])
table(gene2mir@data[["experiment"]])
table(gene2mir@data[["mature_mirna_id"]])
table(gene2mir@data$database)
example1_result <- example4@data 
example1_sum <- example4@summary                   # 提取summary数据
head(example1_sum)
apply(example4@summary[, 6:10], 2, sum)
example1_Lucifer <- example1_result[grep("Luciferase", example1_result[, "experiment"]), ]
head(example1_Lucifer)

# 提取3个数据检索得到的target_symbol，绘制Veen图
library(VennDiagram)            
library(tidyverse)
db <- unique(example1_result$database)
row1 <- example1_result %>% filter(database == db[1])
row2 <- example1_result %>% filter(database == db[2])
row3 <- example1_result %>% filter(database == db[3])
venn.plot <-venn.diagram(
  x = list(target1=unique(row1$target_symbol),
           target2=unique(row2$target_symbol),
           target3=unique(row3$target_symbol)),
  filename ="Venn1.tif",
  lty ="dotted",
  lwd =0.5,
  col ="black",
  fill =c("dodgerblue", "goldenrod1", "darkorange1"),
  alpha =0.60,
  cat.col =c("dodgerblue", "goldenrod1", "darkorange1"),
  cat.cex =1,
  cat.fontface="bold",
  margin =0.05,
  cex =1)
# 获取交集
venn_list <- list(target1=unique(row1$target_symbol),
                  target2=unique(row2$target_symbol),
                  target3=unique(row3$target_symbol))
for (i in 1:length(venn_list)) {
  if(i == 1){interGenes <- venn_list[[1]]}
  else{interGenes <- intersect(interGenes, venn_list[[i]])}
}
interGenes

write.csv(example1_Lucifer,"example1_Lucifer.csv")



library(igraph)
library(dplyr)
library(magrittr)
network_data <- read.csv("cerna.csv", header = TRUE)
# 创建空的网络对象
g <- graph.empty(n =length(c(unique(network_data$mirna),unique(network_data$lncrna),unique(network_data$mrna))), directed = TRUE)

# 添加节点
g <- g %>%
  set_vertex_attr("name", value = c(unique(network_data$lncrna), unique(network_data$mirna), unique(network_data$mrna))) %>%
  set_vertex_attr("type", value = c(rep("lncrna", length(unique(network_data$lncrna))), 
                                    rep("mirna", length(unique(network_data$mirna))), 
                                    rep("mrna", length(unique(network_data$mrna)))))
g <- set_vertex_attr(g,"color", value = ifelse(V(g)$type == "lncrna", "#fb8072", ifelse(V(g)$type == "mirna", "yellow3", "#80b1d3")))

# 添加边与边长
afedge <- c()
aflength <- c()
for(i in 1:nrow(network_data)) {
  lncrna_node <- which(V(g)$name == network_data[i,1])
  mirna_node <- which(V(g)$name == network_data[i,2])
  mrna_node <- which(V(g)$name == network_data[i,3])
  aflength <- c(aflength,20,10)
  afedge <- c(afedge,lncrna_node,mirna_node,mirna_node,mrna_node)
  
}
g <- g %>% add_edges(afedge) %>% set_edge_attr("edge.length", value = aflength)

# 添加节点大小
lncrna.size=as.vector(scale(as.vector(table(network_data$lncrna)),center = F))+15
mirna.size=as.vector(scale(as.vector(table(network_data$mirna)),center = F))+6
mrna.size=as.vector(scale(as.vector(table(network_data$mrna)),center = F))+10
V(g)$size=c(lncrna.size,mirna.size,mrna.size)
# 使用Graphopt算进行布局，保存为ceRNA.net.pdf文件
pdf(file="ceRNA.net.pdf",height=10,width=10)
plot(g, 
     layout=layout.graphopt(g),  
     vertex.label=V(g)$name,
     vertex.label.family="sans",
     vertex.label.cex=ifelse(V(g)$type == "lncRNA", 0.8, ifelse(V(g)$type == "mirna", 0.8, 0.8)),
     vertex.size=V(g)$size, 
     vertex.color=V(g)$color,
     vertex.label.color="black", 
     edge.arrow.size=0.5, 
     edge.width=0.5
)
dev.off()
