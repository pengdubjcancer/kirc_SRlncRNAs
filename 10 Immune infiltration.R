#免疫浸润 
clin <- read.csv("group.csv",sep = ",",header =T,row.names = 1)
exp <- read.table("TIMER.txt",header = T,row.names = 1,sep = "\t")
exp <- read.table("XCELL.txt",header = T,row.names = 1,sep = "\t") 
exp <- read.table("EPIC.txt",header = T,row.names = 1,sep = "\t")
exp <- read.table("QUANTISEQ.txt",header = T,row.names = 1,sep = "\t")
exp <- read.table("CIBERSORT.txt",header = T,row.names = 1,sep = "\t")
exp <- read.table("MCPCOUNTER.txt",header = T,row.names = 1,sep = "\t")
exp <- read.table("CIBERSORT-ABS.txt",header = T,row.names = 1,sep = "\t")
exp <- t(exp)
#计算显著性
table(clin$Group)
library(limma)
group_list <- factor(rep(c('Low','High'),c(186,185)))
group_list
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp)
contrast.matrix <- makeContrasts('Low-High',levels = design)
fit <- lmFit(exp,design = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
alldiff=topTable(fit2,coef = 1,n = Inf)
alldiff$lab = as.factor(ifelse(alldiff$P.Value>=0.05,"",
                               ifelse(alldiff$P.Value>=0.01&alldiff$P.Value<0.05,"*",
                                      ifelse(alldiff$P.Value>=0.001&alldiff$P.Value<0.01,"**", "***"
                                      ))))
alldiff$new <- paste(rownames(alldiff),alldiff$lab)
alldiff <- alldiff[rownames(exp),]
rownames(exp) <- alldiff$new
pre_heatdata <- t(scale(t(exp)))
pre_heatdata[pre_heatdata> 1] <- 1###这里的数值自己调，哪个好看就怎么来
pre_heatdata[pre_heatdata< -1] <- -1###这里的数值自己调，哪个好看就怎么来
annColors <- list()
annColors[['group']] <- c('low'='steelblue','high'='orange')
library(pheatmap)
pheatmap(pre_heatdata,
         color = colorRampPalette(c("#5bc0eb",'black',"#ECE700"))(1000),
         annotation_col = clin,
         annotation_colors = annColors,
         treeheight_row = 50,
         gaps_col = 186,
         #gaps_row = 15,
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F)