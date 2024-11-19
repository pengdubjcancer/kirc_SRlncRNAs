install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
library(synapser) 
# 登录
synLogin(email = "dupeng1587@proton.me", password = "xxxxxxxx")
# Welcome,  !NULL
# 获取数据
synRNA <- synGet( "syn2701943", downloadLocation = "~/Downloads/PCBC" )

setwd("/home/data/t200405/synapser")
exp <- read.table("rnaseq_norm.tsv", header = TRUE)
exp <- distinct(exp,Gene,.keep_all = T)
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp[1:3,1:3]
#ID转换
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

y <- read.csv("syn3156503.csv",header = T,sep = ",")
y <- t(y)
colnames(y) <- y[1,]
y <- y[-1,]
head(y)
# 对数据进行均值中心化
X <- exp
m <- apply(X, 1, mean)
X <- X-m
# 将样本分为干细胞组和非干细胞组
sc <- which(y == "SC")
X.sc <- X[,sc]
X.or <- X[,-sc]
X <- t(X)
X.sc <- t(X.sc)
X.sc <- as.matrix(X.sc)
model.RNA <- gelnet(X.sc, NULL, 0, 1)
#得到每个基因的权重
head(model.RNA$w)
#TSPAN6          TNMD          DPM1         SCYL3      C1orf112         FUCA2 
#6.828043e-05  3.261994e-03  3.170824e-04 -1.860469e-05  1.114508e-03  2.916178e-04 
save(X, y, model.RNA, file = "model.rda")
library(tidyverse)
library(openxlsx)
#1）读取表达量数据，，或者给绝对路径
expr <- read.table("TCGA-KIRC.htseq_fpkm.tsv",header=T,sep="\t",row.names=1)
#2）读取probeMap文件，转换Ensembl_ID
#ID和Gene symbol对应列表
geneann<-read.table("gencode.v22.annotation.gene.probeMap",header=T,sep="\t",row.names=1)
#二者ID进行匹配，并添加一列gsym
expr$gsym <- geneann[rownames(expr),]$gene
#去除重复的Gene name
expr<-distinct(expr,gsym,.keep_all=T)
#将行名改为Gene name
row.names(expr)<-expr$gsym
#将添加的gsym这一列删除
expr<-subset(expr,select= -gsym)
exp <- expr
rm(expr)
#提取交叠基因的表达谱及权重
common <- intersect(names(model.RNA$w), rownames(exp))
X <- exp[common, ]
w <- model.RNA$w[common]
#对于 RNA 表达数据，使用 spearman 计算权重与表达值之间的相关性来衡量样本的干性指数，并进行标准化使其落在 [0,1] 之间
score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
score <- score - min(score)
score <- score / max(score)
head(score)
write.table(score,"score.txt")