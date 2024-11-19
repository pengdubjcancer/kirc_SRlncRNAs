#清空环境变量
rm(list=ls())   
options(stringsAsFactors = F)
library(psych)
mRNA<-read.table("TCGA_stemness_mRNA.txt",header = T,sep = "\t", check.names=FALSE)
lncRNA<-read.table("TCGA_stemness_lncRNA.txt",header = T,sep = "\t", check.names=FALSE)

res<-corr.test(mRNA,lncRNA,use="pairwise",method="pearson",adjust="holm",alpha=0.05)
warnings()
p<-read.table("mrna_lncrna_p.txt",header = T,sep = "\t")
p<-res$p
r<-res$r
write.table(p,file="mrna_lncrna_p.txt",row.names = TRUE,col.names = TRUE,sep = "\t")
write.table(r,file="mrna_lncrna_r.txt",row.names = TRUE,col.names = TRUE,sep = "\t")
library(pheatmap)
pdat2<-matrix(ifelse(abs(r)<=0.6,0,r),nrow(r))
pdat2<-matrix(ifelse(p>=0.00001,0,pdat2),nrow(r))
bk<-seq(-1,1,by=0.1)

colnames(pdat2)<-colnames(r)
rownames(pdat2)<-rownames(r)

mycol<-c(colorRampPalette(c('blue','white'))(9),
         colorRampPalette(c('white','white'))(5),
          colorRampPalette(c("white","red"))(8))
pheatmap(pdat2,cellwidth = 10,cellheight = 30,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = mycol,
         legend_breaks = seq(-1,1,by=0.2),
         breaks = bk)

