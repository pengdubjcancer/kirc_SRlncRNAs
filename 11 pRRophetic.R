#引用包
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
pFilter=0.05               #pvalue的过滤条件
expFile="HiSeqV2.txt"         #表达数据文件
riskFile="SAMPLE_riskscore.txt"     #风险文件
allDrugs=c("AICAR", "AKT.inhibitor.VIII", "ATRA", "Axitinib", "Bexarotene", "Bicalutamide","Bleomycin",
           "Bortezomib", "Bosutinib","Camptothecin", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin",  "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide","Gefitinib", "Gemcitabine", "GSK269962A", "Imatinib", "JNK.Inhibitor.VIII", "Lapatinib", "Lenalidomide", "Metformin", "Methotrexate",  "Midostaurin", "Mitomycin.C", "Nilotinib", "Obatoclax.Mesylate", "Paclitaxel", "Parthenolide", "Pazopanib", "Pyrimethamine",  "Rapamycin", "Roscovitine", "Salubrinal", "Shikonin",  "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "Vinblastine", "Vinorelbine", "Vorinostat" )

#读取表达输入文件,并对数据进行处理
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)



for(drug in allDrugs){
#预测药物敏感性
senstivity=pRRopheticPredict(data, drug, selection=1)
senstivity=senstivity[senstivity!="NaN"]
#senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#风险文件和药物敏感性合并
sameSample=intersect(row.names(risk), names(senstivity))
risk=risk[sameSample, "Group",drop=F]
senstivity=senstivity[sameSample]
rt=cbind(risk, senstivity)

#设置比较组
rt$risk=factor(rt$Group, levels=c("low", "high"))
type=levels(factor(rt[,"Group"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
			      xlab="Risk",
			      ylab=paste0(drug, " senstivity (IC50)"),
			      legend.title="Risk",
			      palette=c("green", "red")
			     )+ 
	stat_compare_means(comparisons=my_comparisons)

pdf(file=paste0(drug, ".pdf"), width=5, height=4.5)
print(boxplot)
write.table((boxplot[["plot_env"]][["data"]]),file=paste0(drug,".txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
dev.off()
}
drug
a<-boxplot[["plot_env"]][["data"]]
print(a)
###