getwd()
set.seed(1235)
setwd("/home/data/t200405/synapser")
options(stringsAsFactors = FALSE)
# 加载数据集
mRNA <- read.table("mRNA.txt",header = T); 
# lncRNA <- read.table("lncRNA.txt",header = T);
# 查看数据
dim(mRNA)
names(mRNA)

#### 数据清洗
# 提取表达矩阵并行列转置，将行名变为样本名
datExpr0 <- as.data.frame(t(mRNA[, -c(1:1)]))
names(datExpr0) <- mRNA$id
rownames(datExpr0) <- names(mRNA)[-c(1:1)]
# 判断矩阵的样本是否都合格？
gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

traitData = read.table("score.txt",header = T);
dim(traitData)  #每行是一个样本，每列是一种信息
names(traitData)
femaleSamples = rownames(datExpr0);
traitRows = match(femaleSamples, traitData$id);
datTraits = traitData[traitRows, -1];
datTraits <- as.matrix(datTraits)
rownames(datTraits) = traitData[traitRows, 1];
colnames(datTraits) = "RNAss"
collectGarbage()

sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datExpr0, datTraits, file = "FemaleLiver-01-dataInput.RData")

lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red")  #查看位于0.9以上的点，可以改变高度值

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$powerEstimate
#结果为3
net = blockwiseModules(datExpr0, power = 3,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
table(net$colors)

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "FemaleLiver-02-networkConstruction-auto.RData")

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
colnames(datTraits)
# 用热图的形式展示相关系数
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#colors = greenWhiteRed(50)不适用于红绿色盲患者，建议用 blueWhiteRed代替.
#基因与表型数据的关系、重要模块：基因显著性和模块成员

RNAss = as.data.frame(datTraits);
names(RNAss) = "RNAss";
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, RNAss, use = "p"));#和体重性状的关联
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(RNAss), sep="");
names(GSPvalue) = paste("p.GS.", names(RNAss), sep="");
#运行以下代码可视化GS和MM
#module = c("tan","black","lightcyan","orange")
rownames(black)
module = "tan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
table(moduleGenes)
tan_module <-as.data.frame(dimnames(data.frame(datExpr0))[[2]][moduleGenes]) 
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneTraitSignificance[moduleGenes, 1])
tan <-as.data.frame(cbind(MM,GS)) #包含了MM和GS的数据，可以保留一下
rownames(tan)=tan_module[,1]
#green_hub <-abs(c$MM)>0.8&abs(c$GS)>0.2 #筛选hub基因 
write.csv(tan, "hubgene_MMGS_tan.csv") 
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for RNAss",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "cyan")
abline(h=0.3,v=0.5,col="red",lwd=1.5) 

names(datExpr0)#会返回所有在分析中的基因ID
tan<-names(datExpr0)[moduleColors=="tan"]#返回属于棕色模块的基因ID
tan<-write.table(tan,"tan.txt")
modOrder = order(-abs(cor(MEs, BZW2, use = "p")));
##添加模块成员的信息：
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.RNAss));  # 排序
geneInfo = geneInfo0[geneOrder, ]
输出为CSV格式，可用fix(geneInfo)在R中查看：
write.csv(geneInfo, file = "geneInfo.csv")