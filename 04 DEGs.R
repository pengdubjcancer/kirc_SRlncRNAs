getwd()
setwd("/home/data/t200405/DEGs")
data <-read.table('tcga差异表达矩阵.txt',header=T,sep = "\t")
data = as.data.frame(avereps(data[,-1],ID = data$Tag) )
group <- c(rep('normal',178),rep('tumor',167))
group <- factor(group)
table(group_list)
# 构建DESeq2中的对象
colData <- data.frame(row.names = colnames(data),
                      group = group_list)
colData$group <- factor(colData$group, levels = c("normal", "tumor"))
head(colData)
data_int <- 2^(data) - 1
data_int <- apply(data_int, 2, as.integer)
rownames(data_int) <- rownames(data)
dds <- DESeqDataSetFromMatrix(
  countData = data_int,
  colData = colData,
  design = ~ group)
# 进行差异表达分析
dds <- DESeq(dds)
# 查看结果的名称
resultsNames(dds)
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))
resOrdered <- res[order(res$padj), ]
DEG_DESeq2 <- as.data.frame(resOrdered)
write.table(DEG_DESeq2,file = "DEG_DESeq2.txt")
#edgeR
BiocManager::install("edgeR")
d = DGEList(counts = data_int, group = group)
# 根据每个基因在每个样本中的 CPM（Counts Per Million）值去除低表达基因
keep <- rowSums(cpm(d) > 1) >= 2
# 或者自动过滤，去除低表达基因
# keep <- filterByExpr(d)
table(keep)
# 从 DGEList 对象中筛选出符合条件的基因
d <- d[keep, , keep.lib.sizes = FALSE]
# 更新样本的库大小信息
d$samples$lib.size <- colSums(d$counts)
# 归一化，TMM 方法
d <- calcNormFactors(d)
# 将归一化后的数据赋值给 dge 变量
dge = d
# 创建设计矩阵，用于指定差异分析模型
design <- model.matrix(~0 + factor(group))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group))
# 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
# 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
fit <- glmFit(dge, design)
# 使用 LRT（Likelihood Ratio Test）计算差异表达
# 注意这里的 contrast 和 DESeq2 不一样，这里我们只需要输入 c(-1, 1) 即可
# -1 对应 normal，1 对应 tumor
lrt <- glmLRT(fit, contrast = c(-1, 1))
# 从 LRT 计算结果中获取前 nrow(dge) 个顶部差异表达基因
nrDEG <- topTags(lrt, n = nrow(dge))
# 将差异表达基因结果转换为数据框形式
DEG_edgeR <- as.data.frame(nrDEG)
write.table(DEG_edgeR,file = "DEG_edgeR.txt")

#limma
library(limma)
exprSet <- data_int

# 创建设计矩阵，指定组别信息
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exprSet)
# 创建 DGEList 对象
dge <- DGEList(counts = exprSet, group = group)
# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)
# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE, normalize = "quantile")
# 如果是芯片数据、TPM数据或已标准化的数据，不需要再进行标准化，可直接从这里开始进行差异分析
# 使用线性模型进行拟合
fit <- lmFit(v, design)
# 和上面两个包一样，需要说明是谁比谁
con <- paste(rev(levels(group)), collapse = "-")
con
# 创建对比矩阵
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)
head(DEG_limma_voom)
write.table(DEG_limma_voom,file = "DEG_limma_voom.txt")
