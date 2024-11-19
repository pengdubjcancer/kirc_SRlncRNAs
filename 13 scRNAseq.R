单细胞测序
sc.data<-Read10X_h5("LIHC_GSE159115_expression.h5")
data <- CreateSeuratObject(counts = sc.data, project = "GSE159115") 
rm(sc.data)
metadata = read.delim("LIHC_GSE159115_CellMetainfo_table.tsv",row.names = 1)  #单细胞其他信息读取 txt格式
data@meta.data <- metadata
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)  #top2000差异基因
top10 <- head(VariableFeatures(data), 10)   #top10差异基因
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
##图3：展示top2000差异基因中Top10高异质性基因
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:20)
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:8)
.libPaths(c("/refdir/Rlib"))
data <- RunUMAP(data, dims = 1:8) 
Idents(data) <- data@meta.data$Cluster 
Idents(data) <- data@meta.data$Celltype..malignancy.
Idents(data) <- data@meta.data$Celltype..major.lineage.
DimPlot(data, reduction = "umap", label = TRUE) 
library(ggplot2)
# list_genes=split(topN$gene, topN$cluster)
list_genes=list(Immu=c("PTPRC"),
                Tcell=c("CD3D", "CD3E"),
                CD8T=c("CD8A", "CD8B"),
                Plasma=c("JCHAIN","IGKC","SLAMF7"),
                MonoMacro=c("MMP9","CST3","CD68","FABP4","CD14","CSF1R"), 
                Endo=c("CHGB","PECAM1","VWF"),
                Malignant=c("CDH1","MYC","CD24","CA9","NDUFA4L2"),
                Epith=c("EPCAM"),
                Pericyte=c("RGS5"))
DotPlot(data, features=list_genes, 
        cols = c('#330066','#336699','#66CC66','#FFCC33'), 
        cluster.idents = T)+
  RotatedAxis()+
  theme(panel.border = element_rect(color="black"), #面板边框
        panel.spacing = unit(1, "mm"), #面板间距
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 8),
        #axis.text = element_text(size = 14),
        # 分面标题
        strip.background = element_rect(color="red"),
        strip.text = element_text(size = 8, margin=margin(b=3, unit="mm")),
        strip.placement = 'outlet', #
        #axis.line = element_blank(), # 坐标轴线 
        axis.text.x = element_text(angle = 90, hjust = 0.5,vjust=0.5)
  )+labs(x="", y="") + 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
VlnPlot(data, features = c('EMX2OS')) 
VlnPlot(data, features = c('LINC00944')) 

sce_res <- data
for (i in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3,0.4, 0.5,0.8,1)){
  sce_res <- FindClusters(sce_res,resolution = i)
}
clustree(sce_res,prefix = 'RNA_snn_res.') 

#拟时序分析
BiocManager::install("monocle",force = TRUE)
dir.create("pseudotime")
table(Idents(data))
table(data@meta.data$seurat_clusters)
table(data@meta.data$orig.ident)
levels(Idents(data))
data = data[, Idents(data) %in% c("Epithelial","Malignant","Fibroblast")]
levels(Idents(data))
head(mycds@phenoData@data)
mycds <- readRDS("cds.RDS")  #scRNAsub是上一节保存的T细胞子集seurat对象

sample_ann <- data@meta.data
sample_ann$celltype = Idents(data)

gene_ann <- data.frame(gene_short_name = rownames(data@assays$RNA),
                       row.names = rownames(data@assays$RNA))

pd <- new('AnnotatedDataFrame', data = sample_ann)
fd <- new('AnnotatedDataFrame', data = gene_ann)
subdata <- as.data.frame(data@assays$RNA@counts)
mycds <- newCellDataSet(as.matrix(subdata),
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size(),
                        lowerDetectionLimit = 1)
mycds <- detectGenes(mycds, min_expr = 1) 
mycds <- mycds[fData(mycds)$num_cells_expressed > 10,]

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds) 

disp_table <- dispersionTable(mycds)
pData(mycds)$BZW1 = log2(exprs(mycds)["BZW1",]+1)
plot_cell_trajectory(mycds,color_by = "BZW1")+ scale_color_gsea()

pData(mycds)$BZW2 = log2(exprs(mycds)["BZW2",]+1)
plot_cell_trajectory(mycds,color_by = "BZW2")+ scale_color_gsea()

unsup_clustering_genes <- subset(disp_table,mean_expression >= 0.1)

mycds <- setOrderingFilter(mycds,unsup_clustering_genes$gene_id)

mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')

mycds <- orderCells(mycds) 
head(mycds@phenoData@data)
plot_cell_trajectory(mycds)##以state进行着色
plot_cell_trajectory(mycds , color_by = "State")
plot_cell_trajectory(mycds , color_by = "Pseudotime",ce1l_size = 0.75)
##根据拟时间值着色#
#绘制cluster的分面图
plot_cell_trajectory(mycds , color_by = "seurat_clusters" ,ce11_size = 0.75)
plot_cell_trajectory(mycds , color_by = "orig.ident" ,ce11_size = 0.75)
plot_cell_trajectory(mycds , color_by = "celltype" ,ce11_size = 0.75)
pData(mycds)$Clusters = pData(mycds)$celltype
diff_test_res <- differentialGeneTest(mycds,fullModelFormulaStr = "~Cluster")

sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes <- sig_genes[order(sig_genes $ pval), ]
head(sig_genes[,c("gene_short_name","pval","qval")])
cg = as.character(head(sig_genes$gene_short_time))
plot_genes_jitter(mycds[BZW2,], grouping = "Cluster", color_by = "Cluster", nrow =3 ,ncol=NULL)
cg2 = as.character(tail(sig_genes$gene_short_time))
plot_gene_jitter(cds[cg2,], grouping = "Cluster", color_by = "Cluster", nrow =3 ,ncol=NULL)

#infercnv and infercnvpy
#R code
library(AnnoProbe)
library(Seurat)
library(infercnv) 
cells.use=colnames(data)
# cells.use=colnames(data)[which(data$CellType %in% c("Epithelial_cells","B_cells","Oligo&Astrocytes"))]
dat=as.data.frame(GetAssayData(subset(data, cells=cells.use)))
groupinfo=data.frame(v1=colnames(dat),v2=data@active.ident[cells.use]) 
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match(geneInfor[,1], rownames(dat) ),]
dim(dat) 

head(groupinfo)
write.table(groupinfo,file = 'groupFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
## 注意这里的exp和R版不一样，做了转置
write.table(t(dat),file = "expFile_t.txt",sep = '\t',quote = F)
write.table(geneInfor,file = "geneFile.txt",sep = '\t',quote = F,col.names = F,row.names = F)
#python code
import scanpy as sc
import infercnvpy as cnv
import pandas as pd
adata = sc.read_text('expFile_t.txt')
anno = pd.read_table('groupFiles.txt',index_col='Id') ## 这里要对文件做相应更改，添加一个列名
adata.obs['cell_groups'] = anno['cells']
gene = pd.read_table('geneFile.txt',index_col='Gene') ## 同样做相应更改，添加列名
adata.var=gene 
# 后面教程就差不多了
cnv.tl.infercnv(
   adata,
   reference_key="cell_groups",
   reference_cat=["CD8T",
       "Endothelial",
       "Erythroblasts",
       "MonoMacro",
       "Pericyte",
       "Plasma",
  ],
   window_size=250,
) 
sc.settings.set_figure_params(dpi=800, facecolor="white")
cnv.pl.chromosome_heatmap(adata, groupby="cell_groups") 
#leiden聚类
cnv.pl.chromosome_heatmap(adata, groupby="leiden", save="heatmap_leiden.png")
##最后结果在这里
adata.obsm["X_cnv"] 

#Cytotrace 
#install.packages() 安装这个几个包
library(CytoTRACE)
exp1 <- as.matrix(data@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 8)
#可视化 CytoTRACE 结果
#plotCytoTRACE(results, phenotype = exp1)
# saveRDS(results,"results.rds") 
phenot <- data$seurat_clusters
phenot <- as.character(phenot)
names(phenot) <- rownames(data@meta.data)
emb <- data@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb)
plotCytoGenes(results, numOfGenes = 30)
#cytotrace联合monocle
cds <- readRDS("cds.rds")
monocle_meta <- data.frame(t(cds@reducedDimS), 
                           cds$Pseudotime, 
                           cds$State, 
                           cds$seurat_clusters)
colnames(monocle_meta) <- c("C1", "C2", "Pseudotime", "State", "Clusters")
phenot1 <- cds$seurat_clusters
phenot1 <- as.character(phenot1)
names(phenot1) <- rownames(monocle_meta)
emb_monocle <- monocle_meta[,1:2]
plotCytoTRACE(results, phenotype = phenot, emb = emb_monocle)
