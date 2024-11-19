#突变组
library(maftools) 
library(TCGAmutations)
.libPaths(c("/home/data/t200405/R/x86_64-pc-linux-gnu-library/4.4")) 
a <- TCGAmutations::tcga_available() 
lihc.maf <- TCGAmutations::tcga_load(study = "KIRC")
lihc <- lihc.maf@data 
write.table(lihc,"lihc.txt",sep = "\t") 
a <- read.table("senstivity.txt",header = T)
total <- a$id 
# total <- as.matrix(total)
total <- subset(lihc, Tumor_Sample_Barcode_min %in% total)
Low <- subset(a, a$Metformin == "Low") 
Low <- Low$id 
Low <- subset(lihc, Tumor_Sample_Barcode_min %in% Low) 
High <- subset(a, a$Metformin == "High") 
High <- High$id 
High <- subset(lihc, Tumor_Sample_Barcode_min %in% High)
total = read.maf(maf = total,
            vc_nonSyn=names(tail(sort(table(total$Variant_Classification)))))
plotmafSummary(maf = total,
               rmOutlier = TRUE,
               addStat = "median",
               dashboard = TRUE,
               titvRaw = FALSE)
##瀑布图(Oncoplots)：
oncoplot(maf = total,
         #draw_titv = FALSE,
         #anno_height = 1,#样本注释区高度
         #legend_height = 4, #图例绘图区高度
         #drawRowBar = T, #是否显示右侧条形图
         #drawColBar = T, #是否显示顶部条形图
         top = 10) #高频突变的Top10基因
#titv函数将SNP分类为Transitions_vs_Transversions，
#并以各种方式返回汇总表的列表。汇总数据也可以显示为一个箱线图，
#显示六种不同转换的总体分布，并作为堆积条形图显示每个样本中的转换比例
laml.titv = titv(maf = total, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
#rainfall plots，展示超突变的基因组区域。
#detectChangePoints设置为TRUE，rainfall plots可以突出显示潜在变化的区域
#rainfallPlot(maf = total, detectChangePoints = TRUE, pointSize = 0.6)
#somaticInteractions函数使用配对Fisher 's精确检验来分析突变基因之间
#的的co-occurring 或者exclusiveness
#exclusive/co-occurance event analysis on top 10 mutated genes.
Interact <- somaticInteractions(maf = total, top = 15, pvalue = c(0.05, 0.1))
#提取P值结果
Interact$gene_sets
#mafComapre参数比较两个不同队列的差异突变基因，检验方式为fisher检验。

High = read.maf(maf = High,
                 vc_nonSyn=names(tail(sort(table(High$Variant_Classification)))))
plotmafSummary(maf = High,
               rmOutlier = TRUE,
               addStat = "median",
               dashboard = TRUE,
               titvRaw = FALSE)
##瀑布图(Oncoplots)：
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'In_Frame_Ins',
  'Translation_Start_Site',
  'Nonstop_Mutation',
  'Frame_Shift_Ins',
  'Frame_Shift_Del',
  'Splice_Site',
  'In_Frame_Del',
  'Missense_Mutation'
)
oncoplot(maf = High, colors = vc_cols,
         #draw_titv = FALSE,
         #anno_height = 1,#样本注释区高度
         #legend_height = 4, #图例绘图区高度
         #drawRowBar = T, #是否显示右侧条形图
         #drawColBar = T, #是否显示顶部条形图
         top = 10) #高频突变的Top10基因

Low = read.maf(maf = Low,
                 vc_nonSyn=names(tail(sort(table(Low$Variant_Classification)))))
plotmafSummary(maf = Low,
               rmOutlier = TRUE,
               addStat = "median",
               dashboard = TRUE,
               titvRaw = FALSE)
##瀑布图(Oncoplots)：
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'In_Frame_Ins',
  'Translation_Start_Site',
  'Nonstop_Mutation',
  'Frame_Shift_Ins',
  'Frame_Shift_Del',
  'Splice_Site',
  'In_Frame_Del',
  'Missense_Mutation'
)
oncoplot(maf = Low, colors = vc_cols,
         #draw_titv = FALSE,
         #anno_height = 1,#样本注释区高度
         #legend_height = 4, #图例绘图区高度
         #drawRowBar = T, #是否显示右侧条形图
         #drawColBar = T, #是否显示顶部条形图
         top = 10) #高频突变的Top10基因

# 使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1=High, m2=Low, m1Name="High", m2Name="Low", minMut=10)
result <- fvsm$results
forestPlot(mafCompareRes=fvsm, pVal=0.05, color=c("maroon", "royalblue"), geneFontSize=0.8)

tmb_low = tmb(maf = Low)   #默认以log10转化的TMB绘图
tmb_low = tmb_low$total_perMB_log
tmb_high = tmb(maf = High)   #默认以log10转化的TMB绘图
tmb_high = tmb_high$total_perMB_log
t_test_result <- t.test(tmb_high, tmb_low)
data <- data.frame(
  group = c(rep("tmb_high", length(tmb_high)), rep("tmb_low", length(tmb_low))),
  value = c(tmb_high, tmb_low)
)
ggplot(data, aes(x = group, y = value)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Group Comparison")
ggboxplot(data,x = "group",y = "value",color = "group",add = "jitter") + 
  stat_compare_means(method = "t.test")