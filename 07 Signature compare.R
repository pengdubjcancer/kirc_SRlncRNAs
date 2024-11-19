# 加载R包
library(org.Hs.eg.db)
library(survival)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(openxlsx)

# core functions
makeCox <- function(Features, 
                    coefs = NULL, 
                    SIGname,
                    Train_expr, 
                    Train_surv, 
                    unmatchR = 0.2,
                    statusVar = "OS", 
                    timeVar = "OS.time"){
  unmatch.Features <- setdiff(Features, colnames(Train_expr))
  
  if (all(is.na(coefs))) coefs = NULL
  if (!is.null(coefs) & is.null(names(coefs))) names(coefs) = Features
  
  if (length(unmatch.Features)/length(Features) > unmatchR){
    message(SIGname, ": Warnings, ", length(unmatch.Features), " of ", length(Features), 
            " (", round(length(unmatch.Features) / length(Features) * 100) , "%)", 
            " Features were not matched in train set. Skip this signature.")
    return(NULL)
  }else{
    if (length(unmatch.Features) > 0){
      message(SIGname, ": ", length(unmatch.Features), " of ", length(Features), 
              " (", round(length(unmatch.Features) / length(Features) * 100) , "%)", 
              " Features were not matched in train set. ")
    }
    subFeature <- intersect(Features, colnames(Train_expr))
    
    if (!is.null(coefs)){
      model = setNames(object = coefs[subFeature], nm = subFeature)
    }else{
      model <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                     data = as.data.frame(Train_expr[, subFeature]))$coefficients
    }
    return(model)
  }
  
 
}

calCindex <- function(model, name,
                      Test_expr, 
                      Test_surv, 
                      Train_expr = NULL,
                      Train_surv = NULL,
                      Train_name = NULL,
                      CohortVar = "Cohort", 
                      metaCohort = TRUE,
                      timeVar = "OS.time", 
                      statusVar = "OS"){
  
  CohortSet <- list(); CohortSurv <- list()
  CohortSet[["Test"]] <- Test_expr
  Test_surv$Sample <- rownames(Test_surv)
  CohortSurv[["Test"]] <- Test_surv[,c("Sample", CohortVar, timeVar, statusVar)]
  
  if((!is.null(Train_expr)) & (!is.null(Train_surv))){
    newSurv <- Train_surv
    
    if(!is.null(Train_name)) newSurv[[CohortVar]] <- Train_name 
    else newSurv[[CohortVar]] <- "Training"
    newSurv$Sample <- rownames(newSurv)
    
    CohortSet[["Train"]] <- Train_expr
    CohortSurv[["Train"]] <- newSurv[,c("Sample", CohortVar, timeVar, statusVar)]
  }
  
  if(metaCohort){
    CohortSet[["Meta"]] <- Test_expr
    CohortSurv[["Meta"]] <- Test_surv
    CohortSurv[["Meta"]][[CohortVar]] <- "MetaCohort"
  }
  
  Cohortlevel <- c(unique(CohortSurv[["Train"]][[CohortVar]]), unique(CohortSurv[["Test"]][[CohortVar]]))
  if(metaCohort) Cohortlevel <- c(Cohortlevel, unique(CohortSurv[["Meta"]][[CohortVar]]))
  CohortSet <- do.call(rbind, CohortSet)
  CohortSurv <- do.call(rbind, CohortSurv)
  CohortSurv[[CohortVar]] <- factor(CohortSurv[[CohortVar]], levels = Cohortlevel)

  if (is.null(model)){
    return(NULL)
  }else if(length(intersect(names(model), colnames(CohortSet))) == 0){
    message(name, ": estimate C-index using calculated RS score")
    if (length(setdiff(CohortSurv$Sample, names(model)))>0) 
      message("There are no match RS score for ", length(setdiff(CohortSurv$Sample, names(model))), " samples")
    Predict.out <- CohortSurv
    Predict.out$RS <- as.vector(model[Predict.out$Sample])
    Predict.out <- split(x = Predict.out, f = Predict.out[,CohortVar])
    f <- as.formula(paste0("Surv(", timeVar,",",statusVar,")~RS"))
    res <- do.call(rbind, lapply(Predict.out, function(data){
      s = summary(coxph(formula = f, data = data))
      c(s$concordance, n = s$n)
    }))
    res <- as.data.frame(res)
    res$Cohort <- factor(rownames(res), levels = Cohortlevel)
    colnames(res) <- c("C", "se", "n", "Cohort")
    res$method = name
    res$RS <- lapply(Predict.out, function(x) x$RS)
    return(res)
  }else{
    message(name, ": estimate C-index using Cox model fitted score")
    
    model <- model[names(model) %in% colnames(CohortSet)]
    model.mat <- as.matrix(model)
    RS <- as.matrix(CohortSet[, names(model)]) %*% model.mat
    
    Predict.out <- CohortSurv
    Predict.out$RS <- as.vector(RS)
    Predict.out <- split(x = Predict.out, f = Predict.out[,CohortVar])
    f <- as.formula(paste0("Surv(", timeVar,",",statusVar,")~RS"))
    res <- do.call(rbind, lapply(Predict.out, function(data){
      s = summary(coxph(formula = f, data = data))
      c(s$concordance, n = s$n)
    }))
    res <- as.data.frame(res)
    res$Cohort <- factor(rownames(res), levels = Cohortlevel)
    colnames(res) <- c("C", "se", "n", "Cohort")
    res$method = name
    res$RS <- lapply(Predict.out, function(x) x$RS)
    return(res)
  }
}

# other functions
standarize.fun <- function(indata, centerFlag, scaleFlag) {  
  scale(indata, center=centerFlag, scale=scaleFlag)
}

scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)
  }
  
  if (is.null(centerFlags)){
    centerFlags = F; message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags)==1){
    centerFlags = rep(centerFlags, length(data)); message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data); message("match centerFlags with cohort by order\n")
  }
  
  if (is.null(scaleFlags)){
    scaleFlags = F; message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags)==1){
    scaleFlags = rep(scaleFlags, length(data)); message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data); message("match scaleFlags with cohort by order\n")
  }
  
  centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
  outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
  # lapply(out.data, function(x) summary(apply(x, 2, var)))
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

# 加载用于模型比较的脚本
source(file.path(code.path, "compare.R"))

# 读取训练集和测试集 -----------------------------------------------------------

## Training Cohort -------------------------------------------------------------
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，表达谱需有一定变异性，以免建模过程报错）
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 训练集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -----------------------------------------------------------
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型）
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 测试集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) # 注意测试集标准化顺序与此一致
Test_set = scaleData(data = Test_expr, cohort = Test_surv$Cohort, centerFlags = T, scaleFlags = T)

# 文件为txt格式；且至少2或3列信息，具体如下：
# 必须列“Model”：以区别不同签名
# 必须列“SYMBOL”：表示签名中所含的基因
# 可选列“Coef”：代表已知签名的系数，若用户未提供该列则根据多变量Cox计算每个基因的系数
pubSIG <- read.table(file.path(data.path, "public signatures.txt"), header = T)
if (!"Coef" %in% colnames(pubSIG)) pubSIG$Coef <- NA # 若未匹配到“Coef“列，则新建“Coef“列并默认为NA
pubSIG <- split(pubSIG[, c("SYMBOL", "Coef")], pubSIG$Model)

## My Signature ----------------------------------------------------------------
mySIGname = "ours" # 本研究所定义的签名的名字，用于在图形中显示
myAlgorithm = "" # 本研究所定义的最优算法，即热图最顶部的算法名称
mySIG <- read.table(file.path(res.path, "RS_mat.txt"), header = T, check.names = F)
mySIG <- setNames(object = mySIG[[myAlgorithm]], nm = rownames(mySIG))
signatures <- pubSIG
signatures[[mySIGname]] <- mySIG
## 计算C指数 -------------------------------------------------------------------
model <- list(); cinfo <- list() # 初始化变量
log.file <- file.path(res.path, "makeCox.log") # 在Results文件夹下新建log文件
if (file.exists(log.file)) file.remove(log.file) # 此log文件用于存放在进行多变量cox分析时的警告
log.file <- file(log.file, open = "a")
sink(log.file, append = TRUE, type = "message")
for (i in names(signatures)){
  if (class(signatures[[i]]) == "data.frame"){
    model[[i]] <- makeCox(Features = signatures[[i]]$SYMBOL, # 签名的基因名
                          coefs = signatures[[i]]$Coef,      # 公共签名所提供的基因系数（如未提供也不必修改此行代码）
                          SIGname = i,                       # 当前循环的签名
                          unmatchR = 0.2,                    # 基因名不匹配率，高于该比率将被剔除；低于匹配率但大于0时会报警告，并存入log文件
                          Train_expr = Train_set,            # 用于计算cox系数的训练集表达谱
                          Train_surv = Train_surv,           # 用于计算cox系数的训练集生存信息
                          statusVar = "OS",                  # 用于构建cox模型的生存结局
                          timeVar = "OS.time")               # 用于构建cox模型的生存时间
  }else{
    model[[i]] = signatures[[i]]
  }
  
  cinfo[[i]] <- calCindex(model = model[[i]],                # 训练的cox模型，为有名字的向量
                          name = i,                          # 当前循环的签名
                          Test_expr = Test_set,              # 用于计算c指数的测试集表达谱
                          Test_surv = Test_surv,             # 用于计算c指数的测试集生存信息
                          Train_expr = Train_set,            # 用于计算c指数的训练集表达谱
                          Train_surv = Train_surv,           # 用于计算c指数的训练集生存信息
                          Train_name = "TCGA",               # 指定训练集的名称
                          #Train_expr = NULL,                # 若不需要评估训练集，则取消此行注释，并注释掉上方对应行
                          #Train_surv = NULL,                # 若不需要评估训练集，则取消此行注释，并注释掉上方对应行
                          CohortVar = "Cohort",              # 用于指定测试集所来自的队列
                          metaCohort = TRUE,                 # 指示是否将测试集合并生成MetaCohort
                          statusVar = "OS",                  # 用于计算c指数的生存结局
                          timeVar = "OS.time")               # 用于计算c指数的生存时间
  message("")
}
closeAllConnections()
cinfo <- do.call(rbind, cinfo)
write.table(cinfo[,1:5], file = file.path(res.path,"cinfo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F) # 输出不同签名在所有队列中的c指数统计量
cinfo <- split(cinfo, cinfo$Cohort)
install.packages("forestplot")
library(forestplot)
rs_forest <- read.csv('senlin2.csv',header = FALSE)
# 读入数据的时候一定要把header设置成FALSE，确保第一行不被当作列名称。
forestplot(labeltext = as.matrix(rs_forest[,1:4]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
           mean = rs_forest$V5, #设置均值
           
           lower = rs_forest$V6, #设置均值的lowlimits限
           
           upper = rs_forest$V7, #设置均值的uplimits限
           
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           xticks=c(0.5,1,2,3),
           zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           clip=c(0.49,2.51),
           boxsize = 0.4, #设置点估计的方形大小
           
           lineheight = unit(8,'mm'),#设置图形中的行距
           
           colgap = unit(2,'mm'),#设置图形中的列间距
           
           lwd.zero = 2,#设置参考线的粗细
           
           lwd.ci = 2,#设置区间估计线的粗细
           
           col=fpColors(box='red',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
           xlab="The estimates",#设置x轴标签
           
           lwd.xaxis=2,#设置X轴线的粗细
           
           lty.ci = "solid",
           
           graph.pos = 5)#设置森林图的位置，此处设置为4，则出现在第四列
