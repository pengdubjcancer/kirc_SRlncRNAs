# 读取并探索gtf文件
gtf <- rtracklayer::import("gencode.v22.annotation.gtf")
gtf <- as.data.frame(gtf)
head(gtf,6)