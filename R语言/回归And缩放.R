#均值、方差回归
setwd(dir = "F:/研究生文件/小论文/DDPM/数据/ALZHEIMER/AD00202/预处理数据/")
getwd()
rm(list = ls())
" 回归操作 "
#此处要求一行是一个细胞 ， 列是基因 
A_org = read.csv("FD1000/AD00203NormalLabel.csv",header = T,row.names = 1)  #这个文件里面是有基因名和细胞名的（第一列和第二列）
A_gen = read.csv("F:/PythonCoding/Jupyter/DDPM/data/output/DiffusionCell/AD00202-1000/初试生成/AD00202-1000outputws2.0.csv",header = T,row.names = 1)
#该类中细胞的标签个数
labelnum = 6
#缩放后存放的位置
det_file = "../缩放后数据/Final2.0w(2000).csv"

A_org[1:4,1:4]
A_gen[1:4,1:4]

#将所有的标签从浮点型转变为int类型
A_org[, ncol(A_org)] <- as.integer(A_org[, ncol(A_org)]) - 1
A_gen[, ncol(A_gen)] <- as.integer(A_gen[, ncol(A_gen)]) - 1


quantile.prob = 0.00001
cat(sprintf("精确到基因表达的%f位\n", quantile.prob))
# 对每一个类进行缩放
file.remove(det_file)

for (curlabel in 1:labelnum-1) {
  print(curlabel)
  #筛选出当前标签的原本数据 以及 生成数据
  A_org_curlabel = A_org[A_org[, ncol(A_org)] == curlabel, ]
  A_gen_curlabel = A_gen[A_gen[, ncol(A_gen)] == curlabel, ]
  originally_nonzero <- A_org_curlabel >0
  #设置最低值关卡，设置每一个基因的最低值进行过滤
  A_gen_curlabel_mins <- abs(apply(A_gen_curlabel,2,FUN=function(x) quantile(x,quantile.prob)))
  cat("清除最低阈值以下的值\n")
  A_gen_curlabel_cor <- replace(A_gen_curlabel, A_gen_curlabel < A_gen_curlabel_mins[col(A_gen_curlabel)-1], 0)
  
  #计算非零元素的标准差
  sd_nonzero <- function(x) sd(x[!x == 0])
  #对原数据 和 生成数据 的每一列进行标准差sd()操作
  sigma_1 <- apply(A_gen_curlabel_cor, 2, sd_nonzero)
  sigma_2 <- apply(A_org_curlabel, 2, sd_nonzero)
  
  #每一列所有值之和 除以 每一列的非零的个数 （均数）
  mu_1 <- colSums(A_gen_curlabel_cor)/colSums(!!A_gen_curlabel_cor)
  mu_2 <- colSums(A_org_curlabel)/colSums(!!A_org_curlabel)
  
  #检测值是否存在
  toscale <- !is.na(sigma_1) & !is.na(sigma_2) & !(sigma_1 == 0 & sigma_2 == 0) & !(sigma_1 == 0)
  
  cat(sprintf("缩放 除了 %d 列的所有列\n", sum(!toscale)))
  
  sigma_1_2 <- sigma_2/sigma_1
  toadd  <- -1*mu_1*sigma_2/sigma_1 + mu_2
  #记录要缩放的列
  A_gen_curlabel_temp <- A_gen_curlabel_cor[,toscale]
  A_gen_curlabel_temp <- sweep(A_gen_curlabel_temp,2, sigma_1_2[toscale],FUN = "*")
  A_gen_curlabel_temp <- sweep(A_gen_curlabel_temp,2, toadd[toscale],FUN = "+")
  
  A_gen_curlabel_cor_sc <- A_gen_curlabel_cor
  A_gen_curlabel_cor_sc[,toscale] <- A_gen_curlabel_temp
  A_gen_curlabel_cor_sc[A_gen_curlabel_cor==0] = 0
  
  lt0 <- A_gen_curlabel_cor_sc  <0
  A_gen_curlabel_cor_sc[lt0] <- 0 
  
  write.table(A_gen_curlabel_cor_sc, det_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}