rm(list = ls())

setwd("F:/研究生文件/小论文/DDPM/数据/ALZHEIMER/AD00204/预处理数据/FD100/")
#原数据
data <- read.csv("./AD00204NormalLabel.csv",header=TRUE,row.names = 1)
#生成数据
gen_data <- read.csv("../../缩放后数据/Final2.0w.csv",header=TRUE,row.names = 1)
#非0率
library(Seurat)
library(SeuratObject)
library(magrittr)
library(dplyr)
sum(data == 0)
data[1:5,1:5]
dim(data)

#在文件中给细胞名赋值后再转到这执行下面程序
#CreateSeuratObject要求矩阵行为基因 ， 列为细胞
#去除不符合条件的 细胞 和 基因
data = data[,1:100]
data = t(data)
scell <- CreateSeuratObject(counts = data,project = "Alzheimer0205")
#获取基因名和 细胞名，为后面的写入文件使用
dim(scell@assays$RNA@layers$counts)

#获取前1000个高变基因
data1000 <- FindVariableFeatures(scell,selection.method = "vst",nfeatures = 100)
top <- head(VariableFeatures(data1000),100) #读取选取的1000个高变基因的索引
write.csv(top , "../NormalTop100.csv",row.names = FALSE)





gen_data[1:5,1:5]
dim(gen_data)

#在文件中给细胞名赋值后再转到这执行下面程序
#CreateSeuratObject要求矩阵行为基因 ， 列为细胞
#去除不符合条件的 细胞 和 基因
gen_data = gen_data[,1:100]
gen_data = t(gen_data)
gen_scell <- CreateSeuratObject(counts = gen_data,project = "Alzheimer0205")
#获取基因名和 细胞名，为后面的写入文件使用
dim(gen_scell@assays$RNA@layers$counts)


#获取前1000个高变基因
gen_data100 <- FindVariableFeatures(gen_scell,selection.method = "vst",nfeatures = 100)
top <- head(VariableFeatures(gen_data100),100) #读取选取的1000个高变基因的索引
write.csv(top , "../GenTop100.csv",row.names = FALSE)