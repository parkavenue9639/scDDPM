rm(list = ls())
#install.packages("Seurat")
#devtools::install_version("Seurat", version = "4.1.4")
#BiocManager::install("Seurat")
#install.packages("Matrix")

setwd("F:/研究生文件/小论文/DDPM/数据/ALZHEIMER/AD01103/原始数据/")
#如果是表格类型的话
data <- read.delim("./AD01103_expr.txt",header=TRUE,row.names = 1)
#非0率
completed_nz <- sum(data>0)/(nrow(data)*ncol(data))
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
scell <- CreateSeuratObject(counts = data,project = "Alzheimer0205", min.cells = ceiling(0.05*ncol(data)))
data = 0
#获取基因名和 细胞名，为后面的写入文件使用
gene_names <- rownames(scell)
cell_names <- colnames(scell)
dim(scell@assays$RNA@layers$counts)

#对counts矩阵列标准化，乘以1万，然后+1取自然对数。
data <- NormalizeData(scell,normalization.method = "LogNormalize", scale.factor = 10000) # 进行标准化
norm <- data@assays$RNA@layers$data

# 创建包含基因名和细胞名的数据框
normalized_df <- as.data.frame(norm)
rownames(normalized_df) <- gene_names
colnames(normalized_df) <- cell_names


#获取前1000个高变基因
data1000 <- FindVariableFeatures(data,selection.method = "vst",nfeatures = 100)
top <- head(VariableFeatures(data1000),100) #读取选取的1000个高变基因的索引
pb <- normalized_df[top,] 
dim(pb)
pb <- t(pb)
write.csv(pb , "../预处理数据/FD100/AD01103PrePro100.csv")

#读取标签
label <- read.delim("AD01103_cell_label.txt",header=TRUE,row.names = 1)
merged_data <- merge(pb, label, by = "row.names", all.x = TRUE)
table(label)
write.csv(merged_data , "../预处理数据/FD100/AD01103PreProLabel100.csv",row.names = FALSE)
