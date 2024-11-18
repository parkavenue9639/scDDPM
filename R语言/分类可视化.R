#该文件想要做可视化，在视觉上查看生成数据是否和原数据重合
#install.packages("BiocManager")
#BiocManager::install("scater")
#BiocManager::install("SC3")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("umap")
library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library(scater)
library(Rtsne)
library(umap)
getwd()
rm(list =ls())
setwd("F:/研究生文件/小论文/DDPM/数据/ALZHEIMER/AD00202")
org_counts <- read.csv("./预处理数据/FD1000/AD00203NormalLabel.csv",row.names=1)
org_label <- read.csv("./预处理数据/FD1000/label.csv",row.names=1,header=T)
gen_counts <- read.csv("./缩放后数据/Final2.0w(2000).csv",row.names=1,header=T)
gen_label <- read.csv("./缩放后数据/label.csv",row.names=1,header=T)
#将所有的标签从浮点型转变为int类型
org_label[, ncol(org_label)] <- as.integer(org_label[, ncol(org_label)]) - 1
gen_label[, ncol(gen_label)]<- as.integer(gen_counts[, ncol(gen_counts)]) + 5

#统一一列名
colnames(gen_counts) <- colnames(org_counts)
all_counts = rbind(org_counts[,1:ncol(org_counts)] , gen_counts[, 1:ncol(gen_counts)])
all_counts = all_counts[,1:ncol(all_counts)-1]
annotation = rbind(org_label, gen_label)

#原始数据进行PCA视觉化
counts <- t(org_counts[,1:1000])
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(counts)
    ,logcounts = as.matrix(counts)
  )
)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
#用PCA降为2维，然后画图
PCAsce <- runPCA(sce)
plotPCA(PCAsce,colour_by = org_label)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[,1:2], col = org_label[,1]+2,asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)


#生成数据进行PCA视觉化
counts <- t(gen_counts[,1:1000])
counts[1:4,1:4]
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(counts)
    ,logcounts = as.matrix(counts)
  )
)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
#用PCA降为2维，然后画图
PCAsce <- runPCA(sce)
plotPCA(PCAsce,colour_by = gen_label)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[,1:2], col = gen_label[,1]+1,asp=0.4,xlab="",ylab="",pch=19,cex=1.0,axes = FALSE)




#原始数据和生成数据放在一起进行PCA视觉化
counts <- t(all_counts)
##此处要经过归一化+log 列是细胞，行是基因，不包含labelhttp://127.0.0.1:11213/graphics/plot_zoom_png?width=1200&height=900
counts[1:4,1:4]
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(counts)
    ,logcounts = as.matrix(counts)
  )
)

rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

#用PCA降为2维，然后画图
PCAsce <- runPCA(sce)
#全部的数据
colors <- c("#777777", "blue")
annotation[,1] = 1
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1:1122,1:2], col = colors[annotation[1:1122,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1123:13122,1:2], col = colors[annotation[1123:13122,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[,1:2], col = colors[annotation[,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
#第1类
colors <- c("red","#777777")
annotation[,1] = 2
annotation[1:96,1] = 1
annotation[1123:3122,1]=1
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1:1122,1:2], col = colors[annotation[1:1122,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1123:13122,1:2], col = colors[annotation[1123:13122,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[,1:2], col = colors[annotation[,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
#第2类
colors <- c("blue","#777777")
annotation[,1] = 2
annotation[820:880,1] = 1
annotation[5123:7122,1]=1
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1:1122,1:2], col = colors[annotation[1:1122,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1123:13122,1:2], col = colors[annotation[1123:13122,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[,1:2], col = colors[annotation[,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)

#第3类
colors <- c("pink","#777777")
annotation[,1] = 2
annotation[1113:1123,1] = 1
annotation[11123:13122,1]=1
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1:1122,1:2], col = colors[annotation[1:1122,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1123:13122,1:2], col = colors[annotation[1123:13122,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[,1:2], col = colors[annotation[,1]],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)

#第4类
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1:1122,1:2], col = annotation[1:1122,1]+1,asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1123:13122,1:2], col = annotation[1123:13122,1],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)

#第5类
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1:1122,1:2], col = annotation[1:1122,1]+1,asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1123:13122,1:2], col = annotation[1123:13122,1],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)

#第6类
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1:1122,1:2], col = annotation[1:1122,1]+1,asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[1123:13122,1:2], col = annotation[1123:13122,1],asp=0.4,xlab="",ylab="",pch=19,cex=1.8,axes = FALSE)


plotPCA(PCAsce,colour_by = annotation)






#Rtsne降维(对原始的数据) 
#只能对matrix(行是细胞，列是基因)进行降维
counts <- org_counts[,1:ncol(org_counts)]
counts[1:4,1:4]
# 检查并删除重复的行
counts <- counts[!duplicated(counts), ]
label <- counts[,ncol(counts)]
counts <- counts[,1:ncol(counts)-1]
#对其除了label进行降维
tsne_matrix = Rtsne(as.matrix(counts))
tsne_matrix$Y[1:10,1:2]
#col 颜色是按照什么进行分类的？ asp:图像的长宽比  main：图的标题  pch：点的状态（19代表实心的点） cex：调节点的大小
plot(tsne_matrix$Y[,1:2], col = label,asp=0.4,main="Impute05IPSC-tSNE",xlab="Tsne1",ylab="Tsne2",pch=19,cex=0.5)


#Rtsne降维(对生成的数据) 
#只能对matrix(行是细胞，列是基因)进行降维
counts <- gen_counts[,1:ncol(gen_counts)]
counts[1:4,1:4]
# 检查并删除重复的行
counts <- counts[!duplicated(counts), ]
label <- counts[,ncol(counts)]
counts <- counts[,1:ncol(counts)-1]
#对其除了label进行降维
tsne_matrix = Rtsne(as.matrix(counts))
tsne_matrix$Y[1:10,1:2]
#col 颜色是按照什么进行分类的？ asp:图像的长宽比  main：图的标题  pch：点的状态（19代表实心的点） cex：调节点的大小
plot(tsne_matrix$Y[,1:2], col = label,asp=0.4,main="Impute05IPSC-tSNE",xlab="Tsne1",ylab="Tsne2",pch=19,cex=0.5)


#Rtsne降维(对源数据 + 生成的数据) 
#只能对matrix(行是细胞，列是基因)进行降维
all_counts = rbind(org_counts[,1:ncol(org_counts)] , gen_counts[, 1:ncol(gen_counts)])
counts <- all_counts
counts[1:4,1:4]
# 检查并删除重复的行
counts <- counts[!duplicated(counts), ]
label <- counts[,ncol(counts)]
counts <- counts[,1:ncol(counts)-1]
#对其除了label进行降维
tsne_matrix = Rtsne(as.matrix(counts))
tsne_matrix$Y[1:10,1:2]
#col 颜色是按照什么进行分类的？ asp:图像的长宽比  main：图的标题  pch：点的状态（19代表实心的点） cex：调节点的大小
plot(tsne_matrix$Y[,1:2], col = label,asp=0.4,main="Impute05IPSC-tSNE",xlab="Tsne1",ylab="Tsne2",pch=19,cex=0.5)









#UMAP降维(对原始的数据) 
#只能对matrix(行是细胞，列是基因)进行降维
counts <- org_counts[,1:ncol(org_counts)]
counts[1:4,1:4]
# 检查并删除重复的行
counts <- counts[!duplicated(counts), ]
label <- counts[,ncol(counts)]
counts <- counts[,1:ncol(counts)-1]
#对其除了label进行降维
umap1 <- umap(as.matrix(counts) , method = 'naive' , n_components = 20, n_neightbors=20)
tsne_matrix = Rtsne(as.matrix(counts))
tsne_matrix$Y[1:10,1:2]
#col 颜色是按照什么进行分类的？ asp:图像的长宽比  main：图的标题  pch：点的状态（19代表实心的点） cex：调节点的大小
plot(tsne_matrix$Y[,1:2], col = label,asp=0.4,main="Impute05IPSC-tSNE",xlab="Tsne1",ylab="Tsne2",pch=19,cex=0.5)

