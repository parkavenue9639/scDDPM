#BiocManager::install("sacter")
#BiocManager::install("SC3")

library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library(scater)
library(Rtsne)
getwd()
setwd("F:/PythonCoding/Jupyter/DDPM/data/output/DiffusionCell/AD0203/fea100")
rm(list =ls())
counts <- read.csv("F:/PythonCoding/Jupyter/DDPM/data/output/DiffusionCell/AD0203/fea100/outputws2.0.csv",row.names=1)
dim(counts)
counts = counts[,1:100]
annotation <- read.csv("F:/PythonCoding/Jupyter/DDPM/data/output/DiffusionCell/AD0203/fea100/gen_label2.csv",row.names =1)


counts <- t(counts)

##此处要经过归一化+log 列是细胞，行是基因，不包含label
counts[1:4,1:4]
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(counts)
    ,logcounts = as.matrix(counts)
  ), 
  colData = annotation
)

rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]


#开始做SC3聚类
sce <- sc3(sce, ks = 6 , biology = TRUE,gene_filter = FALSE)

write.csv(colData(sce) , "./sc3_2.0_label.csv")