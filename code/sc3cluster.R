# ====================== Package Setup ======================
# Uncomment the lines below if you haven't installed the required packages
# BiocManager::install("scater")
# BiocManager::install("SC3")

library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library(scater)
library(Rtsne)

# ====================== File Path Setup ======================
getwd()
setwd("F:/PythonCoding/Jupyter/DDPM/data/output/DiffusionCell/AD0203/fea100")
rm(list = ls())

# ====================== Load Expression Data ======================
# Load generated gene expression matrix (rows: genes, columns: cells)
counts <- read.csv("outputws2.0.csv", row.names = 1)
cat("✅ Raw matrix dimensions: ", dim(counts)[1], "genes ×", dim(counts)[2], "cells\n")

# Select first 100 cells for visualization/clustering
counts <- counts[, 1:100]

# Load annotation file (e.g., predicted labels)
annotation <- read.csv("gen_label2.csv", row.names = 1)

# Transpose expression matrix (SC3 requires genes as rows, cells as columns)
counts <- t(counts)
cat("✅ Transposed matrix dimensions:", dim(counts)[1], "cells ×", dim(counts)[2], "genes\n")

# ====================== Create SingleCellExperiment Object ======================
# Input matrix must be normalized and log-transformed
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(counts),
    logcounts = as.matrix(counts)
  ),
  colData = annotation
)

# Remove duplicated gene symbols
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
cat("✅ SC3 input matrix after removing duplicates:", dim(sce), "\n")

# ====================== Run SC3 Clustering ======================
# Perform SC3 clustering with k = 6 clusters
sce <- sc3(sce, ks = 6, biology = TRUE, gene_filter = FALSE)
cat("✅ SC3 clustering complete.\n")

# ====================== Save Results ======================
write.csv(colData(sce), "./sc3_2.0_label.csv")
cat("✅ SC3 clustering labels saved to 'sc3_2.0_label.csv'\n")
