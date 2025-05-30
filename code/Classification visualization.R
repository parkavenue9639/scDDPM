# This script performs visualization to assess the overlap between real and generated scRNA-seq data.
# It includes PCA, t-SNE, and UMAP dimensionality reduction on real, generated, and combined datasets.

# Load required packages
# install.packages("BiocManager")
# BiocManager::install("scater")
# BiocManager::install("SC3")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("umap")

library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library(scater)
library(Rtsne)
library(umap)

# Set working directory and clean environment
rm(list = ls())
setwd("F:/Graduate_Projects/DDPM/data/ALZHEIMER/AD00202")

# Load real and generated data
org_counts <- read.csv("./preprocessed/FD1000/AD00203NormalLabel.csv", row.names = 1)
org_label  <- read.csv("./preprocessed/FD1000/label.csv", row.names = 1, header = TRUE)
gen_counts <- read.csv("./scaled/Final2.0w(2000).csv", row.names = 1, header = TRUE)
gen_label  <- read.csv("./scaled/label.csv", row.names = 1, header = TRUE)

# Convert labels to integers
org_label[, ncol(org_label)] <- as.integer(org_label[, ncol(org_label)]) - 1
gen_label[, ncol(gen_label)] <- as.integer(gen_counts[, ncol(gen_counts)]) + 5

# Standardize column names and merge data
colnames(gen_counts) <- colnames(org_counts)
all_counts <- rbind(org_counts[, 1:ncol(org_counts)], gen_counts[, 1:ncol(gen_counts)])
all_counts <- all_counts[, -ncol(all_counts)]
annotation <- rbind(org_label, gen_label)

# -------- PCA visualization on real data --------
counts <- t(org_counts[, 1:1000])
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts = as.matrix(counts)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
PCAsce <- runPCA(sce)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[, 1:2],
     col = org_label[, 1] + 2, asp = 0.4, pch = 19, cex = 1.8, axes = FALSE)

# -------- PCA visualization on generated data --------
counts <- t(gen_counts[, 1:1000])
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts = as.matrix(counts)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
PCAsce <- runPCA(sce)
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[, 1:2],
     col = gen_label[, 1] + 1, asp = 0.4, pch = 19, cex = 1.0, axes = FALSE)

# -------- PCA on combined real and generated data --------
counts <- t(all_counts)
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts = as.matrix(counts)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
PCAsce <- runPCA(sce)

# Example: visualize selected cell types (class 1~6)
colors <- c("red", "#777777")
annotation[, 1] <- 2
annotation[1:96, 1] <- 1
annotation[1123:3122, 1] <- 1
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[, 1:2],
     col = colors[annotation[, 1]], asp = 0.4, pch = 19, cex = 1.8, axes = FALSE)

# -------- t-SNE on real data --------
counts <- org_counts[, -ncol(org_counts)]
label <- org_counts[, ncol(org_counts)]
counts <- counts[!duplicated(counts), ]
tsne_matrix <- Rtsne(as.matrix(counts))
plot(tsne_matrix$Y[, 1:2], col = label, asp = 0.4, pch = 19, cex = 0.5)

# -------- t-SNE on generated data --------
counts <- gen_counts[, -ncol(gen_counts)]
label <- gen_counts[, ncol(gen_counts)]
counts <- counts[!duplicated(counts), ]
tsne_matrix <- Rtsne(as.matrix(counts))
plot(tsne_matrix$Y[, 1:2], col = label, asp = 0.4, pch = 19, cex = 0.5)

# -------- t-SNE on combined data --------
all_counts <- rbind(org_counts, gen_counts)
counts <- all_counts[, -ncol(all_counts)]
label <- all_counts[, ncol(all_counts)]
counts <- counts[!duplicated(counts), ]
tsne_matrix <- Rtsne(as.matrix(counts))
plot(tsne_matrix$Y[, 1:2], col = label, asp = 0.4, pch = 19, cex = 0.5)

# -------- UMAP (20D) + t-SNE (2D) on real data --------
counts <- org_counts[, -ncol(org_counts)]
label <- org_counts[, ncol(org_counts)]
counts <- counts[!duplicated(counts), ]
umap1 <- umap(as.matrix(counts), method = "naive", n_components = 20, n_neighbors = 20)
tsne_matrix <- Rtsne(as.matrix(counts))
plot(tsne_matrix$Y[, 1:2], col = label, asp = 0.4, pch = 19, cex = 0.5)
