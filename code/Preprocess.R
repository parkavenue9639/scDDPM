# ================== Environment Setup ==================
rm(list = ls())

# Required packages (uncomment to install if needed)
# install.packages("Seurat")
# devtools::install_version("Seurat", version = "4.1.4")
# BiocManager::install("Seurat")
# install.packages("Matrix")


# ================== Load Data ==================
# Load gene expression matrix (rows: genes, columns: cells)
data <- read.delim("./AD01103_expr.txt", header = TRUE, row.names = 1)

# Calculate non-zero rate
completed_nz <- sum(data > 0) / (nrow(data) * ncol(data))
cat("✅ Non-zero rate:", round(completed_nz * 100, 2), "%\n")

library(Seurat)
library(SeuratObject)
library(magrittr)
library(dplyr)

# Basic statistics
cat("🔢 Dimensions: ", dim(data)[1], "genes ×", dim(data)[2], "cells\n")
cat("🧪 Total zero entries:", sum(data == 0), "\n")

# ================== Create Seurat Object ==================
# Rows: genes, Columns: cells. Apply minimum cell filter (5%)
scell <- CreateSeuratObject(counts = data, project = "Alzheimer0205",
                            min.cells = ceiling(0.05 * ncol(data)))
data <- 0  # Free memory

# Save gene and cell names
gene_names <- rownames(scell)
cell_names <- colnames(scell)
cat("🧬 Filtered dimensions:", dim(scell@assays$RNA@layers$counts), "\n")

# ================== Normalize Data ==================
# Log-normalization (CPM → log1p)
data <- NormalizeData(scell, normalization.method = "LogNormalize", scale.factor = 10000)
norm <- data@assays$RNA@layers$data
normalized_df <- as.data.frame(norm)
rownames(normalized_df) <- gene_names
colnames(normalized_df) <- cell_names

# ================== Select Top Highly Variable Genes ==================
# Select top 100 variable genes (vst method)
data1000 <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 100)
top <- head(VariableFeatures(data1000), 100)
pb <- normalized_df[top, ]
pb <- t(pb)
cat("✅ Selected top variable genes:", dim(pb), "\n")

# Save processed matrix (expression only)
write.csv(pb, "../FD100/AD01103PrePro100.csv")

# ================== Merge with Cell Labels ==================
label <- read.delim("AD01103_cell_label.txt", header = TRUE, row.names = 1)
merged_data <- merge(pb, label, by = "row.names", all.x = TRUE)

cat("🔖 Label distribution:\n")
print(table(label))

# Save final labeled matrix
write.csv(merged_data, "../FD100/AD01103PreProLabel100.csv", row.names = FALSE)
cat("✅ Processing complete. Labeled data saved.\n")
