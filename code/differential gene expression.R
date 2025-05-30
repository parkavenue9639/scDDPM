# Clear environment
rm(list = ls())

# Set working directory


# Load original and generated data
data <- read.csv("./AD00204NormalLabel.csv", header = TRUE, row.names = 1)
gen_data <- read.csv("../../scaled/Final2.0w.csv", header = TRUE, row.names = 1)

# Load necessary libraries
library(Seurat)
library(SeuratObject)
library(magrittr)
library(dplyr)

# Check sparsity
sum(data == 0)
data[1:5, 1:5]
dim(data)

# Select the first 100 cells and transpose: Seurat requires features as rows and cells as columns
data <- data[, 1:100]
data <- t(data)

# Create Seurat object for original data
scell <- CreateSeuratObject(counts = data, project = "Alzheimer0205")

# Check Seurat object's internal structure (optional)
dim(scell@assays$RNA@layers$counts)

# Identify top 100 highly variable genes using VST method
data1000 <- FindVariableFeatures(scell, selection.method = "vst", nfeatures = 100)
top <- head(VariableFeatures(data1000), 100)

# Save top genes to CSV
write.csv(top, "../NormalTop100.csv", row.names = FALSE)

# ------------------ Repeat for Generated Data ------------------

gen_data[1:5, 1:5]
dim(gen_data)

# Select and transpose generated data
gen_data <- gen_data[, 1:100]
gen_data <- t(gen_data)

# Create Seurat object for generated data
gen_scell <- CreateSeuratObject(counts = gen_data, project = "Alzheimer0205")

# Check structure
dim(gen_scell@assays$RNA@layers$counts)

# Identify top 100 variable genes
gen_data100 <- FindVariableFeatures(gen_scell, selection.method = "vst", nfeatures = 100)
top <- head(VariableFeatures(gen_data100), 100)

# Save top genes
write.csv(top, "../GenTop100.csv", row.names = FALSE)
