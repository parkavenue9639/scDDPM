# README

This repository contains the implementation and downstream analysis of the scDDPM model for single-cell data generation and evaluation.

## 🔧 Step 1: Data Preprocessing

Run the preprocessing script written in R to generate processed input data for each dataset:

```
code/Preprocess.R
```

This script prepares the datasets for generative modeling by applying normalization, filtering, and transformation procedures.

## 🧠 Step 2: Data Generation with scDDPM

Use the following Python script to generate synthetic single-cell data using the scDDPM model:

```
code/scDDPM.py
```

This script takes the preprocessed data as input and performs model training, inference, and post-processing to produce the final generated data.

## 📊 Step 3: Downstream Analysis

After data generation, perform the following downstream tasks:

### 3.1 Clustering

Run SC3 clustering using:

```
code/sc3cluster.R
```

### 3.2 Visualization

Visualize the clustering results using:

```
code/Classification visualization.R
```

### 3.3 Differential Expression Analysis

Identify the top 100 differentially expressed genes (DEGs) using:

```
code/differential gene expression.R
```

### 3.4 Pathway Enrichment (KEGG)

Compare enriched biological pathways between real and generated data using:

```
code/kegg.R
```
