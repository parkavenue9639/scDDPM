# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(enrichplot)
library(VennDiagram)

# ========== Load Data ==========
# Load real and generated expression data
expr_real <- read.csv("AD00203PreProLabel1000.csv", check.names = FALSE)
rownames(expr_real) <- make.unique(as.character(expr_real[, 1]))
expr_real <- expr_real[, -1]

expr_gen <- read.csv("Final2.0w.csv", check.names = FALSE)
rownames(expr_gen) <- make.unique(as.character(expr_gen[, 1]))
expr_gen <- expr_gen[, -1]

# ========== Extract Labels ==========
labels_real <- as.character(expr_real$label)
expr_real$label <- NULL
expr_real <- t(expr_real)

labels_gen <- as.character(expr_gen$label)
expr_gen$label <- NULL
expr_gen <- t(expr_gen)

# Set the label pairs to compare (real vs generated)
real_label1 <- "1"; real_label2 <- "5"
gen_label1  <- "0"; gen_label2  <- "4"

# Find column indices for each group
group1_real <- which(labels_real == real_label1)
group2_real <- which(labels_real == real_label2)
group1_gen  <- which(labels_gen == gen_label1)
group2_gen  <- which(labels_gen == gen_label2)

cat("Real data columns: ", length(group1_real), "vs", length(group2_real), "\n")
cat("Generated data columns: ", length(group1_gen), "vs", length(group2_gen), "\n")

if (any(sapply(list(group1_real, group2_real, group1_gen, group2_gen), length) == 0)) {
  stop("❌ One or more label groups are empty. Please check label assignments.")
}

# ========== Differential Expression Analysis ==========
get_deg <- function(expr_mat, idx1, idx2) {
  pvals <- apply(expr_mat, 1, function(x) {
    tryCatch(wilcox.test(x[idx1], x[idx2])$p.value, error = function(e) NA)
  })
  logfc <- rowMeans(expr_mat[, idx1, drop = FALSE]) - rowMeans(expr_mat[, idx2, drop = FALSE])
  deg_df <- data.frame(gene = rownames(expr_mat), logFC = logfc, pval = pvals)
  return(deg_df)
}

deg_real <- get_deg(expr_real, group1_real, group2_real)
deg_gen  <- get_deg(expr_gen, group1_gen, group2_gen)

# ========== Gene ID Conversion ==========
gene_real <- bitr(deg_real$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_gen  <- bitr(deg_gen$gene,  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cat("🧬 DEGs (real): ", nrow(deg_real), "\n")
cat("🧬 DEGs (generated): ", nrow(deg_gen), "\n")
cat("✅ Mapped ENTREZ IDs (real): ", nrow(gene_real), "\n")
cat("✅ Mapped ENTREZ IDs (generated): ", nrow(gene_gen), "\n")

if (nrow(gene_real) == 0 || nrow(gene_gen) == 0) {
  stop("❌ No valid ENTREZ IDs found. Please verify gene names.")
}

# ========== KEGG Enrichment ==========
kegg_real <- enrichKEGG(gene = gene_real$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
kegg_gen  <- enrichKEGG(gene = gene_gen$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)

if (is.null(kegg_real) || is.null(kegg_gen)) {
  stop("❌ enrichKEGG returned NULL. Check if valid genes are recognized by KEGG.")
}

# ========== Pathway Comparison ==========
real_terms <- kegg_real@result$Description
gen_terms  <- kegg_gen@result$Description

common_terms <- intersect(real_terms, gen_terms)
only_real    <- setdiff(real_terms, gen_terms)
only_gen     <- setdiff(gen_terms, real_terms)

write.csv(data.frame(Common = common_terms), "common_kegg_pathways.csv", row.names = FALSE)

# ========== Visualization ==========
pdf("KEGG_comparison_and_venn.pdf", width = 14, height = 8)

# Side-by-side dotplots
par(mfrow = c(1, 2))
print(dotplot(kegg_real, showCategory = 15, title = "Real Data KEGG"))
print(dotplot(kegg_gen, showCategory = 15, title = "Generated Data KEGG"))

# Venn diagram
grid.newpage()
draw.pairwise.venn(
  area1       = length(real_terms),
  area2       = length(gen_terms),
  cross.area  = length(common_terms),
  category    = c("Real", "Generated"),
  fill        = c("skyblue", "pink"),
  alpha       = 0.6,
  cat.pos     = c(-20, 20),
  cat.dist    = 0.05,
  scaled      = FALSE
)

dev.off()
cat("✅ Analysis complete. Results saved to KEGG_comparison_and_venn.pdf and common_kegg_pathways.csv\n")
