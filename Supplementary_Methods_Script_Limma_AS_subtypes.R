
#----load libraries-----
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("DESeq2")
library("tidyverse")
library("readxl")
library("ggplot2")
library("pheatmap")
library("ComplexHeatmap")
library("EnhancedVolcano")
library("vegan")
library("org.Hs.eg.db")
library("dplyr")
library("heatmap.plus")
library("heatmap.plus")
library("circlize")
library("pheatmap")
library("dplyr")
library("reshape2")
library("ggsignif")
library("Hmisc")
library("cluster")
library("scales")
library("dplyr")
library("viridis")
library("Polychrome")
library("matrixStats")
library("ConsensusClusterPlus")
library("vsn")
library("MSnbase")
library("DEP")
library("Biobase")
library("clusterProfiler")
library("sva")
library("httr")
library("jsonlite")
library("dplyr")
library("tibble")
library("gprofiler2")
library("msigdbr")
library("umap")

counts <- read.csv("AS_count_data_cluster_proteins.csv", sep=";", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))

col<- read.csv("AS_col_data_cluster.csv", sep=";", row.names = 1)

common_samples <- intersect(colnames(counts), rownames(col))
exprs <- counts[, common_samples, drop = FALSE]
meta  <- col[common_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta)) 


#-----------set wd----------
setwd("D:/SARCOMA_RESEARCH/CSF1R_Angiosarcoma_Köln_Fortune/Proteomic_analyses/Results_proteomics")


#---------UMAP proteomics data by AS subtype-------
valid_samples <- rownames(meta[meta$AS_subtype_simple %in% c(1,2,3,5), ]) # include only main subtypes & benign tissue
exprs <- exprs[, valid_samples, drop = FALSE]
meta  <- meta[valid_samples, , drop = FALSE]

X <- t(exprs)
X_scaled <- scale(X)
set.seed(42)
umap_result <- umap(X_scaled, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

umap_df <- as.data.frame(umap_result$layout)

colnames(umap_df) <- c("UMAP1", "UMAP2")
meta <- meta %>% 
  rownames_to_column(var = "Sample")

umap_df <- umap_df %>%
  mutate(Sample = rownames(X)) %>%
  left_join(meta, by = "Sample")

meta$row.names<-rownames(meta)
umap_df <- umap_df %>%
  mutate(Sample = rownames(X)) %>%
  left_join(meta, by = c("Sample" = "row.names"))

umap_df$AS_subtype_simple.x<-as.factor(umap_df$AS_subtype_simple.x)

cluster_colors <- c(
  "1" = "#F4A582", # 1 = primary cutaneous AS
  "2"="#67001F", # 2 = secondary AS
  "3"="#E69F00", # 3 = visceral AS
  "5"="skyblue") # 5 = benign tissue

png("Figure_3_UMAP_plot_all_AS.png", height = 7, width=9, units="cm", res=600)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = umap_df$AS_subtype_simple.x)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  theme_minimal(base_size = 14) +
  labs(title = NULL) +
  theme_light()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )
dev.off()


#------LIMMA Primary AS vs. benign tissue----------
counts <- read.csv("AS_count_data_cluster_proteins.csv", sep=";", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))

col<- read.csv("AS_col_data_cluster.csv", sep=";", row.names = 1)

common_samples <- intersect(colnames(counts), rownames(col))
exprs <- counts[, common_samples, drop = FALSE]
meta  <- col[common_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta)) 

valid_samples <- rownames(meta[meta$AS_subtype_simple %in% c(1, 5), ])

exprs <- exprs[, valid_samples, drop = FALSE]
meta  <- meta[valid_samples, , drop = FALSE]

all(colnames(exprs) == rownames(meta)) 

meta$condition <- factor(meta$AS_subtype_simple,
                         levels = c(1, 5),
                         labels = c("Primary_AS", "Benign"))
design <- model.matrix(~0 + condition, data = meta)
colnames(design) <- c("Primary_AS", "Benign")

contrast.matrix <- makeContrasts(Primary_AS_vs_Benign = Primary_AS - Benign, levels = design)

fit <- lmFit(exprs, design)
fit2 <- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit2, trend = TRUE)      
res <- topTable(fit2, coef="Primary_AS_vs_Benign", number=Inf, adjust.method="BH")

res_unique <- res %>%
  dplyr::filter(!is.na(rownames(res))) %>%     # optional: drop rows with missing protein names
  dplyr::distinct(rownames(res), .keep_all = TRUE)

# Define log2FC thresholds
log2FC_up <- 1.0
log2FC_down <- -1.0
pval_threshold <- 0.05

upregulated <- res_unique[res_unique$logFC >= log2FC_up & res_unique$adj.P.Val < pval_threshold, ]
downregulated <- res_unique[res_unique$logFC <= log2FC_down & res_unique$adj.P.Val < pval_threshold, ]

cat("Number of upregulated proteins (log2FC >= 1.0):", nrow(upregulated), "\n")
cat("Number of downregulated proteins (log2FC <= -1.0):", nrow(downregulated), "\n")

res_unique$Regulation <- "NS" 
res_unique$Regulation[res_unique$logFC >= log2FC_up & res_unique$adj.P.Val < pval_threshold] <- "Up"
res_unique$Regulation[res_unique$logFC <= log2FC_down & res_unique$adj.P.Val < pval_threshold] <- "Down"


#------LIMMA Secondary AS vs. benign tissue ----------
counts <- read.csv("AS_count_data_cluster_proteins.csv", sep=";", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))

col<- read.csv("AS_col_data_cluster.csv", sep=";", row.names = 1)

common_samples <- intersect(colnames(counts), rownames(col))
exprs <- counts[, common_samples, drop = FALSE]
meta  <- col[common_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta)) 
valid_samples <- rownames(meta[meta$AS_subtype_simple %in% c(2, 5), ])
exprs <- exprs[, valid_samples, drop = FALSE]
meta  <- meta[valid_samples, , drop = FALSE]

all(colnames(exprs) == rownames(meta)) 

meta$condition <- factor(meta$AS_subtype_simple,
                         levels = c(2, 5),
                         labels = c("Secondary_AS", "Benign"))
design <- model.matrix(~0 + condition, data = meta)
colnames(design) <- c("Secondary_AS", "Benign")

contrast.matrix <- makeContrasts(Secondary_AS_vs_Benign = Secondary_AS - Benign, levels = design)

fit <- lmFit(exprs, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)
res <- topTable(fit2, coef="Secondary_AS_vs_Benign", number=Inf, adjust.method="BH")

res_unique <- res %>%
  dplyr::filter(!is.na(rownames(res))) %>%  
  dplyr::distinct(rownames(res), .keep_all = TRUE)

#Define log2FC thresholds
log2FC_up <- 1.0
log2FC_down <- -1.0
pval_threshold <- 0.05

upregulated <- res_unique[res_unique$logFC >= log2FC_up & res_unique$adj.P.Val < pval_threshold, ]
downregulated <- res_unique[res_unique$logFC <= log2FC_down & res_unique$adj.P.Val < pval_threshold, ]

cat("Number of upregulated proteins (log2FC >= 1.0):", nrow(upregulated), "\n")
cat("Number of downregulated proteins (log2FC <= -1.0):", nrow(downregulated), "\n")

res_unique$Regulation <- "NS"  # NS = not significant
res_unique$Regulation[res_unique$logFC >= log2FC_up & res_unique$adj.P.Val < pval_threshold] <- "Up"
res_unique$Regulation[res_unique$logFC <= log2FC_down & res_unique$adj.P.Val < pval_threshold] <- "Down"


#------LIMMA Visceral AS vs. benign tissue ----------
counts <- read.csv("AS_count_data_cluster_proteins.csv", sep=";", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))

col<- read.csv("AS_col_data_cluster.csv", sep=";", row.names = 1)

common_samples <- intersect(colnames(counts), rownames(col))
exprs <- counts[, common_samples, drop = FALSE]
meta  <- col[common_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta)) 
valid_samples <- rownames(meta[meta$AS_subtype_simple %in% c(3, 5), ])
exprs <- exprs[, valid_samples, drop = FALSE]
meta  <- meta[valid_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta)) 

meta$condition <- factor(meta$AS_subtype_simple,
                         levels = c(3, 5),
                         labels = c("Visceral_AS", "Benign"))
design <- model.matrix(~0 + condition, data = meta)
colnames(design) <- c("Visceral_AS", "Benign")

contrast.matrix <- makeContrasts(Visceral_AS_vs_Benign = Visceral_AS - Benign, levels = design)

fit <- lmFit(exprs, design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2, trend = TRUE)  
res <- topTable(fit2, coef="Visceral_AS_vs_Benign", number=Inf, adjust.method="BH")

res_unique <- res %>%
  dplyr::filter(!is.na(rownames(res))) %>%     
  dplyr::distinct(rownames(res), .keep_all = TRUE)

#Define log2FC thresholds
log2FC_up <- 1.0
log2FC_down <- -1.0
pval_threshold <- 0.05 

upregulated <- res_unique[res_unique$logFC >= log2FC_up & res_unique$adj.P.Val < pval_threshold, ]
downregulated <- res_unique[res_unique$logFC <= log2FC_down & res_unique$adj.P.Val < pval_threshold, ]

cat("Number of upregulated proteins (log2FC >= 1.0):", nrow(upregulated), "\n")
cat("Number of downregulated proteins (log2FC <= -1.0):", nrow(downregulated), "\n")

res_unique$Regulation <- "NS"  # NS = not significant
res_unique$Regulation[res_unique$logFC >= log2FC_up & res_unique$adj.P.Val < pval_threshold] <- "Up"
res_unique$Regulation[res_unique$logFC <= log2FC_down & res_unique$adj.P.Val < pval_threshold] <- "Down"

