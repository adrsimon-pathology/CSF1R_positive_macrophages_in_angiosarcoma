
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

#--------------set wd----------------
setwd("D:/SARCOMA_RESEARCH/CSF1R_Angiosarcoma_Köln_Fortune/Proteomic_analyses/Results_proteomics")



#-------LIMMA Cluster 1 vs. other------------
counts <- read.csv("AS_count_data_cluster_proteins.csv", sep=";", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))

col<- read.csv("AS_col_data_cluster.csv", sep=";", row.names = 1)

common_samples <- intersect(colnames(counts), rownames(col))
exprs <- counts[, common_samples, drop = FALSE]
meta  <- col[common_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta)) 
table(meta$AS_subtype_simple)
valid_samples <- rownames(meta)[!is.na(meta$Cluster1_vs_other)]
exprs <- exprs[, valid_samples, drop = FALSE]
meta <- meta[valid_samples, , drop = FALSE]

all(colnames(exprs) == rownames(meta)) 

meta$condition <- factor(meta$Cluster)
levels(meta$condition) <- c("Cluster1", "Other")
table(meta$condition)

design <- model.matrix(~0 + condition, data = meta)
colnames(design) <- c("Cluster1", "Other")

contrast.matrix <- makeContrasts(Cluster1_vs_Other = Cluster1 - Other, levels = design)

fit <- lmFit(exprs, design)
fit2 <- contrasts.fit(fit, contrast.matrix)  # apply contrast
fit2 <- eBayes(fit2, trend = TRUE)        # compute moderated stats
res <- topTable(fit2, coef="Cluster1_vs_Other", number=Inf, adjust.method="BH")

#Define log2FC thresholds
log2FC_up <- 1.0
log2FC_down <- -1.0
pval_threshold <- 0.05

upregulated <- res_unique[res_unique$logFC >= log2FC_up & res_unique$P.Value < pval_threshold, ]
downregulated <- res_unique[res_unique$logFC <= log2FC_down & res_unique$P.Value < pval_threshold, ]

cat("Number of upregulated proteins (log2FC >= 1.0):", nrow(upregulated), "\n")
cat("Number of downregulated proteins (log2FC <= -1.0):", nrow(downregulated), "\n")

res_unique$Regulation <- "NS"  # NS = not significant
res_unique$Regulation[res_unique$logFC >= log2FC_up & res_unique$P.Value < pval_threshold] <- "Up"
res_unique$Regulation[res_unique$logFC <= log2FC_down & res_unique$P.Value < pval_threshold] <- "Down"

table(res_unique$Regulation)
table(meta$AS_subtype_simple)


#----------Figure 4B UMAP Cluster 1 vs. other---------------------
counts <- read.csv("AS_count_data_cluster.csv", sep=";", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))
col <- read.csv("AS_col_data_cluster.csv", sep=";", row.names = 1)

common_samples <- intersect(colnames(counts), rownames(col))
exprs <- counts[, common_samples, drop = FALSE]
meta  <- col[common_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta))

table(meta$AS_subtype_simple)
valid_samples <- rownames(meta)[!is.na(meta$Cluster1_vs_other)]
exprs <- exprs[, valid_samples, drop = FALSE]
meta <- meta[valid_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta))
X <- t(exprs)
X_scaled <- scale(X)

set.seed(42)
umap_result <- umap(X_scaled, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
meta <- meta %>% 
  rownames_to_column(var = "Sample")  # convert row names → column "Sample"

umap_df <- umap_df %>%
  mutate(Sample = rownames(X)) %>%
  left_join(meta, by = "Sample")

meta$row.names<-rownames(meta)
umap_df <- umap_df %>%
  mutate(Sample = rownames(X)) %>%
  left_join(meta, by = c("Sample" = "row.names"))

umap_df$AS_subtype_simple.x<-as.factor(umap_df$AS_subtype_simple.x)
umap_df$Cluster.x<-as.factor(umap_df$Cluster.x)

cluster_colors <- c(
  "1" = "black",
  "2" = "darkmagenta",   # orange
  "3" = "darkorange")

png("Figure_3_UMAP_plot.png", height = 7, width=9, units="cm", res=600)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster.x)) +
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


#----------Figure 4A: HEATMAP Clusters-----------
counts <- read.csv("AS_count_data_cluster_proteins.csv", sep=";", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))
col <- read.csv("AS_col_data_cluster.csv", sep=";", row.names = 1)
common_samples <- intersect(colnames(counts), rownames(col))
exprs <- counts[, common_samples, drop = FALSE]
meta  <- col[common_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta))

table(meta$AS_subtype_simple)
valid_samples <- rownames(meta)[!is.na(meta$Cluster1_vs_other)]
exprs <- exprs[, valid_samples, drop = FALSE]
meta <- meta[valid_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta))

var_genes <- apply(exprs, 1, var, na.rm = TRUE)
top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:6885]
exprs_top <- exprs[top_var_genes, ]
-
exprs_scaled <- t(scale(t(exprs_top)))

colsplit <- data.frame(Cluster = meta$Cluster)
rownames(colsplit) <- rownames(meta)

cluster_colors <- c(
  "1" = "black",
  "2" = "darkmagenta",
  "3" = "darkorange"
)
levels_to_colors <- list(Cluster = cluster_colors)

ha <- HeatmapAnnotation(df = colsplit, col = levels_to_colors, na_col = "white")

bottom_annot_df <- data.frame(
  Histology = factor(meta$AS_subtype_simple,
                     levels = sort(unique(meta$AS_subtype_simple)))
)
rownames(bottom_annot_df) <- rownames(meta)

score_col_fun <- list(
  Histology = c(
    "1" = "#F4A582",
    "2" = "#67001F",
    "3" = "#E69F00",
    "4" = "#666666"
  )
)
bottom_annot <- HeatmapAnnotation(df = bottom_annot_df,
                                  col = score_col_fun,
                                  na_col = "white")

table(meta$Cluster)
exprs_col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("navy", "white", "red")
)

png(file = "Figure_Proteomics_Heatmap_Top200.png", 
    width = 21, height = 7, units = "cm", res = 600)

heat <- Heatmap(exprs_scaled,
                name = "Z-score",
                top_annotation = ha,
                bottom_annotation = bottom_annot,
                column_split = colsplit,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                show_row_names = FALSE,
                show_row_dend = FALSE, 
                show_column_names = FALSE,
                row_names_side = "left",
                column_title = NULL,
                row_title = NULL,
                col = exprs_col_fun,
                border = TRUE)

draw(heat)

dev.off()

table(meta$Cluster)


#---------------Figure 4C Volcano plot---------
limma_res <- res_unique %>%
  mutate(Significant = case_when(
    logFC > 1.0 & P.Value < 0.05                  ~ "Up (log2FC > 1 & P < 0.05)",
    logFC <= -1.0 & P.Value < 0.05               ~ "Down (log2FC <= -1 & P < 0.05)",
    logFC > -1.0 & logFC < 1.0 & P.Value < 0.05  ~ "Nominal (|log2FC| < 1 & P < 0.05)",
    TRUE                                         ~ "NS"
  ))

png(file="Volcano_Proteomics_ThreeColors.png", height=8, width=20, units="cm", res=600)
ggplot(limma_res, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Up (log2FC > 1 & P < 0.05)" = "red",
    "Down (log2FC <= -1 & P < 0.05)" = "blue",
    "Nominal (|log2FC| < 1 & P < 0.05)" = "lightgoldenrod",
    "NS" = "grey"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = NULL,
    x = "log2 Fold Change",
    y = "-log10(P-value)",
    color = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )
dev.off()



# --- 1. Identify significant up/downregulated proteins ---
counts <- read.csv("AS_count_data_cluster_proteins.csv", sep=";", row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))
col <- read.csv("AS_col_data_cluster.csv", sep=";", row.names = 1)

# --- Match samples ---
common_samples <- intersect(colnames(counts), rownames(col))
exprs <- counts[, common_samples, drop = FALSE]
meta  <- col[common_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta))

table(meta$AS_subtype_simple)

# --- Keep only samples with non-NA AS_subtype_simple ---
valid_samples <- rownames(meta)[!is.na(meta$Cluster1_vs_other)]
exprs <- exprs[, valid_samples, drop = FALSE]
meta <- meta[valid_samples, , drop = FALSE]
all(colnames(exprs) == rownames(meta))

# --- 1. Identify significant up/downregulated proteins ---
sig_genes <- res_unique %>%
  filter(P.Value <= 0.05 & (logFC >= 1.0 | logFC <= -1.0)) %>%
  arrange(desc(abs(logFC)))

cat("Significant proteins:", nrow(sig_genes), "\n")

# --- 2. Extract expression data for these proteins ---
exprs_sig <- exprs[rownames(exprs) %in% sig_genes$Protein.Groups, ]
matching <- match(rownames(exprs_sig), sig_genes$Protein.Groups)

rownames(exprs_sig) <- sig_genes$protein.name[matching]

# Sanity check
head(rownames(exprs_sig))


# --- 3. Scale expression per protein for visualization ---
exprs_sig_scaled <- t(scale(t(exprs_sig)))

# --- 4. Build annotation dataframe ---
col_annot <- data.frame(Cluster = meta$Cluster1_vs_other)
rownames(col_annot) <- rownames(meta)

# Define color scheme for clusters
cluster_colors <- list(Cluster = c(
  "1" = "black",
  "2" = "grey50"
))

ha_top <- HeatmapAnnotation(df = col_annot, col = cluster_colors, na_col = "white")

# --- 5. Bottom annotation ---
bottom_annot_df <- data.frame(
  Histology = factor(meta$AS_subtype_simple,
                     levels = sort(unique(meta$AS_subtype_simple)))
)
rownames(bottom_annot_df) <- rownames(meta)

score_col_fun <- list(
  Histology = c(
    "1" = "#F4A582",
    "2" = "#67001F",
    "3" = "#E69F00",
    "4" = "grey"
  )
)
bottom_annot <- HeatmapAnnotation(df = bottom_annot_df,
                                  col = score_col_fun,
                                  na_col = "white")

# --- 6. Define color function for expression values ---
exprs_col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "red3"))

# --- 7. Draw heatmap ---
png(file = "Heatmap_Significant_Proteins_LeftLabels_Legend.png",
    width = 25, height = 17, units = "cm", res = 600)

heat_sig <- Heatmap(exprs_sig_scaled,
                    name = "Z-score",
                    top_annotation = ha_top,
                    bottom_annotation = bottom_annot,
                    column_split = col_annot,
                    cluster_columns = FALSE,
                    cluster_rows = TRUE,
                    show_row_dend = FALSE,
                    show_row_names = TRUE,     # show gene names
                    row_names_side = "left",   # move names to the left
                    show_column_names = FALSE,
                    column_title = NULL,
                    row_title = NULL,
                    col = exprs_col_fun,
                    border = TRUE,
                    heatmap_legend_param = list(title = NULL))  # optional: empty title

draw(heat_sig,
     show_heatmap_legend = TRUE,
     show_annotation_legend = TRUE,
     padding = unit(c(5, 70, 2, 5), "mm"))  # top, right, bottom, left

dev.off()

# Install if not already installed
install.packages("gprofiler2")

msig_h <- msigdbr(species = "Homo sapiens", category = "H")

# Convert to TERM2GENE format for clusterProfiler::enricher()
msig_h_list <- msig_h %>%
  select(gs_name, gene_symbol) %>%
  as.data.frame()

cluster1_data <- read.csv("Cluster_1_up.csv", header = TRUE, stringsAsFactors = FALSE)

# --- Step 2: Extract the list of genes/proteins ---
cluster_proteins <- cluster1_data$Gene

# Remove any missing values or duplicates
cluster_proteins <- unique(na.omit(cluster_proteins))

# Check what the list looks like
print(cluster_proteins)

gost_results <- gost(
  query = cluster_proteins,
  organism = "hsapiens",                # Change if not human
  sources = c("REAC", "KEGG", "GO:BP"), # Biological Process, KEGG, Reactome
  correction_method = "fdr",
  user_threshold = 0.05
)
# View results
head(gost_results$result)

# Make a copy
gost_results_flat <- gost_results$result

# Identify which columns are lists
list_cols <- sapply(gost_results_flat, is.list)

# Flatten every list column into comma-separated strings
gost_results_flat[list_cols] <- lapply(
  gost_results_flat[list_cols],
  function(col) sapply(col, function(x) paste(x, collapse = ","))
)

# Write to CSV safely
write.csv(gost_results_flat, "cluster1_up_gprofiler_enrichment.csv", row.names = FALSE)

# --- 1️⃣ Separate g:Profiler enrichment results ---
go_bp    <- gost_results_flat %>% filter(source == "GO:BP")
kegg     <- gost_results_flat %>% filter(source == "KEGG")
reactome <- gost_results_flat %>% filter(source == "REAC")

# Convert to TERM2GENE format for enricher()
msig_h_list <- msig_h %>%
  select(gs_name, gene_symbol) %>%
  as.data.frame()

# Run enrichment on your upregulated genes
ego_hallmark <- enricher(
  gene          = cluster_proteins,   # your input vector of gene symbols
  TERM2GENE     = msig_h_list,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Extract results
hallmark_df <- as.data.frame(ego_hallmark) %>%
  mutate(source = "MSigDB Hallmark")

# --- 3️⃣ Select top 3 from each source ---
go_top <- go_bp %>%
  arrange(p_value) %>%
  slice_head(n = 2) %>%
  mutate(source = "GO:BP")

kegg_top <- kegg %>%
  arrange(p_value) %>%
  slice_head(n = 2) %>%
  mutate(source = "KEGG")

reactome_top <- reactome %>%
  arrange(p_value) %>%
  slice_head(n = 2) %>%
  mutate(source = "REAC")

hallmark_top <- hallmark_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 2)

# --- 4️⃣ Combine all sources ---
combined_top <- bind_rows(go_top, kegg_top, reactome_top, hallmark_top) %>%
  mutate(
    term_label = case_when(
      !is.na(term_name) ~ paste0(term_name, " (", source, ")"),
      !is.na(Description) ~ paste0(Description, " (", source, ")"),
      TRUE ~ NA_character_
    ),
    adj_p = ifelse(!is.na(p_value), p_value, p.adjust),
    neg_log10_p = -log10(adj_p)
  ) %>%
  arrange(adj_p)

# --- 5️⃣ Define colors ---
source_colors <- c(
  "GO:BP"            = "darkred",  # firebrick (dark, classic red)
  "KEGG"             = "indianred",  # cinnabar
  "REAC"             = "#FF7F50",  # coral
  "MSigDB Hallmark"  = "red"   # brownish red
)

# --- 6️⃣ Plot unified bar chart ---
png(file="Cluster1_up.png", height=6, width=20, units="cm", res=600)

ggplot(combined_top,
       aes(x = reorder(term_label, neg_log10_p),
           y = neg_log10_p,
           fill = source)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = source_colors) +
  labs(
    x = NULL,
    y = expression(-log[10](adjusted~p)),
    title = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    legend.position = "none"    # ⬅️ hides the legend
  )

dev.off()


cluster2_data <- read.csv("Cluster_2_up.csv", header = TRUE, stringsAsFactors = FALSE)

# --- Step 2: Extract the list of genes/proteins ---
cluster_proteins <- cluster2_data$Gene

# Remove any missing values or duplicates
cluster_proteins <- unique(na.omit(cluster_proteins))

# Check what the list looks like
print(cluster_proteins)

gost_results <- gost(
  query = cluster_proteins,
  organism = "hsapiens",                # Change if not human
  sources = c("REAC", "KEGG", "GO:BP"), # Biological Process, KEGG, Reactome
  correction_method = "fdr",
  user_threshold = 0.05
)
# View results
head(gost_results$result)

# Make a copy
gost_results_flat <- gost_results$result

# Identify which columns are lists
list_cols <- sapply(gost_results_flat, is.list)

# Flatten every list column into comma-separated strings
gost_results_flat[list_cols] <- lapply(
  gost_results_flat[list_cols],
  function(col) sapply(col, function(x) paste(x, collapse = ","))
)

# Write to CSV safely
write.csv(gost_results_flat, "cluster2_up_gprofiler_enrichment.csv", row.names = FALSE)

# --- 1️⃣ Separate g:Profiler enrichment results ---
go_bp    <- gost_results_flat %>% filter(source == "GO:BP")
kegg     <- gost_results_flat %>% filter(source == "KEGG")
reactome <- gost_results_flat %>% filter(source == "REAC")

# Convert to TERM2GENE format for enricher()
msig_h_list <- msig_h %>%
  select(gs_name, gene_symbol) %>%
  as.data.frame()

# Run enrichment on your upregulated genes
ego_hallmark <- enricher(
  gene          = cluster_proteins,   # your input vector of gene symbols
  TERM2GENE     = msig_h_list,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Extract results
hallmark_df <- as.data.frame(ego_hallmark) %>%
  mutate(source = "MSigDB Hallmark")

# --- 3️⃣ Select top 3 from each source ---
go_top <- go_bp %>%
  arrange(p_value) %>%
  slice_head(n = 2) %>%
  mutate(source = "GO:BP")

kegg_top <- kegg %>%
  arrange(p_value) %>%
  slice_head(n = 2) %>%
  mutate(source = "KEGG")

reactome_top <- reactome %>%
  arrange(p_value) %>%
  slice_head(n = 2) %>%
  mutate(source = "REAC")

hallmark_top <- hallmark_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 2)

# --- 4️⃣ Combine all sources ---
combined_top <- bind_rows(go_top, kegg_top, reactome_top, hallmark_top) %>%
  mutate(
    term_label = case_when(
      !is.na(term_name) ~ paste0(term_name, " (", source, ")"),
      !is.na(Description) ~ paste0(Description, " (", source, ")"),
      TRUE ~ NA_character_
    ),
    adj_p = ifelse(!is.na(p_value), p_value, p.adjust),
    neg_log10_p = -log10(adj_p)
  ) %>%
  arrange(adj_p)

# --- 5️⃣ Define colors ---
source_colors <- c(
  "GO:BP"            = "darkred",  # firebrick (dark, classic red)
  "KEGG"             = "indianred",  # cinnabar
  "REAC"             = "#FF7F50",  # coral
  "MSigDB Hallmark"  = "red"   # brownish red
)

# --- 6️⃣ Plot unified bar chart ---
png(file="Cluster2_up.png", height=6, width=20, units="cm", res=600)

ggplot(combined_top,
       aes(x = reorder(term_label, neg_log10_p),
           y = neg_log10_p,
           fill = source)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = source_colors) +
  labs(
    x = NULL,
    y = expression(-log[10](adjusted~p)),
    title = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    legend.position = "none"    # ⬅️ hides the legend
  )

dev.off()


#--------pathway analyses and enrichment annotation: principle------------

cluster1_data <- read.csv("Cluster_1_down.csv",
                          header = TRUE, 
                          stringsAsFactors = FALSE) # insert clusters identified in Cytoscape

cluster_proteins <- cluster1_data$Gene

cluster_proteins <- unique(na.omit(cluster_proteins))
print(cluster_proteins)

gost_results <- gost(
  query = cluster_proteins,
  organism = "hsapiens",              
  sources = c("REAC", "KEGG", "GO:BP"), 
  correction_method = "fdr",
  user_threshold = 0.05
)
head(gost_results$result)

gost_results_flat <- gost_results$result

list_cols <- sapply(gost_results_flat, is.list)
gost_results_flat[list_cols] <- lapply(
  gost_results_flat[list_cols],
  function(col) sapply(col, function(x) paste(x, collapse = ","))
)

go_bp    <- gost_results_flat %>% filter(source == "GO:BP")
kegg     <- gost_results_flat %>% filter(source == "KEGG")
reactome <- gost_results_flat %>% filter(source == "REAC")

msig_h_list <- msig_h %>%
  select(gs_name, gene_symbol) %>%
  as.data.frame()

ego_hallmark <- enricher(
  gene          = cluster_proteins, 
  TERM2GENE     = msig_h_list,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

hallmark_df <- as.data.frame(ego_hallmark) %>%
  mutate(source = "MSigDB Hallmark")

#Select top 3 from each source
go_top <- go_bp %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  mutate(source = "GO:BP")

kegg_top <- kegg %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  mutate(source = "KEGG")

reactome_top <- reactome %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  mutate(source = "REAC")

hallmark_top <- hallmark_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)

combined_top <- bind_rows(go_top, kegg_top, reactome_top, hallmark_top) %>%
  mutate(
    term_label = case_when(
      !is.na(term_name) ~ paste0(term_name, " (", source, ")"),
      !is.na(Description) ~ paste0(Description, " (", source, ")"),
      TRUE ~ NA_character_
    ),
    adj_p = ifelse(!is.na(p_value), p_value, p.adjust),
    neg_log10_p = -log10(adj_p)
  ) %>%
  arrange(adj_p)

source_colors <- c(
  "GO:BP"            = "navy",  # firebrick (dark, classic red)
  "KEGG"             = "steelblue",  # cinnabar
  "REAC"             = "skyblue",  # coral
  "MSigDB Hallmark"  = "blue"   # brownish red
)

png(file="Cluster1_down.png", height=6, width=20, units="cm", res=600)

ggplot(combined_top,
       aes(x = reorder(term_label, neg_log10_p),
           y = neg_log10_p,
           fill = source)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = source_colors) +
  labs(
    x = NULL,
    y = expression(-log[10](adjusted~p)),
    title = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    legend.position = "none")

dev.off()


