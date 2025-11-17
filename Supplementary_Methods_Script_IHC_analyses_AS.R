
#---------Libraries-----------------

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


#--------------Set working directory-------------------------

setwd("D:/SARCOMA_RESEARCH/CSF1R_Angiosarcoma_Köln_Fortune/IHC_analyses")


#---------------Read in clinicopathological data--------------
col_all <- read.csv("AS_data_09_25.csv", header = TRUE, row.names = 1, sep=";")
write.csv(col_all,file="AS_backup.csv")
col_all$follow_up_m[col_all$follow_up_m < 1] <- NA
col_all$DFS_m[col_all$DFS_m < 1] <- NA
col_all$age<-as.numeric(col_all$age)
col_all$sex<-as.factor(col_all$sex)
col_all$primary_size_cm<-as.numeric(col_all$primary_size_cm)
col_all$AS_subtype_simple<-as.factor(col_all$AS_subtype_simple)


#-----------DESCRIPTIVE ANALYSES---------------------------

# Define list of cells for Immune Score Calculation
cell_list_p<- c(
  "p_T_cells_CD8", 
  "p_T_reg_cells",
  "p_T_cells_other",
  "p_M2_macrophages", 
  "p_M2_CSF1R_macrophages",
  "p_M0_M1_macrophages",
  "p_B_cells",
  "p_NK_cells",
  "p_Plasma_cells",
  "p_Mast_cells",
  "p_Neutrophils"
)
cell_list_r<- c(
  "r_T_cells_CD8", 
  "r_T_reg_cells",
  "r_T_cells_other",
  "r_M2_macrophages", 
  "r_M2_CSF1R_macrophages",
  "r_M0_M1_macrophages",
  "r_B_cells",
  "r_NK_cells",
  "r_Plasma_cells",
  "r_Mast_cells",
  "r_Neutrophils"
)

cell_list_m<- c(
  "m_T_cells_CD8", 
  "m_T_reg_cells",
  "m_T_cells_other",
  "m_M2_macrophages", 
  "m_M2_CSF1R_macrophages",
  "m_M0_M1_macrophages",
  "m_B_cells",
  "m_NK_cells",
  "m_Plasma_cells",
  "m_Mast_cells",
  "m_Neutrophils"
)

col_all <- col_all[complete.cases(col_all[, cell_list_p]), ]
col_all$p_ImmuneScore<- rowSums(col_all[,cell_list_p], na.rm = FALSE)
col_all$r_ImmuneScore<- rowSums(col_all[,cell_list_r], na.rm = FALSE)
col_all$m_ImmuneScore<- rowSums(col_all[,cell_list_m], na.rm = FALSE)


#-------- Supplementary Figure 1A: Draw heatmap--------------

list_heatmap <- c(
  "p_B_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_T_cells_CD8_proportions",
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions",
  "p_NK_cells_proportions",
  "p_Neutrophils_proportions",
  "p_Mast_cells_proportions"
)

cp_heatmap <- dplyr::select(col_all, list_heatmap)
cp_heatmap <- t(cp_heatmap)
cp_heatmap <- data.frame(cp_heatmap)
cp_heatmap_sort <- cp_heatmap[, intersect(rownames(col_all), colnames(cp_heatmap))]
cp_heatmap <- cp_heatmap[, colnames(cp_heatmap_sort)]

desired_order <-c(
  "p_B_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_T_cells_CD8_proportions",
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_Neutrophils_proportions",
  "p_Mast_cells_proportions",
  "p_NK_cells_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions"
)

desired_order <- intersect(desired_order, rownames(cp_heatmap))
cp_heatmap <- cp_heatmap[desired_order, ]
cp_heatmap_matrix <- as.matrix(cp_heatmap)

colsplit <- data.frame(Subtype = col_all$AS_subtype_simple)
rownames(colsplit) <- rownames(col_all)

levels_to_colors <- list(Subtype = c(
  "1" = "#F4A582", #1 = primary cutaneous AS
  "2"="#67001F", #2 = secondary AS
  "3"="#E69F00", #3 = visceral AS
  "4"="#666666" #4 = AS, other
))

ha <- HeatmapAnnotation(df = colsplit, col = levels_to_colors, na_col = "white")

bottom_annot_df <- data.frame(
  ImmuneScore = col_all$p_ImmuneScore
)

rownames(bottom_annot_df) <- rownames(col_all)
immune_range   <- range(col_all$p_ImmuneScore, na.rm = TRUE)
score_col_fun <- list(
  ImmuneScore   = colorRamp2(seq(immune_range[1], immune_range[2], length.out = 4),
                             c("black", "white", "lightgoldenrod", "darkorange")

))

bottom_annot <- HeatmapAnnotation(df = bottom_annot_df, col = score_col_fun, na_col = "white")

png(file="Supplementary_Figure1A_heatmap_.png", width=15, height=9, units="cm", res=1200)
heat <- Heatmap(cp_heatmap_matrix,
                top_annotation = ha,
                column_split = colsplit,
                name = "Relative Proportion",
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                column_title = "Sample",
                row_title = NULL,
                show_row_names = FALSE,
                show_column_names = FALSE,
                row_names_side = "left",
                col = colorRamp2(c(0, 0.02, 0.3, 0.5, 0.8),
                                 c("navy", "steelblue", "yellow", "red", "darkred")),
                bottom_annotation = bottom_annot)
print(heat)
dev.off()



#-------------Supplementary Figure 1B Stacked bar plot----------------

immune_cell_cols <-c(
  "p_T_cells_CD8_proportions", 
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_B_cells_proportions",
  "p_NK_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_Mast_cells_proportions",
  "p_Neutrophils_proportions"
)

df_long <- col_all %>%
  dplyr::select(all_of(immune_cell_cols), AS_subtype_simple) %>%  
  pivot_longer(cols = all_of(immune_cell_cols), 
               names_to = "Immune_Cell", 
               values_to = "Proportion")

desired_order<-c("p_B_cells_proportions",
                 "p_Plasma_cells_proportions",
                 "p_T_cells_CD8_proportions",
                 "p_T_reg_cells_proportions",
                 "p_T_cells_other_proportions",
                 "p_Neutrophils_proportions",
                 "p_Mast_cells_proportions",
                 "p_NK_cells_proportions",
                 "p_M0_M1_macrophages_proportions",
                 "p_M2_macrophages_proportions", 
                 "p_M2_CSF1R_macrophages_proportions")
immune_levels <- levels(factor(df_long$Immune_Cell)) 

df_long$Immune_Cell <- factor(df_long$Immune_Cell, levels = desired_order)

palette_22 <- c(
  "navy", 
  "mediumorchid2", 
  "red", 
  "#67001F", 
  "salmon", 
  "pink1", 
  "gold", 
  "cyan", 
  "dodgerblue", 
  "darkmagenta", 
  "black"  
)


df_long$Immune_Cell <- factor(df_long$Immune_Cell, levels = desired_order)

png(file="Supplementary_Figure_1B_stackedbar.png", width=11, height=7, units="cm", res=600)
ggplot(df_long, aes(x = AS_subtype_simple, y = Proportion, fill = Immune_Cell)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = NULL,
       x = NULL,
       y = NULL,
       fill = "Immune Cell Type") +
  scale_fill_manual(values = palette_22)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw()
dev.off()


#-------------compare immune cell parameters across AS subtypes-----------
cell_columns <- c("p_T_cells_CD8", 
  "p_T_reg_cells",
  "p_T_cells_other",
  "p_M2_macrophages", 
  "p_M2_CSF1R_macrophages",
  "p_M0_M1_macrophages",
  "p_B_cells",
  "p_NK_cells",
  "p_Plasma_cells",
  "p_Mast_cells",
  "p_Neutrophils",
  "p_ImmuneScore",
  "p_T_cells_CD8_proportions", 
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_B_cells_proportions",
  "p_NK_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_Mast_cells_proportions",
  "p_Neutrophils_proportions"
  )

kruskal_results <- data.frame(
  Cell_Type = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (cell in cell_columns) {
  test_result <- kruskal.test(
    formula = as.formula(paste(cell, "~ AS_subtype_simple")),
    data = col_all
  )
  
  kruskal_results <- rbind(kruskal_results, data.frame(
    Cell_Type = cell,
    p_value = test_result$p.value
  ))
}

kruskal_results$fdr_p_value <- p.adjust(kruskal_results$p_value, method = "fdr")
kruskal_results$bonferroni_p_value <- p.adjust(kruskal_results$p_value, method = "bonferroni")

print(kruskal_results)

immune_summary <- col_all %>%
  group_by(AS_subtype_simple) %>%
  summarise(
    count = n(),
    mean = mean(p_Neutrophils, na.rm = TRUE),
    median = median(p_Neutrophils, na.rm = TRUE),
    sd = sd(p_Neutrophils, na.rm = TRUE),
    min = min(p_Neutrophils, na.rm = TRUE),
    max = max(p_Neutrophils, na.rm = TRUE)
  )
print(immune_summary)


#----------Supplementary Figure 1C: Plot Neutrophils------------------
library(ggbreak)
my_colors <- c("1" = "#F4A582",
               "2"="#67001F",
               "3"="#E69F00",
               "4"="#666666")
png(file="S1_Neutrophil_counts.png", height=10, width=14, units="cm", res=600)
ggplot(col_all, aes(x = factor(AS_subtype_simple), y = p_Neutrophils)) +
  geom_boxplot(aes(fill = factor(AS_subtype_simple)), width = 0.4, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = "black") +
  scale_y_break(c(1000, 1300), scales = c(0.3, 0.1)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # major grid
    panel.grid.minor = element_line(color = "gray90", size = 0.3),  # minor grid
    axis.text = element_text(color = "grey40", size = 12),
    axis.ticks = element_line(color = "grey40"),
    axis.line.x = element_line(color = "grey40"),
    axis.line.y = element_line(color = "grey40")
  )
dev.off()
summary(col_all$p_ImmuneScore)


#-------------- Supplementary Figure 1D & 1E: Plot NK cells and Immune Score------------------

png(file="S1_ImmuneScore_cells_counts.png", height=10, width=14, units="cm", res=600)
ggplot(col_all, aes(x = factor(AS_subtype_simple), y = p_ImmuneScore)) +
  geom_boxplot(aes(fill = factor(AS_subtype_simple)), width = 0.4, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = "black") +
  scale_y_break(c(7500, 8000), scales = c(0.1, 0.1)) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # major grid
    panel.grid.minor = element_line(color = "gray90", size = 0.3),  # minor grid
    axis.text = element_text(color = "grey40", size = 12),
    axis.ticks = element_line(color = "grey40"),
    axis.line.x = element_line(color = "grey40"),
    axis.line.y = element_line(color = "grey40")
  )
dev.off()


#-------------Supplementary Figure 1F: Correlation of immune cells in AS---------

selected_columns <- col_all[, c("p_T_cells_CD8", 
  "p_T_reg_cells",
  "p_T_cells_other",
  "p_M2_macrophages", 
  "p_M2_CSF1R_macrophages",
  "p_M0_M1_macrophages",
  "p_B_cells",
  "p_NK_cells",
  "p_Plasma_cells",
  "p_Mast_cells",
  "p_Neutrophils")]

cor_test <- rcorr(as.matrix(selected_columns))
cor_matrix <- cor_test$r
p_matrix <- cor_test$P

melted_cor_matrix <- melt(cor_matrix)
melted_p_matrix <- melt(p_matrix)
melted_cor_matrix$p_value <- melted_p_matrix$value

melted_cor_matrix$Var1_num <- as.numeric(melted_cor_matrix$Var1)
melted_cor_matrix$Var2_num <- as.numeric(melted_cor_matrix$Var2)

lower_tri <- subset(melted_cor_matrix, Var1_num >= Var2_num)

var_levels <- levels(melted_cor_matrix$Var1)
lower_tri$Var1 <- factor(lower_tri$Var1, levels = rev(var_levels))
lower_tri$Var2 <- factor(lower_tri$Var2, levels = var_levels)

lower_tri$fontface <- ifelse(!is.na(lower_tri$p_value) & lower_tri$p_value < 0.05, "bold", "plain")

png(file="correlation_plot_immune cells.png", height=15, width=15, units="cm", res=1200)
ggplot(data = lower_tri, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value), fontface = fontface),
            color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  theme_minimal() +
  labs(x = "", y = "", title = NULL) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # rotate x labels
  )

dev.off()


#------Supplementary Figure 1G: correlation of immune cell proportions in AS----------

selected_columns <- col_all[, c("p_T_cells_CD8_proportions", 
                                "p_T_reg_cells_proportions",
                                "p_T_cells_other_proportions",
                                "p_M2_macrophages_proportions", 
                                "p_M2_CSF1R_macrophages_proportions",
                                "p_M0_M1_macrophages_proportions",
                                "p_B_cells_proportions",
                                "p_NK_cells_proportions",
                                "p_Plasma_cells_proportions",
                                "p_Mast_cells_proportions",
                                "p_Neutrophils_proportions")]

cor_test <- rcorr(as.matrix(selected_columns))
cor_matrix <- cor_test$r
p_matrix <- cor_test$P

melted_cor_matrix <- melt(cor_matrix)
melted_p_matrix <- melt(p_matrix)
melted_cor_matrix$p_value <- melted_p_matrix$value

melted_cor_matrix$Var1_num <- as.numeric(melted_cor_matrix$Var1)
melted_cor_matrix$Var2_num <- as.numeric(melted_cor_matrix$Var2)

lower_tri <- subset(melted_cor_matrix, Var1_num >= Var2_num)

var_levels <- levels(melted_cor_matrix$Var1)
lower_tri$Var1 <- factor(lower_tri$Var1, levels = rev(var_levels))
lower_tri$Var2 <- factor(lower_tri$Var2, levels = var_levels)

lower_tri$fontface <- ifelse(!is.na(lower_tri$p_value) & lower_tri$p_value < 0.05, "bold", "plain")

png(file="correlation_plot_immune cell_proportions.png", height=20, width=20, units="cm", res=1200)
ggplot(data = lower_tri, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value), fontface = fontface),
            color = "black", size = 3) +
  scale_fill_gradient2(low = "navy", high = "darkorange3", mid = "white", midpoint = 0,
                       limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  theme_minimal() +
  labs(x = "", y = "", title = NULL) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # rotate x labels
  )

dev.off()



#--------------Desriptive statistics-----------------

test_res <- kruskal.test(age ~ AS_subtype_simple, data = col_all)
print(test_res)
test_res <- kruskal.test(primary_size_cm ~ AS_subtype_simple, data = col_all)
print(test_res)
chisq.test(col_all$AS_subtype_simple, col_all$sex)
chisq.test(col_all$AS_subtype_simple, col_all$AS_epithelioid)
chisq.test(col_all$AS_subtype_simple, col_all$localization)
chisq.test(col_all$AS_subtype_simple, col_all$AS_multifocal)
chisq.test(col_all$AS_subtype_simple, col_all$R)

col_all$M_stage <- ifelse(col_all$metastasis_yes_no == 1 & col_all$metastasis_synchrone_metachrone == 1, 1,
                          ifelse(is.na(col_all$metastasis_yes_no), NA, 0))
chisq.test(col_all$AS_subtype_simple, col_all$M_stage)

#Adjust p with Bonferroni corrections
p.adjust(c(0.004810525, #age
           0.5490727, #tumor_size
           4.987e-06, #sex
           0.1512, #epithelioid_morphology
           2.2e-16, #localization
           0.008377, #multifocal
           0.001135, #R status
           0.007904), #M stage
         method="bonferroni")


#---------SURVIVAL ANALYSES for macrophages ------------

#set  dichotomizations
col_surv<-col_all

col_surv$M0_M1_prop_median <- as.factor(ifelse(col_surv$p_M0_M1_macrophages_proportions >
                                                median(col_surv$p_M0_M1_macrophages_proportions, na.rm = TRUE), 1, 2))
col_surv$M2_CD163_prop_median <- as.factor(ifelse(col_surv$p_M2_macrophages_proportions >
                                                   median(col_surv$p_M2_macrophages_proportions, na.rm = TRUE), 1, 2))
col_surv$M2_CD163_CSF1R_prop_median <- as.factor(ifelse(col_surv$p_M2_CSF1R_macrophages_proportions >
                                                         median(col_surv$p_M2_CSF1R_macrophages_proportions, na.rm = TRUE), 1, 2))
col_surv$M2_all_prop_median <- as.factor(ifelse(col_surv$p_M2_macrophages_all_proportions >
                                                 median(col_surv$p_M2_macrophages_all_proportions , na.rm = TRUE), 1, 2))

median(col_surv$p_M0_M1_macrophages_proportions)
median(col_surv$p_M2_macrophages_proportions)
median(col_surv$p_M2_CSF1R_macrophages_proportions)

#Overall_Survival
surv_object <- Surv(col_surv$follow_up_m, col_surv$vital_status)
fit_M0_M1_prop_median <- survfit(surv_object ~ M0_M1_prop_median, data = col_surv)

png(file="Fig1_M0_M1_OS.png", height=6, width=8, units="cm", res=600)
ggsurvplot(fit_M0_M1_prop_median, data = col_surv, pval = F, risk.table = F,
           legend = "none", 
           title = NULL, 
           palette = c("dodgerblue","grey40"),
           xlab = "Time (months)", ylab= "OS Probability")
dev.off()

png(file="Fig1_M0_M1_OS_table.png", height=14, width=12, units="cm", res=600)
ggsurvplot(fit_M0_M1_prop_median, data = col_surv, pval = T, risk.table = T,
           legend = "none", 
           title = NULL, 
           palette = c( "dodgerblue","grey40"),
           xlab = "Time (months)", ylab= "OS Probability")
dev.off()

immune_summary <- col_surv %>%
  group_by(M0_M1_prop_median) %>%
  summarise(
    count = n(),
    mean = mean(follow_up_m, na.rm = TRUE),
    median = median(follow_up_m, na.rm = TRUE),
    sd = sd(follow_up_m, na.rm = TRUE),
    min = min(follow_up_m, na.rm = TRUE),
    max = max(follow_up_m, na.rm = TRUE)
  )
print(immune_summary)

surv_object <- Surv(col_surv$follow_up_m, col_surv$vital_status)
fit_M2_CD163_prop_median <- survfit(surv_object ~ M2_CD163_prop_median, data = col_surv)

png(file="Fig1_CD163_OS.png", height=6, width=8, units="cm", res=600)
ggsurvplot(fit_M2_CD163_prop_median, data = col_surv, pval = F, risk.table = F,
           legend = "none", 
           title = NULL, 
           palette = c("darkmagenta","grey40"), xlab = "Time (months)", ylab= "OS Probability")
dev.off()

png(file="Fig1_CD163_OS_table.png", height=14, width=12, units="cm", res=600)
ggsurvplot(fit_M2_CD163_prop_median, data = col_surv, pval = T, risk.table = T,
           legend = "none", 
           title = NULL, 
           palette = c( "darkmagenta","grey40"), xlab = "Time (months)", ylab= "OS Probability")
dev.off()

immune_summary <- col_surv %>%
  group_by(M2_CD163_prop_median) %>%
  summarise(
    count = n(),
    mean = mean(follow_up_m, na.rm = TRUE),
    median = median(follow_up_m, na.rm = TRUE),
    sd = sd(follow_up_m, na.rm = TRUE),
    min = min(follow_up_m, na.rm = TRUE),
    max = max(follow_up_m, na.rm = TRUE)
  )
print(immune_summary)

surv_object <- Surv(col_surv$follow_up_m, col_surv$vital_status)
fit_M2_CD163_CSF1R_prop_median <- survfit(surv_object ~ M2_CD163_CSF1R_prop_median, data = col_surv)
png(file="Fig1_CD163_CSF1R_OS.png", height=6, width=8, units="cm", res=600)
ggsurvplot(fit_M2_CD163_CSF1R_prop_median, data = col_surv, pval = F, risk.table = F,
           legend = "none", 
           title = NULL, 
           palette = c("black","grey60"), xlab = "Time (months)", ylab= "OS Probability")
dev.off()

png(file="Fig1_CD163_CSF1R_OS_table.png", height=14, width=15, units="cm", res=600)
ggsurvplot(fit_M2_CD163_CSF1R_prop_median, data = col_surv, pval = T, risk.table = T,
           legend = "none", 
           title = NULL, 
           palette = c("black","grey60"), xlab = "Time (months)", ylab= "OS Probability")
dev.off()

write.csv(col_surv, file="col_surv_median_CSF1R.csv")
immune_summary <- col_surv %>%
  group_by(M2_CD163_CSF1R_prop_median) %>%
  summarise(
    count = n(),
    mean = mean(follow_up_m, na.rm = TRUE),
    median = median(follow_up_m, na.rm = TRUE),
    sd = sd(follow_up_m, na.rm = TRUE),
    min = min(follow_up_m, na.rm = TRUE),
    max = max(follow_up_m, na.rm = TRUE)
  )
print(immune_summary)


#Disease-free survival
surv_object <- Surv(col_surv$DFS_m, col_surv$DFS_status)
fit_M0_M1_prop_median <- survfit(surv_object ~ M0_M1_prop_median, data = col_surv)
png(file="Fig1_M0_M1_DFS.png", height=6, width=8, units="cm", res=600)
ggsurvplot(fit_M0_M1_prop_median, data = col_surv, pval = F, risk.table = F,
           legend = "none", 
           title = NULL, 
           palette = c("dodgerblue","grey40"),
           xlab = "Time (months)", ylab="DFS Probability")
dev.off()

png(file="Fig1_M0_M1_DFS_table.png", height=14, width=12, units="cm", res=600)
ggsurvplot(fit_M0_M1_prop_median, data = col_surv, pval = T, risk.table = T,
           legend = "none", 
           title = NULL, 
           palette = c("dodgerblue","grey40"),
           xlab = "Time (months)", ylab="DFS Probability")
dev.off()

immune_summary <- col_surv %>%
  group_by(M0_M1_prop_median) %>%
  summarise(
    count = n(),
    mean = mean(DFS_m, na.rm = TRUE),
    median = median(DFS_m, na.rm = TRUE),
    sd = sd(DFS_m, na.rm = TRUE),
    min = min(DFS_m, na.rm = TRUE),
    max = max(DFS_m, na.rm = TRUE)
  )
print(immune_summary)

surv_object <- Surv(col_surv$DFS_m, col_surv$DFS_status)
fit_M2_CD163_prop_median <- survfit(surv_object ~ M2_CD163_prop_median, data = col_surv)

png(file="Fig1_CD163_DFS.png", height=6, width=8, units="cm", res=600)
ggsurvplot(fit_M2_CD163_prop_median, data = col_surv, pval = F,risk.table = F,
           legend = "none", 
           title = NULL, 
           palette = c("darkmagenta","grey40"), xlab = "Time (months)", ylab="DFS Probability")
dev.off()

png(file="Fig1_CD163_DFS_table.png", height=14, width=12, units="cm", res=600)
ggsurvplot(fit_M2_CD163_prop_median, data = col_surv, pval = T,risk.table = T,
           legend = "none", 
           title = NULL, 
           palette = c("darkmagenta","grey40"), xlab = "Time (months)",ylab="DFS Probability")
dev.off()

immune_summary <- col_surv %>%
  group_by(M2_CD163_prop_median) %>%
  summarise(
    count = n(),
    mean = mean(DFS_m, na.rm = TRUE),
    median = median(DFS_m, na.rm = TRUE),
    sd = sd(DFS_m, na.rm = TRUE),
    min = min(DFS_m, na.rm = TRUE),
    max = max(DFS_m, na.rm = TRUE)
  )
print(immune_summary)

surv_object <- Surv(col_surv$DFS_m, col_surv$DFS_status)
fit_M2_CD163_CSF1R_prop_median <- survfit(surv_object ~ M2_CD163_CSF1R_prop_median, data = col_surv)

png(file="Fig1_CSF1R_DFS.png", height=6, width=8, units="cm", res=600)
ggsurvplot(fit_M2_CD163_CSF1R_prop_median, data = col_surv, pval = F,risk.table = F,
           legend = "none", 
           title = NULL, 
           palette = c("black","grey46"), xlab = "Time (months)", ylab="DFS Probability")
dev.off()

png(file="Fig1_CSF1R_DFS_table.png", height=14, width=15, units="cm", res=600)
ggsurvplot(fit_M2_CD163_CSF1R_prop_median, data = col_surv, pval = T,risk.table = T,
           legend = "none", 
           title = NULL, 
           palette = c("black","grey46"), xlab = "Time (months)",ylab="DFS Probability")
dev.off()


immune_summary <- col_surv %>%
  group_by(M2_CD163_CSF1R_prop_median) %>%
  summarise(
    count = n(),
    mean = mean(DFS_m, na.rm = TRUE),
    median = median(DFS_m, na.rm = TRUE),
    sd = sd(DFS_m, na.rm = TRUE),
    min = min(DFS_m, na.rm = TRUE),
    max = max(DFS_m, na.rm = TRUE)
  )
print(immune_summary)



#---------K-means clustering of immune cells -------------------

prop_list <- c(
  "p_T_cells_CD8_proportions", 
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_B_cells_proportions",
  "p_NK_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_Mast_cells_proportions",
  "p_Neutrophils_proportions",
  "p_ImmuneScore"
)

col_all_cluster <- col_all[complete.cases(col_all[, prop_list]), ]
col_all_cluster<-subset(col_all_cluster, AS_subtype_simple != 4)
cell_data <- col_all_cluster[, prop_list]

df_scaled <- scale(cell_data)

k_range <- 2:6
survival_pvalues <- c()
for (k in k_range) {
  set.seed(42)
  kmeans_model <- kmeans(df_scaled, centers = k, nstart = 25)
  col_all_cluster$Cluster <- as.factor(kmeans_model$cluster)
  
  surv_object <- Surv(col_all_cluster$follow_up_m, col_all_cluster$vital_status)
  fit <- survfit(surv_object ~ Cluster, data = col_all_cluster)

  surv_test <- survdiff(surv_object ~ Cluster, data = col_all_cluster)
  p_value <- 1 - pchisq(surv_test$chisq, length(surv_test$n) - 1)
  
  survival_pvalues <- c(survival_pvalues, p_value)
}

best_k <- k_range[which.min(survival_pvalues)]
cat("Optimal number of clusters based on survival analysis:", best_k, "\n")

set.seed(321)
final_kmeans <- kmeans(df_scaled, centers = best_k, nstart = 25)
col_all_cluster$Cluster <- as.factor(final_kmeans$cluster)

new_labels <- c(`3` = 1, `1` = 2, `2` = 3)
col_all_cluster$Cluster <- as.factor(new_labels[as.character(col_all_cluster$Cluster)])

table(col_all_cluster$Cluster, col_all_cluster$histology)

table(col_all_cluster$Cluster)

#Figure 2A: PCA plot
PCA_cols <- c(
  "p_B_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_T_cells_CD8_proportions",
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions",
  "p_NK_cells_proportions",
  "p_Neutrophils_proportions",
  "p_Mast_cells_proportions",
  "p_ImmuneScore"
)

immune_data <- col_all_cluster[, PCA_cols]

pca_result <- prcomp(immune_data, scale. = TRUE)
pca_data   <- as.data.frame(pca_result$x)
pca_data$Cluster <- as.factor(col_all_cluster$Cluster)

pc_var      <- pca_result$sdev^2
pc_var_exp  <- pc_var / sum(pc_var) * 100

custom_colors <- c("black", "darkmagenta", "darkorange")

png("Figure_2A_PCA_plot.png", height = 7, width=9, units="cm", res=600)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 1.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, size = 0.8) +
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(pc_var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(pc_var_exp[2], 1), "%)"),
    color = "Cluster"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_grey()
dev.off()

#Figure 2B: Heatmap
list_heatmap <- c(
  "p_B_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_T_cells_CD8_proportions",
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions",
  "p_NK_cells_proportions",
  "p_Neutrophils_proportions",
  "p_Mast_cells_proportions"
)

cp_heatmap <- dplyr::select(col_all_cluster, list_heatmap)
cp_heatmap <- t(cp_heatmap)
cp_heatmap <- data.frame(cp_heatmap)
cp_heatmap_sort <- cp_heatmap[, intersect(rownames(col_all_cluster), colnames(cp_heatmap))]
cp_heatmap <- cp_heatmap[, colnames(cp_heatmap_sort)]
desired_order <-c(
  "p_B_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_T_cells_CD8_proportions",
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_Neutrophils_proportions",
  "p_Mast_cells_proportions",
  "p_NK_cells_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions"
)

desired_order <- intersect(desired_order, rownames(cp_heatmap))
cp_heatmap <- cp_heatmap[desired_order, ]
cp_heatmap_matrix <- as.matrix(cp_heatmap)
colsplit <- data.frame(Cluster = col_all_cluster$Cluster)
rownames(colsplit) <- rownames(col_all_cluster)

levels_to_colors <- list(Cluster = c(
  "1" = "black",
  "2" = "darkmagenta",
  "3" = "darkorange"
))

ha <- HeatmapAnnotation(df = colsplit, col = levels_to_colors, na_col = "white")
table(col_all_cluster$AS_subtype_simple)
bottom_annot_df <- data.frame(
  ImmuneScore = col_all_cluster$p_ImmuneScore,
  Histology = factor(col_all_cluster$AS_subtype_simple, levels = sort(unique(col_all_cluster$AS_subtype_simple)))
)
rownames(bottom_annot_df) <- rownames(col_all_cluster)

immune_range   <- range(col_all_cluster$p_ImmuneScore, na.rm = TRUE)

score_col_fun <- list(
  ImmuneScore   = colorRamp2(seq(immune_range[1], immune_range[2], length.out = 4), c("black", "white", "lightgoldenrod", "darkorange")),
  Histology     = c("1" = "#F4A582",
                    "2"="#67001F",
                    "3"="#E69F00",
                    "4"="#666666")
)

bottom_annot <- HeatmapAnnotation(df = bottom_annot_df, col = score_col_fun, na_col = "white")

png(file="Figure2B_heatmap_2.png", width=15, height=15, units="cm", res=1200)
heat <- Heatmap(cp_heatmap_matrix,
                top_annotation = ha,
                column_split = colsplit,
                name = "Relative Proportion",
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                column_title = "Sample",
                row_title = NULL,
                show_row_names = FALSE,
                show_column_names = FALSE,
                row_names_side = "left",
                col = colorRamp2(c(0, 0.02, 0.3, 0.5, 0.8),
                                 c("navy", "steelblue", "yellow", "red", "darkred")),
                bottom_annotation = bottom_annot)
print(heat)
dev.off()

# Figure 2C Stacked bar plot

immune_cell_cols <-c(
  "p_T_cells_CD8_proportions", 
  "p_T_reg_cells_proportions",
  "p_T_cells_other_proportions",
  "p_M2_macrophages_proportions", 
  "p_M2_CSF1R_macrophages_proportions",
  "p_M0_M1_macrophages_proportions",
  "p_B_cells_proportions",
  "p_NK_cells_proportions",
  "p_Plasma_cells_proportions",
  "p_Mast_cells_proportions",
  "p_Neutrophils_proportions"
)

df_long <- col_all_cluster %>%
  dplyr::select(all_of(immune_cell_cols), Cluster) %>%  
  pivot_longer(cols = all_of(immune_cell_cols), 
               names_to = "Immune_Cell", 
               values_to = "Proportion")

desired_order<-c("p_B_cells_proportions",
                 "p_Plasma_cells_proportions",
                 "p_T_cells_CD8_proportions",
                  "p_T_reg_cells_proportions",
                   "p_T_cells_other_proportions",
                   "p_Neutrophils_proportions",
                   "p_Mast_cells_proportions",
                   "p_NK_cells_proportions",
                   "p_M0_M1_macrophages_proportions",
                   "p_M2_macrophages_proportions", 
                   "p_M2_CSF1R_macrophages_proportions")
immune_levels <- levels(factor(df_long$Immune_Cell)) 

df_long$Immune_Cell <- factor(df_long$Immune_Cell, levels = desired_order)

palette_22 <- c(
  "navy", 
  "mediumorchid2", 
  "red", 
  "#67001F", 
  "salmon", 
  "pink1", 
  "gold", 
  "cyan", 
  "dodgerblue", 
  "darkmagenta", 
  "black"  
)


df_long$Immune_Cell <- factor(df_long$Immune_Cell, levels = desired_order)

png(file="Figure_2C_stackedbar.png", width=11, height=7, units="cm", res=600)
ggplot(df_long, aes(x = Cluster, y = Proportion, fill = Immune_Cell)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = NULL,
       x = NULL,
       y = NULL,
       fill = "Immune Cell Type") +
  scale_fill_manual(values = palette_22)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw()
dev.off()


#--------------Survival analyses for Clusters----------------------

# Overall survival
surv_object <- Surv(col_all_cluster$follow_up_m, col_all_cluster$vital_status)
fit_final <- survfit(surv_object ~ Cluster, data = col_all_cluster)

png(file="Fig2_Cluster_OS.png", width=10, height = 6, units="cm", res=600)
ggsurvplot(fit_final, data = col_all_cluster, pval = F, risk.table = F,
           legend = "none", 
           title = NULL, 
           palette = c("black", "darkmagenta","darkorange"), 
           xlab = "Time (months)", ylab= "OS Probability")
dev.off()

png(file="Fig2_Cluster_OS_table.png", width=12, height = 17, units="cm", res=600)
ggsurvplot(fit_final, data = col_all_cluster, pval = T, risk.table = T,
           legend = "none", 
           title = NULL, 
           palette = c("black", "darkmagenta","darkorange"), 
           xlab = "Time (months)", ylab= "OS Probability")
dev.off()

immune_summary <- col_all_cluster %>%
  group_by(Cluster) %>%
  summarise(
    count = n(),
    mean = mean(follow_up_m, na.rm = TRUE),
    median = median(follow_up_m, na.rm = TRUE),
    sd = sd(follow_up_m, na.rm = TRUE),
    min = min(follow_up_m, na.rm = TRUE),
    max = max(follow_up_m, na.rm = TRUE)
  )
print(immune_summary)

# Disease-free survival
surv_object <- Surv(col_all_cluster$DFS_m, col_all_cluster$DFS_status)
fit_final <- survfit(surv_object ~ Cluster, data = col_all_cluster)

png(file="Fig2_Cluster_DFS.png", width=10, height = 6, units="cm", res=600)
ggsurvplot(fit_final, data = col_all_cluster, pval = F, risk.table = F,
           legend = "none", 
           title = NULL, 
           palette = c("black", "darkmagenta","darkorange"),
           xlab = "Time (months)", ylab= "DFS Probability")
dev.off()


png(file="Fig2_Cluster_DFS_table.png", width=12, height = 17, units="cm", res=600)
ggsurvplot(fit_final, data = col_all_cluster, pval = T, risk.table = T,
           legend = "none", 
           title = NULL, 
           palette = c("black", "darkmagenta","darkorange"), 
           xlab = "Time (months)", ylab= "DFS Probability")
dev.off()

immune_summary <- col_all_cluster %>%
  group_by(Cluster) %>%
  summarise(
    count = n(),
    mean = mean(DFS_m, na.rm = TRUE),
    median = median(DFS_m, na.rm = TRUE),
    sd = sd(DFS_m, na.rm = TRUE),
    min = min(DFS_m, na.rm = TRUE),
    max = max(DFS_m, na.rm = TRUE)
  )
print(immune_summary)



#-----------Correlations of cluster with histopathological parameters---------
test_res <- kruskal.test(age ~ Cluster, data = col_all_cluster)
print(test_res)
age <- col_all_cluster %>%
  group_by(Cluster) %>%
  summarise(
    count = n(),
    mean = mean(age, na.rm = TRUE),
    median = median(age, na.rm = TRUE),
    sd = sd(age, na.rm = TRUE),
    min = min(age, na.rm = TRUE),
    max = max(age, na.rm = TRUE)
  )
print(age)

test_res <- kruskal.test(primary_size_cm ~ Cluster, data = col_all_cluster)
print(test_res)
size <- col_all_cluster %>%
  group_by(Cluster) %>%
  summarise(
    count = n(),
    mean = mean(primary_size_cm, na.rm = TRUE),
    median = median(primary_size_cm, na.rm = TRUE),
    sd = sd(primary_size_cm, na.rm = TRUE),
    min = min(primary_size_cm, na.rm = TRUE),
    max = max(primary_size_cm, na.rm = TRUE)
  )
print(size)
chisq.test(col_all_cluster$Cluster, col_all_cluster$sex)
chisq.test(col_all_cluster$Cluster, col_all_cluster$AS_subtype_simple)
chisq.test(col_all_cluster$Cluster, col_all_cluster$AS_epithelioid)
chisq.test(col_all_cluster$Cluster, col_all_cluster$localization)
chisq.test(col_all_cluster$Cluster, col_all_cluster$AS_multifocal)
chisq.test(col_all_cluster$Cluster, col_all_cluster$R)
table(col_all_cluster$Cluster, col_all_cluster$R)
chisq.test(col_all_cluster$Cluster, col_all_cluster$M_stage)
table(col_all_cluster$Cluster, col_all_cluster$M_stage)
test_res <- kruskal.test(p_ImmuneScore ~ Cluster, data = col_all_cluster)
print(test_res)
immune <- col_all_cluster %>%
  group_by(Cluster) %>%
  summarise(
    count = n(),
    mean = mean(p_ImmuneScore, na.rm = TRUE),
    median = median(p_ImmuneScore, na.rm = TRUE),
    sd = sd(p_ImmuneScore, na.rm = TRUE),
    min = min(p_ImmuneScore, na.rm = TRUE),
    max = max(p_ImmuneScore, na.rm = TRUE)
  )
print(immune)
#Adjust p with Bonferroni corrections
p.adjust(c(0.3462, #age
           0.5794, #tumor_size
           0.2218, #sex
           0.3524, #subtype
           0.7848, #epithelioid_morphology
           0.06206, #localization
           0.6553, #multifocal
           0.04763, #R status
           0.04851,#M stage
           0.002214), #Immune Score
         method="bonferroni")

table(subset(col_all_cluster, Cluster==1)$M_stage)
table(subset(col_all_cluster, Cluster==1)$M_stage)/20

col_all_cluster$cluster_CSF1R_vs_other<-ifelse(col_all_cluster$Cluster ==1, 1, 0)



#-------------identify univariate survival factors--------------------
variables <- c("cluster_CSF1R_vs_other", 
               "age", 
               "sex", 
               "R", 
               "localization",
               "AS_epithelioid",
               "AS_subtype_simple",
               "M_stage",
               "primary_size_cm",
               "AS_multifocal",
               "p_ImmuneScore")
significant_vars <- c()

surv_object <- Surv(col_all_cluster$follow_up_m, col_all_cluster$vital_status)

for (var in variables) {
  formula_str <- as.formula(paste("surv_object ~", var))
  cox_model <- try(coxph(formula_str, data = col_all_cluster), silent = TRUE)
  if (!inherits(cox_model, "try-error")) {
    pval <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
    if (!is.na(pval) && pval <= 0.05) {
      significant_vars <- c(significant_vars, var)
    }
  }
}

cat("Signifikante Variablen (p <= 0.05):", paste(significant_vars, collapse = ", "), "\n")


#----------Multivariable analyses--------------
surv_object <- Surv(col_all_cluster$follow_up_months, col_all_cluster$vital_status)
cox_multi <- coxph(surv_object ~ cluster_CSF1R_vs_other+
                     sex+M_stage,
                   data = col_all_cluster
)
summary(cox_multi)

surv_object <- Surv(col_all_cluster$DFS_m, col_all_cluster$DFS_status)
cox_multi <- coxph(surv_object ~ cluster_CSF1R_vs_other+
                     sex+M_stage,
                   data = col_all_cluster
)
summary(cox_multi)


