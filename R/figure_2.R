# Figure 2: Immune and Methylation Profiling
# Figure 2A
# VP to add relevant scripts and processed data

# Figure 2B
# VP to add relevant scripts and processed data

# Figure 2C
# VP to add relevant scripts and processed data

# Figure 2D
# VP to add relevant scripts and processed data
# Figure 2E: ROC Analysis - ssGSEA Score Validation
library(pROC)

d1 <- read.delim('data/ValidationCohort_ssGSEA_scores.txt')

d1$response_binary <- ifelse(d1$group == "PD", 1, 0)
d1$PD <- as.numeric(d1$PD)

roc_obj <- roc(
  response = d1$response_binary,
  predictor = d1$PD,
  ci = TRUE,
  boot.n = 2000
)

plot(
  roc_obj,
  col = "red",
  lwd = 3,
  ci = TRUE,
  ci.type = "shape",
  ci.col = "#cce5ff",
  legacy.axes = TRUE,
  print.auc = TRUE,
  print.auc.cex = 1.5,
  print.auc.x = 0.6,
  print.auc.y = 0.4,
  grid = FALSE,
  auc.polygon = FALSE
)
# Figure 2F: Kaplan-Meier - PFS by Predicted Response Groups
library(survival)
library(survminer)

collapsed_data <- read.csv('data/baseline_data.csv')

collapsed_data <- collapsed_data[order(collapsed_data$score_max),]

collapsed_data$Category <- "Non-responder"
collapsed_data$Color <- "grey"

collapsed_data$Category[collapsed_data$solid_patient == "Solid 19"] <- "Complete Responder"
collapsed_data$Color[collapsed_data$solid_patient == "Solid 19"] <- "forestgreen"

collapsed_data$Category[collapsed_data$solid_patient %in% c("Solid 13", "Solid 27", "Solid 26")] <- "Partial Responder"
collapsed_data$Color[collapsed_data$solid_patient %in% c("Solid 13", "Solid 27", "Solid 26")] <- "goldenrod"

collapsed_data <- collapsed_data[order(collapsed_data$mean_score),]
collapsed_data$group <- ifelse(collapsed_data$median_score < 0.3, "<0.3", "≥0.3")

surv_object <- Surv(time = collapsed_data$pfs, event = rep(1, nrow(collapsed_data)))

fit <- survfit(surv_object ~ group, data = collapsed_data)

ggsurvplot(
  fit,
  data = collapsed_data,
  pval = F,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.title = "Group",
  legend.labs = c("Predicted non-responder", "Predicted responder"),
  xlab = "Time (months)",
  ylab = "Progression-Free Survival",
  palette = c("darkred", "darkblue")
)
# Figure 2G: Kaplan-Meier - PD Stratification in IDH-mutant Gliomas
library(survminer)
library(survival)
library(survivalROC)

survival <- read.csv('data/medip_idhm_OS.csv')
survival$overall_survival_years <- as.numeric(as.character(survival$overall_survival)) / 365

table(cut(survival$overall_survival_years, breaks = c(0, 2, 5, 10, 20)), 1 - survival$alive)
roc <- survivalROC(
  Stime = survival$overall_survival_years,
  status = 1 - survival$alive,
  marker = survival$PD,
  predict.time = 5,
  method = "KM"
)

youden_index <- which.max(roc$TP - roc$FP)
optimal_cutoff <- roc$cut.values[youden_index]
cat("Optimal PD cutoff (Youden's Index):", optimal_cutoff, "\n")

survival$PD_group <- ifelse(survival$PD > optimal_cutoff, "High PD", "Low PD")

surv_object <- Surv(time = survival$overall_survival_years, event = 1 - survival$alive)
fit <- survfit(surv_object ~ PD_group, data = survival)

plot_obj <- ggsurvplot(
  fit,
  data = survival,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.title = NULL,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  palette = c("#3E4989", "#9DB3DE"),
  break.time.by = 5,
  xlab = "Years",
  ylab = "Overall Survival Probability",
  surv.median.line = "none",
  legend.title = NULL,
  ggtheme = theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "top",
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
)

print(plot_obj)
# Figure 2H: CIBERSORT Immune Deconvolution Pie Charts
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

data <- read.csv("data/estimate_scores.csv")

cibersort_data <- data[grepl("_CIBERSORT$", data$cell_type), ]
cibersort_data$cell_type_clean <- gsub("_CIBERSORT$", "", cibersort_data$cell_type)

cibersort_data$cell_type_simplified <- cibersort_data$cell_type_clean
cibersort_data$cell_type_simplified <- gsub("T cell CD4\\+ memory resting", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell CD4\\+ memory activated", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell CD4\\+ naive", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell follicular helper", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell regulatory \\(Tregs\\)", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell gamma delta", "T cell CD8+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("B cell naive", "B cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("B cell memory", "B cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("B cell plasma", "B cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("NK cell resting", "NK cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("NK cell activated", "NK cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Macrophage M0", "Macrophage", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Macrophage M1", "Macrophage", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Macrophage M2", "Macrophage", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Myeloid dendritic cell resting", "mDC", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Myeloid dendritic cell activated", "mDC", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Mast cell activated", "Mast cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Mast cell resting", "Mast cell", cibersort_data$cell_type_simplified)

responder_data <- aggregate(Responders ~ cell_type_simplified, data = cibersort_data, sum)
non_responder_data <- aggregate(Non_Responders ~ cell_type_simplified, data = cibersort_data, sum)

responder_data <- responder_data[responder_data$Responders > 0, ]
non_responder_data <- non_responder_data[non_responder_data$Non_Responders > 0, ]

colors <- c(
  "T cell CD8+" = "#E31A1C",
  "T cell CD4+" = "#1F78B4",
  "Macrophage" = "#33A02C",
  "Neutrophil" = "#6A3D9A",
  "mDC" = "#FF7F00",
  "Monocyte" = "#FFFF99",
  "NK cell" = "#A6761D",
  "B cell" = "#FB9A99",
  "Mast cell" = "#999999"
)

create_compact_pie_chart <- function(data, value_col, title_text) {
  ggplot(data, aes_string(x = "''", y = value_col, fill = "cell_type_simplified")) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.8) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        face = "bold",
        family = "Times",
        margin = margin(b = 5)
      ),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    ggtitle(title_text)
}

responder_plot <- create_compact_pie_chart(responder_data, "Responders", "Responders")
non_responder_plot <- create_compact_pie_chart(non_responder_data, "Non_Responders", "Non-Responders")

legend_data <- data.frame(
  cell_type = names(colors),
  value = 1
)

legend_order <- c("T cell CD8+", "T cell CD4+", "B cell", "NK cell",
                 "Macrophage", "Monocyte", "mDC", "Mast cell", "Neutrophil")
legend_data$cell_type <- factor(legend_data$cell_type, levels = legend_order)

legend_plot <- ggplot(legend_data, aes(x = 1, y = value, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors, name = "") +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10, family = "Times", color = "black"),
    legend.key.size = unit(0.5, "cm"),
    legend.key = element_rect(color = "white", linewidth = 0.5),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 5, 0)
  ) +
  guides(fill = guide_legend(
    nrow = 2,
    byrow = TRUE,
    override.aes = list(color = "white", linewidth = 0.5)
  ))

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- g_legend(legend_plot)

title_grob <- textGrob(
  "Bulk RNA Immune Deconvolution",
  gp = gpar(
    fontsize = 16,
    fontface = "bold",
    fontfamily = "Times"
  )
)

compact_plot <- arrangeGrob(
  title_grob,
  legend,
  arrangeGrob(responder_plot, non_responder_plot, nrow = 1, widths = c(1, 1)),
  nrow = 3,
  heights = c(0.6, 1.0, 3.5)
)

grid.draw(compact_plot)
# Figure 2I: Mirrored Bar Plot - Immune Cell Proportions
library(ggplot2)
library(dplyr)
library(tidyr)

data <- read.csv("data/estimation_matrix_samples.csv")

responders <- c("solid_19", "solid_13", "solid_26", "solid_27")
all_samples <- colnames(data)[-1]
non_responders <- setdiff(all_samples, responders)

cibersort_rows <- grepl("_CIBERSORT$", data$cell_type)
cibersort_data <- data[cibersort_rows, ]

cibersort_data$cell_type_clean <- gsub("_CIBERSORT$", "", cibersort_data$cell_type)

cibersort_data$cell_type_simplified <- cibersort_data$cell_type_clean
cibersort_data$cell_type_simplified <- gsub("T cell CD4\\+ memory resting", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell CD4\\+ memory activated", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell CD4\\+ naive", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell follicular helper", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell regulatory \\(Tregs\\)", "T cell CD4+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("T cell gamma delta", "T cell CD8+", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("B cell naive", "B cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("B cell memory", "B cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("B cell plasma", "B cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("NK cell resting", "NK cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("NK cell activated", "NK cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Macrophage M0", "Macrophage", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Macrophage M1", "Macrophage", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Macrophage M2", "Macrophage", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Myeloid dendritic cell resting", "mDC", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Myeloid dendritic cell activated", "mDC", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Mast cell activated", "Mast cell", cibersort_data$cell_type_simplified)
cibersort_data$cell_type_simplified <- gsub("Mast cell resting", "Mast cell", cibersort_data$cell_type_simplified)

sample_columns <- colnames(cibersort_data)[2:(ncol(cibersort_data)-2)]

aggregated_data <- cibersort_data %>%
  group_by(cell_type_simplified) %>%
  summarise(across(all_of(sample_columns), \(x) sum(x, na.rm = TRUE)), .groups = 'drop')

long_data <- aggregated_data %>%
  pivot_longer(cols = -cell_type_simplified, names_to = "Sample", values_to = "Proportion") %>%
  mutate(Response = ifelse(Sample %in% responders, "Responder", "Non-Responder"))

long_data <- long_data %>%
  group_by(cell_type_simplified) %>%
  filter(sum(Proportion) > 0) %>%
  ungroup()

mirror_data <- long_data %>%
  group_by(cell_type_simplified, Response) %>%
  summarise(Mean = mean(Proportion), .groups = 'drop')

mirror_data <- mirror_data %>%
  mutate(Mean_Mirror = ifelse(Response == "Non-Responder", -Mean, Mean))

cell_order <- mirror_data %>%
  group_by(cell_type_simplified) %>%
  summarise(Total = sum(abs(Mean_Mirror))) %>%
  arrange(desc(Total)) %>%
  pull(cell_type_simplified)

mirror_data$cell_type_simplified <- factor(mirror_data$cell_type_simplified, levels = cell_order)

mirror_plot <- ggplot(mirror_data, aes(x = cell_type_simplified, y = Mean_Mirror, fill = Response)) +
  geom_bar(stat = "identity", position = "identity", width = 0.7) +
  scale_fill_manual(values = c("Responder" = "#4A90E2", "Non-Responder" = "#8B1A1A")) +
  scale_y_continuous(labels = function(x) abs(x),
                     breaks = seq(-0.4, 0.4, 0.1),
                     limits = c(-max(abs(mirror_data$Mean_Mirror)) * 1.1,
                               max(abs(mirror_data$Mean_Mirror)) * 1.1)) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11, family = "Times"),
    axis.title = element_text(size = 12, face = "bold", family = "Times"),
    plot.title = element_text(size = 16, face = "bold", family = "Times", hjust = 0.5),
    plot.subtitle = element_text(size = 12, family = "Times", hjust = 0.5),
    legend.title = element_text(size = 11, face = "bold", family = "Times"),
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Immune Cell Proportions in Responders vs Non-Responders",
    x = "Cell Type",
    y = "Mean Proportion",
    fill = "Response Status"
  ) +
  annotate("text", x = 0.5, y = max(abs(mirror_data$Mean_Mirror)) * 0.9,
           label = "Responders", hjust = 1, size = 4, fontface = "italic", family = "Times") +
  annotate("text", x = 0.5, y = -max(abs(mirror_data$Mean_Mirror)) * 0.9,
           label = "Non-Responders", hjust = 0, size = 4, fontface = "italic", family = "Times")

print(mirror_plot)
# Figure 2J: Spatial Transcriptomics Visualization
library(Seurat)
library(ggplot2)

# Seurat object available in GSEXXXXX
solid <- readRDS('data/solid_all10samples_b2c_malignant_only.rds')

celltype_colors <- c(
  "Malignant" = "#E41A1C",
  "Macrophage" = "#377EB8",
  "Endothel" = "#4DAF4A",
  "Tcell" = "#984EA3",
  "OPC" = "#FF7F00",
  "Pericyte" = "#FFFF33",
  "Astrocyte" = "#A65628",
  "Oligodendrocyte" = "#F781BF",
  "Excitatory neuron" = "#999999",
  "Bcell" = "#66C2A5",
  "Inhibitory neuron" = "#FC8D62"
)

obj <- solid[["solid_19"]]
p <- SpatialDimPlot(obj,
                    group.by = "RCTD_dominant_celltype",
                    pt.size.factor = 10,
                    alpha = 1,
                    image.alpha = 0,
                    cols = celltype_colors) +
     theme(legend.position = "none")

print(p)

obj1 <- solid[["solid_14"]]
p1 <- SpatialDimPlot(obj1,
                     group.by = "RCTD_dominant_celltype",
                     pt.size.factor = 10,
                     alpha = 1,
                     image.alpha = 0,
                     cols = celltype_colors) +
      theme(legend.position = "none")

print(p1)
# Figure 2K: Macrophage and Tumor Proportion Analysis
library(Seurat)
library(dplyr)
library(ggplot2)

# Seurat object available in GSEXXXXX
solid <- readRDS('data/solid_all10samples_b2c_with_rctd_filtered.rds')

responders <- c("solid_13", "solid_19", "solid_26", "solid_27")
non_responders <- c("solid_14", "solid_24", "solid_25", "solid_33", "solid_39", "solid_40")

proportion_data <- list()

for (sample in c(responders, non_responders)) {
  sample_data <- solid[[sample]]
  rctd_cells <- sample_data$RCTD_dominant_celltype
  cell_proportions <- table(rctd_cells) / length(rctd_cells)

  cell_proportions_df <- as.data.frame(cell_proportions)
  colnames(cell_proportions_df) <- c("Cell_Type", "Proportion")
  cell_proportions_df$Sample <- sample
  cell_proportions_df$Group <- ifelse(sample %in% responders, "Responder", "Non-responder")

  proportion_data[[sample]] <- cell_proportions_df
}

final_proportions_df <- do.call(rbind, proportion_data)

final_proportions_df <- final_proportions_df %>%
  filter(Cell_Type %in% c("Macrophage", "Malignant"))

non_responder_means <- final_proportions_df %>%
  filter(Group == "Non-responder") %>%
  group_by(Cell_Type) %>%
  summarize(non_responder_mean = mean(Proportion))

normalized_df <- final_proportions_df %>%
  left_join(non_responder_means, by = "Cell_Type") %>%
  mutate(Normalized_Proportion = Proportion / non_responder_mean)

normalized_df <- normalized_df %>%
  mutate(Cell_Type = recode(Cell_Type,
                            "Malignant" = "Tumor",
                            "Macrophage" = "Macrophage"))

custom_boxplot_stats <- function(y) {
  m <- mean(y)
  s <- sd(y)
  data.frame(
    ymin = m - 1.5 * s,
    lower = m - s,
    middle = m,
    upper = m + s,
    ymax = m + 1.5 * s
  )
}

p <- ggplot(normalized_df, aes(x = Cell_Type, y = Normalized_Proportion, fill = Group)) +
  stat_summary(geom = "boxplot", fun.data = custom_boxplot_stats,
               position = position_dodge(0.8), width = 0.6, alpha = 0.7) +
  geom_hline(yintercept = 1.2, linetype = "dotted", color = "black") +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "black") +
  scale_fill_manual(values = c("Responder" = "#4A90E2", "Non-responder" = "#8B1A1A")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "Normalized Proportion (relative to non-responders)",
       title = "Macrophage & Tumor: Responders vs Non-Responders") +
  theme(axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

print(p)
# Figure 2L: M1/M2 Macrophage Ratio Analysis
library(Seurat)
library(ggplot2)
library(dplyr)

# Seurat object available in GSEXXXXX
all_samples <- readRDS("data/solid_all10samples_b2c_with_rctd_filtered_m1m2.rds")

responders <- c("solid_13", "solid_19", "solid_26", "solid_27")
non_responders <- c("solid_14", "solid_24", "solid_25", "solid_33", "solid_39", "solid_40")

results <- data.frame()

for (sample_name in names(all_samples)) {
  sample_obj <- all_samples[[sample_name]]

  macro_data <- sample_obj@meta.data %>%
    filter(RCTD_dominant_celltype == "Macrophage")

  counts <- table(macro_data$macro_subtype)

  total <- sum(counts)
  m1_count <- ifelse("M1-like" %in% names(counts), counts["M1-like"], 0)
  m2_count <- ifelse("M2-like" %in% names(counts), counts["M2-like"], 0)

  response <- ifelse(sample_name %in% responders, "Responder", "Non-responder")

  results <- rbind(results, data.frame(
    sample = sample_name,
    response = response,
    M1_pct = as.numeric(m1_count) / total * 100,
    M2_pct = as.numeric(m2_count) / total * 100
  ))
}

results$M1M2_ratio <- results$M1_pct / (results$M2_pct + 0.01)

responder_ratios <- results$M1M2_ratio[results$response == "Responder"]
non_responder_ratios <- results$M1M2_ratio[results$response == "Non-responder"]
test_result <- wilcox.test(responder_ratios, non_responder_ratios)
p_value <- test_result$p.value

p_text <- ifelse(p_value < 0.001, "p < 0.001",
                ifelse(p_value < 0.01, sprintf("p = %.3f", p_value),
                      sprintf("p = %.2f", p_value)))

y_max <- max(results$M1M2_ratio) * 1.1

p <- ggplot(results, aes(x = response, y = M1M2_ratio, fill = response)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  scale_fill_manual(values = c("Responder" = "#4A90E2", "Non-responder" = "#8B1A1A")) +
  annotate("text", x = 1.5, y = y_max, label = p_text, size = 5, fontface = "bold") +
  labs(title = "M1/M2 Ratio: Responder vs Non-responder",
       x = "Response Status", y = "M1/M2 Ratio", fill = "Status") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

print(p)
