# Extended Data Figure 4

# Extended Data Figure 4A: LIMMA Differential Expression Analysis
# RNA counts data available in GSEXXXXX (download to data/ folder)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

count_data <- read.table("data/rna_counts.txt", header = TRUE, sep = "\t",
                         check.names = FALSE)

gene_names <- count_data$gene_symbol
gene_names <- make.unique(gene_names, sep = "_")
rownames(count_data) <- gene_names
count_data <- count_data[, -1]

responders <- c("solid_19", "solid_13", "solid_26", "solid_27")
all_samples <- colnames(count_data)
non_responders <- setdiff(all_samples, responders)

condition <- factor(ifelse(all_samples %in% responders, "Responder", "Non_responder"),
                   levels = c("Non_responder", "Responder"))
design <- model.matrix(~condition)
colnames(design) <- c("Intercept", "ResponderVsNonResponder")

count_matrix <- as.matrix(count_data[, all_samples])

keep <- rowSums(count_matrix) >= 10
count_matrix_filtered <- count_matrix[keep, ]

dge <- DGEList(counts = count_matrix_filtered)
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = "ResponderVsNonResponder",
                   number = Inf, sort.by = "P", p.value = 1)

results$gene <- rownames(results)
results$neg_log10_pvalue <- -log10(results$P.Value)
results$neg_log10_adj_pvalue <- -log10(results$adj.P.Val)

padj_threshold <- 0.1
pvalue_threshold <- 0.05
lfc_threshold <- 1.0

results$color_category <- "Not Significant"
results$color_category[results$adj.P.Val < padj_threshold & results$logFC > lfc_threshold] <- "Upregulated (padj < 0.1)"
results$color_category[results$adj.P.Val < padj_threshold & results$logFC < -lfc_threshold] <- "Downregulated (padj < 0.1)"
results$color_category[results$P.Value < pvalue_threshold & results$adj.P.Val >= padj_threshold & results$logFC > lfc_threshold] <- "Upregulated (p < 0.05)"
results$color_category[results$P.Value < pvalue_threshold & results$adj.P.Val >= padj_threshold & results$logFC < -lfc_threshold] <- "Downregulated (p < 0.05)"

immune_genes <- c(
  "IFNG", "IFNGR1", "IFNGR2", "STAT1", "STAT2", "IRF1", "IRF2", "IRF7", "IRF8", "IRF9",
  "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL2", "CCL3", "CCL4", "CCL5", "CCL19", "CCL21",
  "IL2", "IL6", "IL12A", "IL12B", "IL15", "IL18", "IL21", "TNF", "TNFA", "TNFRSF9",
  "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CD28", "CTLA4", "PDCD1", "LAG3",
  "TIGIT", "HAVCR2", "BTLA", "CD274", "PDCD1LG2", "ICOS", "ICOSLG",
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L",
  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6",
  "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
  "HLA-DPA1", "HLA-DPB1", "HLA-DPB2",
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  "B2M", "TAP1", "TAP2", "TAPBP", "CALR", "CANX", "PDIA3",
  "CD47", "SIRPA", "CD200", "CD200R1", "VSIR", "NECTIN2", "PVRIG",
  "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1",
  "KLRK1", "KLRD1", "KLRC1", "KLRC2", "KLRC3", "KLRC4", "NCAM1", "FCGR3A",
  "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CFB", "CFD",
  "NLRP3", "TLR4", "TLR7", "TLR8", "TLR9", "MYD88", "TRAF6",
  "FOXP3", "IL2RA", "IL10", "TGFB1", "IDO1", "ARG1"
)

parp_sensitive_genes <- c(
  "BRCA1", "BRCA2", "PALB2", "BARD1", "BRIP1", "RAP80", "ABRAXAS1",
  "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAD52", "RAD54L", "RAD54B",
  "XRCC2", "XRCC3", "RAD55", "RAD57", "RAD59",
  "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG",
  "FANCI", "FANCJ", "FANCL", "FANCM", "FANCN", "FANCO", "FANCP", "FANCQ",
  "FANCR", "FANCS", "FANCT", "FANCU", "FANCV", "FANCW", "SLX4", "ERCC4",
  "ATM", "ATR", "ATRX", "CHEK1", "CHEK2", "WEE1", "CDC25A", "CDC25C",
  "TP53", "TP53BP1", "MDC1", "H2AX", "RIF1", "53BP1", "TOPBP1",
  "CDK1", "CDK2", "CDK4", "CDK6", "CDK12", "CCNE1", "CCND1", "CCNA2",
  "PARP1", "PARP2", "PARP3", "PARP4", "PARP5A", "PARP5B", "PARP6", "PARP7",
  "PARP8", "PARP9", "PARP10", "PARP11", "PARP12", "PARP13", "PARP14", "PARP15", "PARP16",
  "XRCC1", "POLB", "POLD1", "POLE", "LIG1", "LIG3", "APEX1", "NEIL1", "NEIL2", "NEIL3",
  "OGG1", "MUTYH", "UNG", "SMUG1", "TDG", "MBD4",
  "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2", "PCNA", "RFC1",
  "XPA", "XPC", "XPD", "XPF", "XPG", "ERCC1", "ERCC2", "ERCC3", "ERCC5", "ERCC6", "ERCC8",
  "RPA1", "RPA2", "RPA3", "REPA1", "PCNA", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5",
  "POLD1", "POLD2", "POLD3", "POLD4", "POLE", "POLE2", "POLE3", "POLE4",
  "NBN", "MRE11A", "RAD50", "RBBP8", "EXO1", "DNA2", "BLM", "WRN", "RECQL4", "RECQL5",
  "SLFN11", "SCHIP1", "CCDC6", "EZH2", "MYC", "RB1", "PTEN", "PIK3CA", "AKT1"
)

results$gene_category <- "Other"
results$gene_category[results$gene %in% parp_sensitive_genes] <- "PARP Sensitive"
results$gene_category[results$gene %in% immune_genes] <- "Immune"

volcano_plot_limma <- ggplot(results, aes(x = logFC, y = neg_log10_pvalue)) +
  geom_point(data = subset(results, gene_category == "Other"),
             color = "gray70", alpha = 0.4, size = 0.8) +
  geom_point(data = subset(results, gene_category == "PARP Sensitive"),
             color = "blue", alpha = 0.8, size = 1.5) +
  geom_point(data = subset(results, gene_category == "Immune"),
             color = "red", alpha = 0.8, size = 1.5) +
  geom_point(data = subset(results, gene_category == "Other" & adj.P.Val < padj_threshold),
             color = "darkgray", alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray30", alpha = 0.3, linewidth = 0.5) +
  annotate("text", x = 3.5, y = max(results$neg_log10_pvalue, na.rm = TRUE) * 0.95,
           label = "Immunotherapy responsive",
           size = 6, hjust = 1, vjust = 1, fontface = "bold",
           color = "red", family = "Times") +
  annotate("text", x = 3.5, y = max(results$neg_log10_pvalue, na.rm = TRUE) * 0.88,
           label = "PARP responsive",
           size = 6, hjust = 1, vjust = 1, fontface = "bold",
           color = "blue", family = "Times") +
  labs(x = expression(log[2]*"(Fold Change)"),
       y = expression(-log[10]*"(P-value)"),
       title = "Differential Gene Expression of Responders") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1)) +
  scale_y_continuous(limits = c(0, max(results$neg_log10_pvalue, na.rm = TRUE) * 1.05)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11, family = "Times"),
    axis.title = element_text(color = "black", size = 13, family = "Times"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Times"),
    legend.position = "none",
    text = element_text(family = "Times")
  )

print(volcano_plot_limma)

# Extended Data Figure 4B: Immunotherapy vs PARP Pathways Comparison
# Note: May require processed GSEA data in data folder
library(ggplot2)
library(dplyr)

gsea_results <- read.csv("data/gsea_results.csv", stringsAsFactors = FALSE)

positive_pathways <- gsea_results[gsea_results$NES > 0, ]

positive_pathways$pathway_clean <- gsub("HALLMARK_", "", positive_pathways$pathway)
positive_pathways$pathway_clean <- gsub("GOBP_", "", positive_pathways$pathway_clean)
positive_pathways$pathway_clean <- gsub("GOMF_", "", positive_pathways$pathway_clean)
positive_pathways$pathway_clean <- gsub("KEGG_MEDICUS_REFERENCE_", "", positive_pathways$pathway_clean)
positive_pathways$pathway_clean <- gsub("REACTOME_", "", positive_pathways$pathway_clean)
positive_pathways$pathway_clean <- gsub("_", " ", positive_pathways$pathway_clean)

immunotherapy_keywords <- c("INTERFERON", "IMMUNE", "INFLAMMATION", "CYTOKINE", "CHEMOKINE",
                           "DENDRITIC", "T CELL", "B CELL", "NK", "NATURAL KILLER",
                           "ANTIGEN", "MHC", "HLA", "TUMOR", "ALLOGRAFT", "REJECTION",
                           "LYMPHOCYTE", "MACROPHAGE", "MONOCYTE", "NEUTROPHIL")

parp_keywords <- c("DNA REPAIR", "DNA DAMAGE", "CELL CYCLE CHECKPOINT", "APOPTOSIS", "DNA REPLICATION",
                  "CHECKPOINT", "HOMOLOGOUS RECOMBINATION", "BASE EXCISION",
                  "NUCLEOTIDE EXCISION", "MISMATCH REPAIR", "FANCONI", "BRCA", "PARP",
                  "DOUBLE STRAND BREAK", "DSB", "ATM", "ATR", "P53", "TP53",
                  "MITOTIC", "S PHASE", "G1 S", "G2 M", "CENTROSOME", "SPINDLE")

positive_pathways$pathway_category <- "Other"

for(keyword in immunotherapy_keywords) {
  positive_pathways$pathway_category[grepl(keyword, toupper(positive_pathways$pathway_clean))] <- "Immunotherapy"
}

exclude_keywords <- c("EXTRACELLULAR MATRIX", "MATRIX STRUCTURAL", "COLLAGEN", "BASEMENT MEMBRANE",
                     "TUMOR", "IMMUNE", "INFLAMMATION", "DENDRITIC", "CHEMOTAXIS")
for(keyword in parp_keywords) {
  is_parp <- grepl(keyword, toupper(positive_pathways$pathway_clean)) & positive_pathways$pathway_category == "Other"
  is_excluded <- sapply(positive_pathways$pathway_clean, function(x) any(sapply(exclude_keywords, function(excl) grepl(excl, toupper(x)))))
  positive_pathways$pathway_category[is_parp & !is_excluded] <- "PARP"
}

immunotherapy_pathways <- positive_pathways[positive_pathways$pathway_category == "Immunotherapy", ]
top_immune <- immunotherapy_pathways[immunotherapy_pathways$pval < 0.05, ]
top_immune <- top_immune[order(-top_immune$NES), ][1:15, ]

parp_pathways <- positive_pathways[positive_pathways$pathway_category == "PARP", ]
parp_pathways <- parp_pathways[order(-parp_pathways$NES), ]
top5_parp <- head(parp_pathways, 5)

plot_data <- rbind(top_immune, top5_parp)
plot_data <- plot_data[order(-plot_data$NES), ]

min_pval <- min(plot_data$pval)
max_pval <- max(plot_data$pval)
min_nes <- min(plot_data$NES)
max_nes <- max(plot_data$NES)

label_colors <- rep("black", nrow(plot_data))

comparison_plot <- ggplot(plot_data, aes(x = NES, y = reorder(pathway_clean, NES))) +
  geom_point(aes(size = abs(NES), color = pval), alpha = 0.8) +
  scale_color_gradient(low = "red", high = "pink",
                      name = "p-value",
                      trans = "log10",
                      breaks = c(0.001, 0.01, 0.1),
                      labels = c("0.001", "0.01", "0.1"),
                      limits = c(min_pval, max_pval)) +
  scale_size_continuous(name = "NES", range = c(3, 8),
                       breaks = c(1.0, 1.2, 1.4, 1.6, 1.8),
                       labels = c("1.0", "1.2", "1.4", "1.6", "1.8"),
                       limits = c(1.0, 1.8),
                       guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  labs(x = "Normalized Enrichment Score (NES)", y = "",
       title = "Gene Pathway Enrichment") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = label_colors, size = 11, hjust = 1, family = "Times"),
    axis.text.x = element_text(color = "black", size = 11, family = "Times"),
    axis.title.x = element_text(color = "black", size = 13, face = "bold", family = "Times"),
    axis.title.y = element_text(color = "black", size = 13, face = "bold", family = "Times"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "Times", color = "black"),
    text = element_text(family = "Times", color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(color = "black", face = "bold"),
    legend.text = element_text(color = "black"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_x_continuous(limits = c(min_nes - 0.05, max_nes + 0.1))

print(comparison_plot)

# Extended Data Figure 4C
# JL to add script and processed data

# Extended Data Figure 4D: Immune Signatures Comparison
# RNA counts data available in GSEXXXXX (download to data/ folder)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

count_data <- read.table("data/rna_counts.txt", header = TRUE, sep = "\t",
                         check.names = FALSE)

gene_names <- count_data$gene_symbol
gene_names <- make.unique(gene_names, sep = "_")
rownames(count_data) <- gene_names
count_matrix <- count_data[, -1]

responders <- c("solid_19", "solid_13", "solid_26", "solid_27")
all_samples <- colnames(count_matrix)
non_responders <- setdiff(all_samples, responders)

tpm_matrix <- apply(count_matrix, 2, function(x) (x / sum(x)) * 1e6)

stat1_genes <- c(
  "TAP1", "GBP1", "IFIH1", "PSMB9", "CXCL9", "IRF1",
  "CXCL11", "CXCL10", "IDO1", "STAT1"
)

tis_genes <- c(
  "TIGIT", "CD27", "CD8A", "PDCD1LG2", "CXCR6", "LAG3",
  "CD274", "CMKLR1", "NKG7", "CCL5", "PSMB10", "IDO1",
  "PPBP", "HLA-DQA1", "CD276", "STAT1", "HLA-DRB1", "HLA-E"
)

parpi7_genes <- c(
  "BRCA1", "CHEK2", "MAPKAPK2", "MRE11A", "NBN", "TDG", "XPA"
)

calculate_zscore_signature <- function(expression_matrix, signature_genes) {
  available_genes <- intersect(signature_genes, rownames(expression_matrix))

  if(length(available_genes) == 0) {
    stop("No signature genes found!")
  }

  gene_data <- log2(expression_matrix[available_genes, , drop = FALSE] + 1)
  gene_data_centered <- t(scale(t(gene_data), center = TRUE, scale = FALSE))
  mean_scores <- colMeans(gene_data_centered, na.rm = TRUE)
  z_scores <- scale(mean_scores)[,1]

  return(z_scores)
}

stat1_scores <- calculate_zscore_signature(tpm_matrix, stat1_genes)
tis_scores <- calculate_zscore_signature(tpm_matrix, tis_genes)
parpi7_scores <- calculate_zscore_signature(tpm_matrix, parpi7_genes)

results_df <- data.frame(
  Sample = colnames(tpm_matrix),
  Group = factor(ifelse(colnames(tpm_matrix) %in% responders, "Responder", "Non-Responder"),
                levels = c("Non-Responder", "Responder")),
  STAT1 = stat1_scores,
  TIS = tis_scores,
  PARPi7 = parpi7_scores,
  stringsAsFactors = FALSE
)

stat1_test <- wilcox.test(STAT1 ~ Group, data = results_df)
tis_test <- wilcox.test(TIS ~ Group, data = results_df)
parpi7_test <- wilcox.test(PARPi7 ~ Group, data = results_df)

calculate_cohens_d <- function(x, y) {
  pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))
  (mean(x) - mean(y)) / pooled_sd
}

responder_data <- results_df[results_df$Group == "Responder", ]
nonresponder_data <- results_df[results_df$Group == "Non-Responder", ]

stat1_d <- calculate_cohens_d(responder_data$STAT1, nonresponder_data$STAT1)
tis_d <- calculate_cohens_d(responder_data$TIS, nonresponder_data$TIS)
parpi7_d <- calculate_cohens_d(responder_data$PARPi7, nonresponder_data$PARPi7)

colors <- c("Non-Responder" = "#3498DB", "Responder" = "#E74C3C")

create_clean_box_plot <- function(data, y_var, title_text, p_value, effect_size, y_label) {
  ggplot(data, aes_string(x = "Group", y = y_var, fill = "Group")) +
    geom_boxplot(alpha = 0.8, width = 0.6, color = "black", linewidth = 0.5, outlier.size = 2) +
    scale_fill_manual(values = colors) +
    labs(
      title = title_text,
      x = "",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "Times"),
      axis.text.x = element_text(size = 12, family = "Times", color = "black"),
      axis.text.y = element_text(size = 11, family = "Times", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold", family = "Times"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", linewidth = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    annotate("text", x = 2, y = min(data[[y_var]]) * 0.9,
             label = paste("p =", formatC(p_value, format = "g", digits = 3)),
             hjust = 1, vjust = 0, size = 3.5, color = "black", family = "Times")
}

stat1_plot <- create_clean_box_plot(
  results_df, "STAT1",
  "STAT1 signature (Rody et al)",
  stat1_test$p.value, stat1_d,
  "Z-score"
)

tis_plot <- create_clean_box_plot(
  results_df, "TIS",
  "TIS signature (Ayers et al)",
  tis_test$p.value, tis_d,
  "Z-score"
)

parpi7_plot <- create_clean_box_plot(
  results_df, "PARPi7",
  "PARPi7 signature (Pusztai et al)",
  parpi7_test$p.value, parpi7_d,
  "Z-score"
)

immune_signatures_plot <- grid.arrange(
  stat1_plot, tis_plot, parpi7_plot,
  nrow = 1,
  widths = c(1, 1, 1),
  top = textGrob("Response Signatures in ICI+PARP Responders",
                gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "Times"))
)

grid.draw(immune_signatures_plot)

# Extended Data Figure 4E: RNA Signature Blind Prediction
# Large files available in GSEXXXXX (download to data/ folder):
# - rna_counts.txt (3.9MB)
# - TCGA_clinical_and_log_matrix.RData (156MB)
library(edgeR)
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
library(GSVA)
library(gridExtra)

signature <- read.csv("data/rna_signature_response.csv", stringsAsFactors = FALSE)
signature <- signature[order(abs(signature$logFC), decreasing = TRUE), ]

count_data <- read.table("data/rna_counts.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(count_data) <- make.unique(count_data$gene_symbol)
count_data <- count_data[, -1]

clinical <- read.csv("data/baseline_data.csv", stringsAsFactors = FALSE)
clinical$solid_patient_clean <- tolower(gsub("-", "_", clinical$solid_patient))
rownames(clinical) <- clinical$solid_patient_clean

dge <- DGEList(counts = count_data)
keep <- rowSums(count_data > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
log_matrix <- cpm(dge, log = TRUE)

signature_genes <- head(signature$gene, 30)
signature_genes_found <- intersect(signature_genes, rownames(log_matrix))

geneSets <- list(signature = signature_genes_found)
ssgsea_param <- ssgseaParam(log_matrix, geneSets)
scores_matrix <- gsva(ssgsea_param, verbose=FALSE)
scores <- scores_matrix[1,]

score_df <- data.frame(
  Sample = names(scores),
  Score = as.numeric(scores),
  stringsAsFactors = FALSE
)
score_df$sample_clean <- tolower(gsub("-", "_", score_df$Sample))

merged_data <- merge(score_df, clinical,
                    by.x = "sample_clean",
                    by.y = "solid_patient_clean",
                    all.x = FALSE)

median_cutoff <- median(merged_data$Score)
merged_data$predicted_group <- ifelse(merged_data$Score > median_cutoff, "High Score", "Low Score")

surv_obj <- Surv(time = merged_data$pfs, event = rep(1, nrow(merged_data)))
fit <- survfit(surv_obj ~ predicted_group, data = merged_data)

p_solid <- ggsurvplot(
  fit, data = merged_data, pval = TRUE, pval.coord = c(25, 0.2),
  conf.int = FALSE, risk.table = TRUE, palette = c("darkred", "darkblue"),
  title = "30-Gene RNA Signature: BLIND Prediction (SOLID Trial)",
  xlab = "Time (months)", ylab = "Progression-Free Survival",
  legend.title = "Predicted Group", legend.labs = c("High Score", "Low Score"),
  risk.table.height = 0.25, tables.theme = theme_cleantable(),
  ggtheme = theme_minimal() + theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 11), axis.title = element_text(size = 12)
  )
)

load("data/TCGA_clinical_and_log_matrix.RData")
clinical_tcga <- TCGA_GBM_clinical_unique
log_matrix_tcga <- TCGA_GBM_TPM_gene_log_unique

signature_genes_found_tcga <- intersect(signature_genes, rownames(log_matrix_tcga))
common_samples_tcga <- intersect(rownames(clinical_tcga), colnames(log_matrix_tcga))
clinical_matched_tcga <- clinical_tcga[common_samples_tcga, ] %>%
  filter(!is.na(OS_months) & !is.na(OS_status))

geneSets_tcga <- list(signature = signature_genes_found_tcga)
ssgsea_param_tcga <- ssgseaParam(log_matrix_tcga, geneSets_tcga)
scores_matrix_tcga <- gsva(ssgsea_param_tcga, verbose=FALSE)
scores_tcga <- scores_matrix_tcga[1,]

median_cutoff_tcga <- median(scores_tcga)
groups_tcga <- ifelse(scores_tcga > median_cutoff_tcga, "High Score", "Low Score")

score_names <- names(scores_tcga)
clinical_names <- rownames(clinical_matched_tcga)
matched_names <- intersect(score_names, clinical_names)
clinical_matched_tcga <- clinical_matched_tcga[matched_names, ]
groups_matched <- groups_tcga[matched_names]
clinical_matched_tcga$predicted_group <- groups_matched

surv_obj_tcga <- Surv(time = clinical_matched_tcga$OS_months, event = clinical_matched_tcga$OS_status)
fit_tcga <- survfit(surv_obj_tcga ~ predicted_group, data = clinical_matched_tcga)

p_tcga <- ggsurvplot(
  fit_tcga, data = clinical_matched_tcga, pval = TRUE, pval.coord = c(25, 0.2),
  conf.int = FALSE, risk.table = TRUE, palette = c("darkred", "darkblue"),
  title = "30-Gene RNA Signature: BLIND Prediction (TCGA GBM)",
  xlab = "Time (months)", ylab = "Overall Survival Probability",
  legend.title = "Predicted Group", legend.labs = c("High Score", "Low Score"),
  risk.table.height = 0.25, tables.theme = theme_cleantable(),
  ggtheme = theme_minimal() + theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 11), axis.title = element_text(size = 12)
  )
)

combined_plot <- arrange_ggsurvplots(list(p_solid, p_tcga), ncol = 1, nrow = 2, risk.table.height = 0.25)
print(combined_plot)
