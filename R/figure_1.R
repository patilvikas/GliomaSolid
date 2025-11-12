# Figure 1: Clinical and Molecular Characterization

# Figure 1A
# No code - no new code was used in the generation of this figure

# Figure 1B
# No code - no new code was used in the generation of this figure

# Figure 1C: Overall Survival Kaplan-Meier
library(survival)
library(survminer)
library(ggsci)

d1 <- read.csv('data/solid_surv_withOS.csv')

if (!all(c("os", "OS_status") %in% names(d1))) {
  stop("Dataframe 'd1' must contain columns named 'os' and 'OS_status'.")
}

surv_obj <- Surv(time = d1$os, event = d1$OS_status)
km_fit <- survfit(surv_obj ~ 1, data = d1)

km_fit_summary_table <- summary(km_fit)$table
median_val <- km_fit_summary_table["median"]
ci_lower <- km_fit_summary_table["0.95LCL"]
ci_upper <- km_fit_summary_table["0.95UCL"]

if (is.na(median_val)) {
  median_os_text <- "Median OS: Not Reached"
} else {
  median_os_text <- paste0("Median OS: ", round(median_val, 1), " months")
  ci_parts <- c()
  if (!is.na(ci_lower)) {
    ci_parts <- c(ci_parts, round(ci_lower, 1))
  } else {
    ci_parts <- c(ci_parts, "NR")
  }
  if (!is.na(ci_upper)) {
    ci_parts <- c(ci_parts, round(ci_upper, 1))
  } else {
    ci_parts <- c(ci_parts, "NR")
  }
  median_os_text <- paste0(median_os_text, " (95% CI: ", ci_parts[1], " - ", ci_parts[2], ")")
}

ggsurv_final_plot <- ggsurvplot(
  km_fit,
  data = d1,
  xlab = "Time (Months)",
  ylab = "Overall Survival Proportion",
  title = "Overall Survival",
  subtitle = median_os_text,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.20,
  palette = pal_npg("nrc")(1),
  size = 1,
  censor = TRUE,
  censor.shape = "|",
  censor.size = 4,
  surv.median.line = 'v',
  risk.table = TRUE,
  risk.table.title = "Number at Risk",
  risk.table.col = "black",
  risk.table.y.text = FALSE,
  risk.table.height = 0.25,
  pval = FALSE,
  legend = "none",
  ggtheme = theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1.5), face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = rel(1.1)),
      axis.title = element_text(size = rel(1.2)),
      axis.text = element_text(size = rel(1.0)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      panel.border = element_blank()
    )
)

ggsurv_final_plot$table <- ggsurv_final_plot$table +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.0, size = rel(1.2))
  )

print(ggsurv_final_plot)

# Figure 1D: Waterfall Plot
library(ggplot2)
library(dplyr)

df <- read.csv('data/MRI_volumetric_data.csv')

df_clean <- df %>% filter(!is.na(percent_original))
df_clean$percent_original <- as.numeric(df_clean$percent_original)
df_clean$percent_change_from_100 <- df_clean$percent_original - 100
df_clean <- df_clean %>% filter(Response != "SD")
df_clean$Response <- factor(df_clean$Response, levels = c("CR", "PR", "PD"))
df_clean <- df_clean %>% arrange(desc(percent_change_from_100))

custom_colors <- c("CR" = "#1F78B4", "PR" = "#33A02C", "PD" = "#E31A1C")
custom_breaks <- c(-200, -100, -50, 0, 50, 100, 200, 500, 1000)

p <- ggplot(df_clean, aes(
  y = reorder(Patient, percent_change_from_100),
  x = percent_change_from_100,
  fill = Response
)) +
  geom_bar(stat = "identity", color = "black", size = 1, width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -50, linetype = "dashed", color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 50),
    breaks = custom_breaks,
    labels = custom_breaks
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  ) +
  labs(y = "Patients", x = "Percent change from baseline (BOR)", fill = "Response", title = "")

print(p)

# Figure 1E: Swimmer Plot
library(ggplot2)

collapsed_data <- read.csv('data/pfs_data.csv')
collapsed_data$SampleNumber <- as.numeric(sub("Solid ", "", collapsed_data$Sample))

collapsed_data$Category <- "Progressive Disease"
collapsed_data$Category[collapsed_data$Sample %in% c("Solid 26", "Solid 27", "Solid 13")] <- "Partial Responder"
collapsed_data$Category[collapsed_data$Sample == "Solid 19"] <- "Complete Responder"

collapsed_data$Subject <- sprintf("SOLID-%03d", as.numeric(sub("Solid ", "", collapsed_data$Sample)))
alive_subjects <- c("SOLID-009", "SOLID-019", "SOLID-024", "SOLID-027", "SOLID-033", "SOLID-034")

collapsed_data$marker_type <- ifelse(collapsed_data$Subject %in% alive_subjects, "Alive", "Deceased")
collapsed_data$marker_type <- factor(collapsed_data$marker_type, levels = c("Deceased", "Alive"))

ggplot(collapsed_data, aes(x = reorder(Sample, pfs), y = pfs)) +
  geom_bar(stat = "identity", aes(fill = Category)) +
  geom_point(data = subset(collapsed_data, marker_type == "Deceased"),
             aes(y = pfs, shape = marker_type), color = "black", size = 3) +
  geom_point(data = subset(collapsed_data, marker_type == "Alive"),
             aes(y = pfs, shape = marker_type), color = "black", size = 4) +
  coord_flip() +
  scale_fill_manual(values = c("#1F78B4", "#33A02C", "#E31A1C")) +
  scale_shape_manual(name = "Status", values = c("Deceased" = 16, "Alive" = 17)) +
  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(x = "Sample", y = "PFS (Months)",
       title = "Swimmer Plot: Progression-Free Survival and Probability of Response",
       fill = "Response") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = c(0.85, 0.7)) +
  ylim(0, 35)

# Figure 1F
# No code - no new code was used in the generation of this figure

# Figure 1G: Oncoplot
library(maftools)
library(dplyr)

d1 <- read.maf('data/exonic_variants.maf')

sample_numbers <- as.character(1:28)
solid_labels <- c("solid_24", "solid_32", "solid_40", "solid_37", "solid_39",
                  "solid_35", "solid_28", "solid_30", "solid_29", "solid_27",
                  "solid_25", "solid_26", "solid_38", "solid_33", "solid_34",
                  "solid_36", "solid_31", "solid_17", "solid_2", "solid_22",
                  "solid_5", "solid_6", "solid_23", "solid_8", "solid_11",
                  "solid_14", "solid_19", "solid_4")

sample_mapping <- data.frame(MAF_Label = sample_numbers, Solid_Label = solid_labels, stringsAsFactors = FALSE)

d1@data <- d1@data %>%
  left_join(sample_mapping, by = c("Tumor_Sample_Barcode" = "MAF_Label")) %>%
  mutate(Tumor_Sample_Barcode = Solid_Label) %>%
  select(-Solid_Label)

responders <- c("solid_19", "solid_13", "solid_26", "solid_27")
d1@data$Response_Status <- ifelse(d1@data$Tumor_Sample_Barcode %in% responders, "Responder", "Non-responder")

tmb_results <- tmb(maf = d1)
sample_mapping <- data.frame(Tumor_Sample_Barcode = sample_numbers, Solid_Label = solid_labels, stringsAsFactors = FALSE)

tmb_results_mapped <- tmb_results %>%
  left_join(sample_mapping, by = "Tumor_Sample_Barcode") %>%
  select(Solid_Label, total, total_perMB, total_perMB_log) %>%
  rename(Tumor_Sample_Barcode = Solid_Label)

tmb_results_mapped$Response_Status <- ifelse(tmb_results_mapped$Tumor_Sample_Barcode %in% responders,
                                             "Responder", "Non-responder")

clinical_annotation <- data.frame(
  Tumor_Sample_Barcode = tmb_results_mapped$Tumor_Sample_Barcode,
  Response_Status = ifelse(tmb_results_mapped$Tumor_Sample_Barcode %in% responders, "Responder", "Non-responder"),
  stringsAsFactors = FALSE
)

genes_updated <- c("TP53", "ATRX", "CIC", "NF1", "PIK3CA", "PIK3R1", "PTEN", "SMARCA4",
                   "EGFR", "IDH2", "RB1", "RBM47", "PDGFRA", "GRIN2A", "BCOR", "RELN",
                   "TRRAP", "HUWE1", "BPTF", "NOTCH1", "NOTCH2", "MECOM")

oncoplot(maf = d1,
         clinicalFeatures = "Response_Status",
         sortByAnnotation = T,
         annotationColor = list(Response_Status = c("Responder" = "darkgreen", "Non-responder" = "darkred")),
         annotationDat = clinical_annotation,
         genes = genes_updated,
         draw_titv = FALSE,
         removeNonMutated = F,
         drawRowBar = FALSE,
         drawColBar = FALSE,
         topPathways = 1,
         showTumorSampleBarcodes = F,
         keepGeneOrder = F,
         sortByMutation = F,
         showPct = FALSE,
         showTitle = FALSE)

# Figure 1H: TMB Comparison
library(maftools)
library(dplyr)
library(ggplot2)
library(ggpubr)

d1 <- read.maf('data/exonic_variants.maf')

sample_numbers <- as.character(1:28)
solid_labels <- c("solid_24", "solid_32", "solid_40", "solid_37", "solid_39",
                  "solid_35", "solid_28", "solid_30", "solid_29", "solid_27",
                  "solid_25", "solid_26", "solid_38", "solid_33", "solid_34",
                  "solid_36", "solid_31", "solid_17", "solid_2", "solid_22",
                  "solid_5", "solid_6", "solid_23", "solid_8", "solid_11",
                  "solid_14", "solid_19", "solid_4")

sample_mapping <- data.frame(MAF_Label = sample_numbers, Solid_Label = solid_labels, stringsAsFactors = FALSE)

d1@data <- d1@data %>%
  left_join(sample_mapping, by = c("Tumor_Sample_Barcode" = "MAF_Label")) %>%
  mutate(Tumor_Sample_Barcode = Solid_Label) %>%
  select(-Solid_Label)

tmb_results <- tmb(maf = d1)

sample_mapping <- data.frame(Tumor_Sample_Barcode = sample_numbers, Solid_Label = solid_labels, stringsAsFactors = FALSE)

tmb_results_mapped <- tmb_results %>%
  left_join(sample_mapping, by = "Tumor_Sample_Barcode") %>%
  select(Solid_Label, total, total_perMB, total_perMB_log) %>%
  rename(Tumor_Sample_Barcode = Solid_Label)

responders <- c("solid_19", "solid_13", "solid_26", "solid_27")
tmb_results_mapped$Response_Status <- ifelse(tmb_results_mapped$Tumor_Sample_Barcode %in% responders,
                                             "Responder", "Non-responder")

ggplot(tmb_results_mapped, aes(x = Response_Status, y = total_perMB, fill = Response_Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.15, size = 3, alpha = 1, color = "black") +
  labs(title = "TMB Comparison: Responders vs Non-responders",
       x = "Response Status", y = "Tumor Mutation Burden (Mut/Mb)") +
  theme_minimal() +
  scale_fill_manual(values = c("Responder" = "darkgreen", "Non-responder" = "darkred")) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, label.y = max(tmb_results_mapped$total_perMB) + 0.1) +
  theme(legend.position = "right")
