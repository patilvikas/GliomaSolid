# Extended Data Figure 3

# Extended Data Figure 3A: CRAWDAD Colocalization Analysis - Responders
# CRAWDAD results available in GSEXXXXX

# Extended Data Figure 3B: CRAWDAD Colocalization Analysis - Non-Responders
# CRAWDAD results available in GSEXXXXX

# Extended Data Figure 3C: Spatial MP Visualization - Responder and Non-Responder
library(Seurat)
library(ggplot2)

# Seurat object available in GSEXXXXX
solid <- readRDS('data/solid_all10samples_b2c_malignant_only_with_MP_filtered.rds')

mp_colors <- c(
  "Stress1" = "#E41A1C",
  "Cilia" = "#377EB8",
  "OPC" = "#4DAF4A",
  "CC" = "#984EA3",
  "AC" = "#FF7F00",
  "Hypoxia" = "#FFFF33",
  "MES" = "#A65628",
  "NPC" = "#F781BF",
  "GPC" = "#999999",
  "ExN" = "#66C2A5"
)

obj <- solid[["solid_19"]]
p <- SpatialDimPlot(obj,
                    group.by = "MP_dominant",
                    pt.size.factor = 10,
                    alpha = 1,
                    image.alpha = 0,
                    cols = mp_colors) +
     theme(legend.position = "none")

print(p)

obj1 <- solid[["solid_14"]]
p1 <- SpatialDimPlot(obj1,
                     group.by = "MP_dominant",
                     pt.size.factor = 10,
                     alpha = 1,
                     image.alpha = 0,
                     cols = mp_colors) +
      theme(legend.position = "none")

print(p1)

# Extended Data Figure 3D: Entropy Boxplot - Malignant MP Microenvironment
# File available in GSEXXXXX (download to data/ folder):
# - hoodscan_malignant_MP_metrics.csv (38MB)
library(ggplot2)

metrics_df <- read.csv("data/hoodscan_malignant_MP_metrics.csv")

test_result <- wilcox.test(Entropy ~ Label, data = metrics_df)
p_value <- test_result$p.value

p_text <- ifelse(p_value < 0.001, "p < 0.001",
                ifelse(p_value < 0.01, sprintf("p = %.3f", p_value),
                      sprintf("p = %.2f", p_value)))

y_max <- max(metrics_df$Entropy) * 1.05

p7 <- ggplot(metrics_df, aes(x = "", y = Entropy, fill = Label)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.size = 1.5) +
  scale_fill_manual(values = c("non_responder" = "#E57373", "responder" = "#4DB6AC"),
                    labels = c("non_responder" = "Non responder", "responder" = "Responder")) +
  annotate("text", x = 1, y = y_max, label = p_text, size = 5, fontface = "bold") +
  labs(title = "Malignant MP microenvironment",
       x = NULL,
       y = "Entropy",
       fill = "") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "cm"),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank()
  )

print(p7)

# Extended Data Figure 3E: Entropy Boxplot - Non-Malignant Microenvironment
# File available in GSEXXXXX (download to data/ folder):
# - hoodscan_metrics_non_malignant.csv (16MB)
library(ggplot2)

metrics_df <- read.csv("data/hoodscan_metrics_non_malignant.csv")

test_result <- wilcox.test(Entropy ~ Label, data = metrics_df)
p_value <- test_result$p.value

p_text <- ifelse(p_value < 0.001, "p < 0.001",
                ifelse(p_value < 0.01, sprintf("p = %.3f", p_value),
                      sprintf("p = %.2f", p_value)))

y_max <- max(metrics_df$Entropy) * 1.05

p <- ggplot(metrics_df, aes(x = "", y = Entropy, fill = Label)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.size = 1.5) +
  scale_fill_manual(values = c("non_responder" = "#E57373", "responder" = "#4DB6AC"),
                    labels = c("non_responder" = "Non responder", "responder" = "Responder")) +
  annotate("text", x = 1, y = y_max, label = p_text, size = 5, fontface = "bold") +
  labs(title = "Non-malignant microenvironment",
       x = NULL,
       y = "Entropy",
       fill = "") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "cm"),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank()
  )

print(p)
