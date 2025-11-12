# Extended Data Figure 2

# Extended Data Figure 2A
# VP to add processed data and script

# Extended Data Figure 2B
# VP to add processed data and script

# Extended Data Figure 2C: TMB across cohorts
# Large file available in GSEXXXXX (download to data/ folder):
# - glass_snv.tsv (169MB, optional)
library(ggplot2)
library(dplyr)
library(scales)

input_file <- "data/tmb_data.txt"

tmb_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tmb_data$total <- pmax(0.1, round(tmb_data$total / 35.8, 1))

glass_file <- "data/glass_snv.tsv"
if(file.exists(glass_file)) {
  glass_raw <- read.table(glass_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  glass_tmb <- glass_raw %>%
    group_by(case_barcode) %>%
    summarise(variant_count = n(), .groups = 'drop') %>%
    mutate(
      total = pmax(0.1, round(variant_count / 3000, 1)),
      Tumor_Sample_Barcode = case_barcode,
      cohort = "GLASS"
    ) %>%
    select(Tumor_Sample_Barcode, total, cohort)

  tmb_data <- rbind(tmb_data, glass_tmb)
}

cohort_order <- tmb_data %>%
  group_by(cohort) %>%
  summarise(median_tmb = median(total, na.rm = TRUE)) %>%
  arrange(median_tmb) %>%
  pull(cohort)

tmb_data$cohort <- ifelse(tmb_data$cohort == "LGG", "TCGA-LGG",
                   ifelse(tmb_data$cohort == "GBM", "TCGA-GBM", tmb_data$cohort))

cohort_order <- gsub("^LGG$", "TCGA-LGG", cohort_order)
cohort_order <- gsub("^GBM$", "TCGA-GBM", cohort_order)

tmb_data$cohort <- factor(tmb_data$cohort, levels = cohort_order)

tmb_data <- tmb_data %>%
  arrange(cohort, total)

cohort_width <- 300
cohort_info <- tmb_data %>%
  group_by(cohort) %>%
  summarise(
    count = n(),
    median_tmb = median(total),
    .groups = 'drop'
  ) %>%
  mutate(
    cohort_start = (as.numeric(cohort) - 1) * cohort_width,
    cohort_end = as.numeric(cohort) * cohort_width,
    mid_index = (cohort_start + cohort_end) / 2,
    min_index = cohort_start + 1,
    max_index = cohort_end
  )

tmb_data <- tmb_data %>%
  group_by(cohort) %>%
  mutate(
    cohort_rank = row_number(),
    cohort_n = n(),
    index = cohort_info$cohort_start[match(cohort, cohort_info$cohort)] +
            (cohort_rank - 1) * (cohort_width - 1) / pmax(cohort_n - 1, 1),
    point_color = ifelse(cohort == "SOLID", "black", "grey60")
  ) %>%
  ungroup()

cohort_info$band_alpha <- rep(c(0.15, 0), length.out = nrow(cohort_info))

p <- ggplot(tmb_data, aes(x = index, y = total)) +
  geom_rect(data = cohort_info,
            aes(xmin = min_index - 0.5, xmax = max_index + 0.5,
                ymin = 0.05, ymax = 1000, alpha = band_alpha),
            fill = "lightblue", inherit.aes = FALSE) +
  geom_point(aes(color = point_color), size = 0.3, alpha = 0.8) +
  scale_color_identity() +
  geom_segment(data = cohort_info,
               aes(x = min_index, xend = max_index,
                   y = median_tmb, yend = median_tmb),
               color = "red",
               linewidth = 0.7,
               inherit.aes = FALSE) +
  scale_y_log10(
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1000"),
    limits = c(0.05, 1000)
  ) +
  scale_alpha_identity() +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "gray90", linetype = "dashed", linewidth = 0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "none"
  ) +
  labs(x = "", y = "TMB") +
  scale_x_continuous(
    breaks = cohort_info$mid_index,
    labels = cohort_info$cohort,
    expand = c(0, 0)
  )

solid_index <- which(cohort_info$cohort == "SOLID")
if(length(solid_index) > 0) {
  p <- p +
    annotate("text",
             x = cohort_info$mid_index[solid_index],
             y = 0.04,
             label = "SOLID",
             color = "red",
             fontface = "bold",
             size = 8 / 2.8,
             angle = 90,
             vjust = 0,
             hjust = 0.5)
}

print(p)
