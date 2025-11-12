# Extended Data Figure 1

# Extended Data Figure 1A
# No code - no new code was used in the generation of this figure

# Extended Data Figure 1B
# No code - no new code was used in the generation of this figure

# Extended Data Figure 1C
# No code - no new code was used in the generation of this figure

# Extended Data Figure 1D: gammaH2AX Positive Cells
library(ggplot2)
library(dplyr)

ihc <- read.csv('data/ihc_data.csv')

t_test_result <- t.test(
  X..pH2AX.Positive.Cells ~ pre_or_post,
  data = ihc,
  var.equal = FALSE
)

p_value <- t_test_result$p.value

ihc_summary <- ihc %>%
  group_by(pre_or_post) %>%
  summarize(
    mean_value = mean(X..pH2AX.Positive.Cells, na.rm = TRUE),
    se_value = sd(X..pH2AX.Positive.Cells, na.rm = TRUE) / sqrt(n())
  )

ggplot(ihc, aes(x = factor(pre_or_post, levels = c("pre", "post")), y = X..pH2AX.Positive.Cells, fill = pre_or_post)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = pre_or_post), width = 0.2, size = 2, alpha = 0.7) +
  labs(
    title = paste0("gammaH2AX Positive Cells (p = ", signif(p_value, 3), ")"),
    x = "Treatment Group",
    y = "Percent gammaH2AX Positive Cells"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("pre" = "#56B4E9", "post" = "#E69F00")) +
  scale_color_manual(values = c("pre" = "#56B4E9", "post" = "#E69F00")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none"
  )

# Extended Data Figure 1E
# No code - no new code was used in the generation of this figure

# Extended Data Figure 1F
# No code - no new code was used in the generation of this figure

# Extended Data Figure 1G
# TBD - YE to add
