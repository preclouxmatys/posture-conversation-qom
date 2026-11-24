# =========================================================
# 06_figures.R
# Figure 2: boxplots of head and upper-limb QoM by posture
# (statistical comparisons are reported in Tables 2-4, not on the figure).
# =========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
source(here("src", "config_r.R"))  # RANDOM_SEED
set.seed(RANDOM_SEED)  # reproducible jitter of the individual points

FINAL_DIR <- here("data", "final")
FIG_DIR <- here("data", "final", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read_excel(file.path(FINAL_DIR, "posture_qom_long.xlsx"))
df$posture <- factor(df$condition, levels = c("SEATED", "SEMI", "STANDING"),
                      labels = c("Seated", "Semi-standing", "Standing"))

plot_box <- function(data, value_col, title, ylab) {
  ggplot(data, aes(x = posture, y = .data[[value_col]], fill = posture)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.8, color = "grey30") +
    scale_fill_manual(values = c("Seated" = "#4C72B0", "Semi-standing" = "#DD8452", "Standing" = "#55A868")) +
    labs(title = title, x = "Posture", y = ylab) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))
}

p_head <- plot_box(df, "QoM_head_norm", "Head", "QoM (shoulder widths)")
p_wrist <- plot_box(df, "QoM_wrist_norm", "Upper limbs (wrists)", "QoM (shoulder widths)")

ggsave(file.path(FIG_DIR, "figure2_head_qom.png"), p_head, width = 5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "figure2_wrist_qom.png"), p_wrist, width = 5, height = 5, dpi = 300)

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- p_head + p_wrist
  ggsave(file.path(FIG_DIR, "figure2_combined.png"), combined, width = 10, height = 5, dpi = 300)
  cat("Combined figure saved: figure2_combined.png\n")
} else {
  cat("Package 'patchwork' not installed: figures saved separately ",
      "(figure2_head_qom.png, figure2_wrist_qom.png). ",
      "Install patchwork for the side-by-side version.\n")
}

cat("Figures saved in", FIG_DIR, "\n")
