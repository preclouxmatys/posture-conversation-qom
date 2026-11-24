# =========================================================
# 04_friedman_robustness.R
# Robustness check: nonparametric Friedman / Wilcoxon tests.
#
# These tests do not model clustering within dyads; they are reported
# alongside the main analysis (03_lmm_analysis.R) as a robustness check.
# =========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

FINAL_DIR <- here("data", "final")

df <- read_excel(file.path(FINAL_DIR, "posture_qom_long.xlsx"))
df$posture <- factor(df$condition, levels = c("SEATED", "SEMI", "STANDING"))

run_friedman <- function(data, value_col, label) {
  wide <- data %>%
    select(subject, posture, value = all_of(value_col)) %>%
    pivot_wider(names_from = posture, values_from = value)

  mat <- as.matrix(wide[, c("SEATED", "SEMI", "STANDING")])
  ft <- friedman.test(mat)

  # Kendall's W = chi2 / (n * (k-1))
  n <- nrow(mat)
  k <- ncol(mat)
  kendall_w <- unname(ft$statistic) / (n * (k - 1))

  cat(sprintf("\n=== Friedman: %s (n=%d) ===\n", label, n))
  print(ft)
  cat(sprintf("Kendall's W = %.3f\n", kendall_w))

  pairs_list <- list(c("SEATED", "SEMI"), c("SEATED", "STANDING"), c("SEMI", "STANDING"))
  posthoc <- lapply(pairs_list, function(p) {
    x <- mat[, p[1]]; y <- mat[, p[2]]
    wt <- wilcox.test(x, y, paired = TRUE, exact = FALSE)
    data.frame(comparison = paste(p[1], "vs", p[2]),
               V = unname(wt$statistic),
               p_raw = wt$p.value)
  }) %>% bind_rows() %>%
    mutate(p_bonferroni = pmin(p_raw * 3, 1))

  print(posthoc)

  list(
    summary = data.frame(label = label, n = n, chi2 = unname(ft$statistic),
                          df = unname(ft$parameter), p = ft$p.value, kendall_w = kendall_w),
    posthoc = posthoc
  )
}

res_head <- run_friedman(df, "QoM_head_norm", "Head")
res_wrist <- run_friedman(df, "QoM_wrist_norm", "Upper limbs")

write_xlsx(
  list(
    friedman_summary = bind_rows(res_head$summary, res_wrist$summary),
    posthoc_head = res_head$posthoc,
    posthoc_wrist = res_wrist$posthoc
  ),
  file.path(FINAL_DIR, "friedman_robustness_results.xlsx")
)
cat("\nSaved:", file.path(FINAL_DIR, "friedman_robustness_results.xlsx"), "\n")
