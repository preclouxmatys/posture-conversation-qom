# =========================================================
# 03_lmm_analysis.R
# Main analysis: linear mixed-effects models (LMM).
#
# The 40 participants are nested within 20 dyads and are not independent
# observations. Random intercepts for dyad and for participant (nested within
# dyad); QoM is log-transformed (right-skewed, see Table 1).
#
# 1) Combined model (head and upper limbs stacked) with a Posture x Segment
#    interaction, testing directly whether the effect of posture differs
#    between segments.
# 2) Per-segment models (head, upper limbs), estimating the effect of posture
#    within each segment.
#
# Input : data/final/posture_qom_long.xlsx
# Output: data/final/lmm_results.xlsx
#         data/final/figures/lmm_diagnostics_combined.png
# =========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(writexl)

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

FINAL_DIR <- here("data", "final")

df <- read_excel(file.path(FINAL_DIR, "posture_qom_long.xlsx"))
df$posture <- factor(df$condition, levels = c("SEATED", "SEMI", "STANDING"))
df$dyad <- factor(df$dyad)
df$subject <- factor(df$subject)

# ---------------------------------------------------------
# 1) Combined model: log(QoM) ~ posture * segment + (1|dyad) + (1|subject)
# ---------------------------------------------------------
long <- df %>%
  select(dyad, subject, posture, QoM_head_norm, QoM_wrist_norm) %>%
  pivot_longer(cols = c(QoM_head_norm, QoM_wrist_norm),
               names_to = "segment", values_to = "QoM") %>%
  mutate(
    segment = recode(segment, QoM_head_norm = "Head", QoM_wrist_norm = "Wrist"),
    segment = factor(segment, levels = c("Head", "Wrist")),
    logQoM = log(QoM)
  )

cat("\n=== Combined model (Posture x Segment interaction) ===\n")
m_combined <- lmer(
  logQoM ~ posture * segment + (1 | dyad) + (1 | subject),
  data = long, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
print(summary(m_combined))

anova_combined <- anova(m_combined)  # Type III, Satterthwaite df (lmerTest)
print(anova_combined)

# Variance components (share of variance attributable to dyad membership)
vc <- as.data.frame(VarCorr(m_combined))
vc$pct_variance <- 100 * vc$vcov / sum(vc$vcov)
cat("\nVariance components:\n")
print(vc)

# ---------------------------------------------------------
# 2) Post hoc contrasts
#    - If the interaction is significant: simple effects of posture, BY segment
#    - Otherwise: main effect of posture, across segments
# ---------------------------------------------------------
interaction_p <- anova_combined["posture:segment", "Pr(>F)"]
cat(sprintf("\nPosture x Segment interaction: p = %.4f\n", interaction_p))

emm_by_segment <- emmeans(m_combined, ~ posture | segment)
posthoc_by_segment <- pairs(emm_by_segment, adjust = "bonferroni")
cat("\nPost hoc contrasts (posture, by segment):\n")
print(posthoc_by_segment)

emm_main <- emmeans(m_combined, ~ posture)
posthoc_main <- pairs(emm_main, adjust = "bonferroni")
cat("\nPost hoc contrasts (posture, main effect, across segments):\n")
print(posthoc_main)

# ---------------------------------------------------------
# 3) Per-segment models (Table 3)
# ---------------------------------------------------------
fit_segment_model <- function(data, label) {
  cat(sprintf("\n=== Model %s : log(QoM) ~ posture + (1|dyad) + (1|subject) ===\n", label))
  m <- lmer(log(QoM) ~ posture + (1 | dyad) + (1 | subject), data = data,
            control = lmerControl(optimizer = "bobyqa"))
  a <- anova(m)
  print(a)
  emm <- emmeans(m, ~posture)
  ph <- pairs(emm, adjust = "bonferroni")
  print(ph)
  list(model = m, anova = a, posthoc = ph)
}

head_data <- df %>% transmute(dyad, subject, posture, QoM = QoM_head_norm)
wrist_data <- df %>% transmute(dyad, subject, posture, QoM = QoM_wrist_norm)

res_head <- fit_segment_model(head_data, "Head")
res_wrist <- fit_segment_model(wrist_data, "Upper limbs")

# ---------------------------------------------------------
# 4) Export
# ---------------------------------------------------------
anova_combined_df <- as.data.frame(anova_combined) %>%
  tibble::rownames_to_column("term")

posthoc_by_segment_df <- as.data.frame(posthoc_by_segment)
posthoc_main_df <- as.data.frame(posthoc_main)

anova_head_df <- as.data.frame(res_head$anova) %>% tibble::rownames_to_column("term")
anova_wrist_df <- as.data.frame(res_wrist$anova) %>% tibble::rownames_to_column("term")
posthoc_head_df <- as.data.frame(res_head$posthoc)
posthoc_wrist_df <- as.data.frame(res_wrist$posthoc)

write_xlsx(
  list(
    combined_anova = anova_combined_df,
    combined_posthoc_by_segment = posthoc_by_segment_df,
    combined_posthoc_main = posthoc_main_df,
    variance_components = vc,
    head_anova = anova_head_df,
    head_posthoc = posthoc_head_df,
    wrist_anova = anova_wrist_df,
    wrist_posthoc = posthoc_wrist_df
  ),
  file.path(FINAL_DIR, "lmm_results.xlsx")
)
cat("\nSaved:", file.path(FINAL_DIR, "lmm_results.xlsx"), "\n")

# ---------------------------------------------------------
# 5) Diagnostics (normality of residuals, homoscedasticity)
# ---------------------------------------------------------
png(file.path(FINAL_DIR, "figures", "lmm_diagnostics_combined.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))
qqnorm(resid(m_combined), main = "Q-Q plot of residuals (combined model)")
qqline(resid(m_combined))
plot(fitted(m_combined), resid(m_combined),
     main = "Residuals vs fitted values", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)
dev.off()
cat("Diagnostics saved: data/final/figures/lmm_diagnostics_combined.png\n")
