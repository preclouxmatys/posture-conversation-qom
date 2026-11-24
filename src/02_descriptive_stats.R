# =========================================================
# 02_descriptive_stats.R
# Descriptive statistics of head and upper-limb QoM by posture (Table 1).
# =========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

FINAL_DIR <- here("data", "final")

df <- read_excel(file.path(FINAL_DIR, "posture_qom_long.xlsx"))
df$posture <- factor(df$condition, levels = c("SEATED", "SEMI", "STANDING"), ordered = TRUE)

desc_long <- df %>%
  select(posture, QoM_head_norm, QoM_wrist_norm) %>%
  pivot_longer(cols = c(QoM_head_norm, QoM_wrist_norm),
               names_to = "segment", values_to = "QoM") %>%
  mutate(segment = recode(segment, QoM_head_norm = "Head", QoM_wrist_norm = "Upper limbs"))

descriptives <- desc_long %>%
  group_by(segment, posture) %>%
  summarise(
    n = sum(!is.na(QoM)),
    median = median(QoM, na.rm = TRUE),
    Q1 = quantile(QoM, .25, na.rm = TRUE),
    Q3 = quantile(QoM, .75, na.rm = TRUE),
    mean = mean(QoM, na.rm = TRUE),
    sd = sd(QoM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(segment, posture)

print(descriptives)

write_xlsx(list(descriptives = descriptives),
           file.path(FINAL_DIR, "table1_descriptives.xlsx"))
cat("Saved:", file.path(FINAL_DIR, "table1_descriptives.xlsx"), "\n")
