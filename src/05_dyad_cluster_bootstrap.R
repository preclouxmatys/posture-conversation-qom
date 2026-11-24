# =========================================================
# 05_dyad_cluster_bootstrap.R
#
# Dyad-level cluster bootstrap (5,000 resamples) giving 95% CIs on the posture
# contrasts that respect the dependency structure of the data: the two members
# of a dyad are always resampled TOGETHER. Complements the LMM
# (03_lmm_analysis.R), which models the same dependency with random effects.
#
# Resampled statistic: difference in mean log(QoM) between pairs of postures
# (Standing-Seated, Standing-Semi, Semi-Seated), per segment (head, upper limbs).
# Vectorized implementation (base R matrices).
#
# Input : data/final/posture_qom_long.xlsx
# Output: data/final/dyad_bootstrap_results.xlsx
# =========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

source(here("src", "config_r.R"))  # RANDOM_SEED, N_BOOT

FINAL_DIR <- here("data", "final")

df <- read_excel(file.path(FINAL_DIR, "posture_qom_long.xlsx"))
df$posture <- factor(df$condition, levels = c("SEATED", "SEMI", "STANDING"))

dyads <- sort(unique(df$dyad))
n_dyads <- length(dyads)

# ---------------------------------------------------------
# For a given segment, builds a matrix (40 participants x 3 postures) of
# log(QoM), with a parallel vector giving the dyad of each row.
# ---------------------------------------------------------
build_subject_matrix <- function(data, value_col) {
  wide <- data %>%
    select(dyad, subject, posture, value = all_of(value_col)) %>%
    mutate(value = log(value)) %>%
    pivot_wider(names_from = posture, values_from = value) %>%
    arrange(dyad, subject)
  list(
    dyad = wide$dyad,
    mat = as.matrix(wide[, c("SEATED", "SEMI", "STANDING")])
  )
}

contrasts_from_colmeans <- function(cm) {
  c(
    Standing_minus_Seated = unname(cm["STANDING"] - cm["SEATED"]),
    Standing_minus_Semi   = unname(cm["STANDING"] - cm["SEMI"]),
    Semi_minus_Seated     = unname(cm["SEMI"] - cm["SEATED"])
  )
}

bootstrap_dyad_cluster <- function(data, value_col, label) {
  built <- build_subject_matrix(data, value_col)
  mat <- built$mat            # 40 x 3, log(QoM)
  subj_dyad <- built$dyad     # dyad of each row (length 40)

  # observed
  obs <- contrasts_from_colmeans(colMeans(mat))

  # row indices (participants) belonging to each dyad, computed once
  rows_by_dyad <- split(seq_len(nrow(mat)), subj_dyad)

  set.seed(RANDOM_SEED)
  boot_mat <- matrix(NA_real_, nrow = N_BOOT, ncol = 3,
                      dimnames = list(NULL, names(obs)))

  for (i in seq_len(N_BOOT)) {
    sampled_dyads <- sample(dyads, size = n_dyads, replace = TRUE)
    row_idx <- unlist(rows_by_dyad[sampled_dyads], use.names = FALSE)
    cm <- colMeans(mat[row_idx, , drop = FALSE])
    boot_mat[i, ] <- contrasts_from_colmeans(cm)
  }

  ci <- apply(boot_mat, 2, quantile, probs = c(.025, .975))
  se <- apply(boot_mat, 2, sd)

  out <- data.frame(
    segment = label,
    contrast = names(obs),
    estimate_logQoM_diff = unname(obs),
    boot_SE = unname(se),
    CI95_lower = ci[1, ],
    CI95_upper = ci[2, ]
  )
  out$pct_boot_same_sign <- sapply(seq_along(obs), function(k) {
    100 * mean(sign(boot_mat[, k]) == sign(obs[k]))
  })
  rownames(out) <- NULL
  out
}

cat("Dyad cluster bootstrap -", N_BOOT, "resamples,", n_dyads, "dyads\n")

t0 <- Sys.time()
res_head <- bootstrap_dyad_cluster(df, "QoM_head_norm", "Head")
res_wrist <- bootstrap_dyad_cluster(df, "QoM_wrist_norm", "Upper limbs")
cat("Elapsed time:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

cat("\n--- Head ---\n"); print(res_head)
cat("\n--- Upper limbs ---\n"); print(res_wrist)

all_results <- bind_rows(res_head, res_wrist)

write_xlsx(list(dyad_bootstrap = all_results),
           file.path(FINAL_DIR, "dyad_bootstrap_results.xlsx"))
cat("\nSaved:", file.path(FINAL_DIR, "dyad_bootstrap_results.xlsx"), "\n")
