# Sitting or standing? Recording posture and conversational movement

Data and analysis code for the article:

> Precloux, M., Guyot, P., & Marin, L. *Sitting or standing? Recording posture as a
> methodological parameter in the study of conversational movement.* Manuscript submitted
> for publication.

The study tests whether the posture in which a face-to-face conversation is recorded
(seated, semi-standing on a high stool, or standing) changes the amount of head and
upper-limb movement produced by interlocutors. Twenty dyads (N = 40 participants) held
three-minute conversations in each posture while their movements were recorded with an
eight-camera Vicon motion capture system (100 Hz), giving 120 participant × posture
recordings. Movement is quantified with the Quantity of Movement (QoM): the cumulative
three-dimensional displacement of the head (temple midpoint) and of both wrists,
normalized by each participant's shoulder width.

## Repository structure

```
posture-conversation-qom/
├── src/
│   ├── config.py                     # paths and processing parameters (Python)
│   ├── config_r.R                    # constants used by the R scripts
│   ├── 01_compute_qom_vicon.py       # raw Vicon CSV -> QoM per participant and posture
│   ├── 02_descriptive_stats.R        # Table 1
│   ├── 03_lmm_analysis.R             # main analysis: linear mixed models (Tables 2-3)
│   ├── 04_friedman_robustness.R      # robustness check: Friedman / Wilcoxon (Table 4)
│   ├── 05_dyad_cluster_bootstrap.R   # robustness check: dyad-level bootstrap (Table 4)
│   ├── 06_figures.R                  # Figure 2
│   └── 07_order_effect.R             # exploratory order effect (Section 3.4)
├── data/
│   ├── raw/                          # raw Vicon exports (not shared, git-ignored)
│   └── final/                        # analysis table, results and figures (versioned)
│       ├── posture_qom_long.xlsx     # QoM per participant and posture (120 rows)
│       ├── table1_descriptives.xlsx
│       ├── lmm_results.xlsx
│       ├── friedman_robustness_results.xlsx
│       ├── dyad_bootstrap_results.xlsx
│       ├── order_effect_results.xlsx
│       └── figures/
├── posture-conversation-qom.Rproj
├── requirements.txt
└── LICENSE
```

## Data availability

The raw Vicon marker trajectories and the videos are **not** shared, because they allow
participants to be identified (ethics approval: Research Ethics Committee of the
University of Montpellier, CER-UM, approval no. UM 2025-101bis).

The analysis table `data/final/posture_qom_long.xlsx` contains everything the statistical
analyses need, with no identifying data: dyad, posture, participant, median shoulder width,
and head and wrist QoM (in mm and in shoulder widths). All results reported in the article
can therefore be reproduced from this repository alone, starting at step 2 below.
`01_compute_qom_vicon.py` documents how this table was computed from the raw exports.

## Processing (`01_compute_qom_vicon.py`)

For each participant and recording:

1. Brief marker-label swaps are detected from aberrant segment lengths (shoulder–elbow or
   elbow–wrist length deviating by more than 25% from its median); the distal marker is set
   to missing over the affected samples (± 5 samples).
2. Gaps of up to 2 s are filled by linear interpolation.
3. Trajectories are low-pass filtered with a fourth-order Butterworth filter (6 Hz),
   applied forward and backward (zero phase lag).
4. The head is represented by the midpoint of the two temple markers, computed at each
   sample.
5. QoM is the cumulative sum of 3D Euclidean displacements between successive samples,
   summed over both wrists for the upper limbs.
6. QoM is divided by the participant's median shoulder width in the recording, so that it
   is dimensionless and expressed in shoulder widths.

Positions are expressed in the laboratory reference frame.

## Statistical analysis

All inferential analyses use log-transformed QoM (right-skewed distributions).

- **Main analysis** (`03_lmm_analysis.R`): linear mixed-effects models with posture as a
  fixed effect and random intercepts for dyad and for participant (nested within dyad).
  A combined model tests the Posture × Segment interaction
  (`log(QoM) ~ posture * segment + (1|dyad) + (1|subject)`); per-segment models estimate
  the effect of posture within each segment (`log(QoM) ~ posture + (1|dyad) + (1|subject)`).
  Type III F tests with Satterthwaite degrees of freedom; Bonferroni-corrected post hoc
  contrasts on estimated marginal means.
- **Robustness checks**: dyad-level cluster bootstrap (5,000 resamples, whole dyads
  resampled together; `05_dyad_cluster_bootstrap.R`) and nonparametric Friedman / Wilcoxon
  tests (`04_friedman_robustness.R`).
- **Order effect** (`07_order_effect.R`): Friedman tests with presentation rank as the
  within-subject factor. Postures were counterbalanced with a Latin square of three orders
  applied in turn to successive dyads (D01, D04, ...: seated, semi-standing, standing;
  D02, D05, ...: semi-standing, standing, seated; D03, D06, ...: standing, seated,
  semi-standing). The conversation scenarios were always presented in the same order, so
  presentation rank and scenario are confounded.

## Setup

Python ≥ 3.10 (step 1 only) and R ≥ 4.0.

```bash
pip install -r requirements.txt
Rscript -e 'install.packages(c("readxl","dplyr","tidyr","writexl","lme4","lmerTest","emmeans","ggplot2","here","patchwork"))'
```

## Reproducing the analysis

```bash
cd src
python3 01_compute_qom_vicon.py      # requires the raw exports in data/raw/vicon_csv/
Rscript 02_descriptive_stats.R
Rscript 03_lmm_analysis.R
Rscript 04_friedman_robustness.R
Rscript 05_dyad_cluster_bootstrap.R
Rscript 06_figures.R
Rscript 07_order_effect.R
```

Without the raw exports, skip step 1: the R scripts read the versioned
`data/final/posture_qom_long.xlsx`. Random seeds are fixed, so the bootstrap and Figure 2
are identical from one run to the next. The R scripts locate the repository root with the
`here` package (anchored by the `.Rproj` file).

## Main results

**Descriptive statistics (QoM in shoulder widths, n = 40 per cell)**

| Segment | Posture | Median [Q1–Q3] | Mean ± SD |
|---|---|---|---|
| Head | Seated | 21.9 [18.2–29.6] | 23.9 ± 9.7 |
| Head | Semi-standing | 22.1 [19.2–28.8] | 23.5 ± 8.0 |
| Head | Standing | 27.5 [19.9–34.1] | 27.6 ± 9.8 |
| Upper limbs | Seated | 61.3 [37.9–75.8] | 65.7 ± 51.7 |
| Upper limbs | Semi-standing | 58.2 [33.1–83.0] | 67.0 ± 47.0 |
| Upper limbs | Standing | 122.7 [82.1–154.5] | 125.3 ± 67.2 |

**Combined mixed model**

| Term | F | df | p |
|---|---|---|---|
| Posture | 34.85 | 2, 195 | < .001 |
| Segment | 403.75 | 1, 195 | < .001 |
| Posture × Segment | 14.76 | 2, 195 | < .001 |

Variance components: dyad 31.0%, participant 22.5%, residual 46.5%.

**Per-segment mixed models** (Bonferroni-corrected contrasts)

| Segment | Posture effect | Standing vs. seated | Standing vs. semi-standing | Semi-standing vs. seated |
|---|---|---|---|---|
| Head | F(2, 78) = 6.35, p = .0028 | p = .0079 | p = .0090 | p = 1.000 |
| Upper limbs | F(2, 78) = 41.81, p < .001 | p < .001 | p < .001 | p = 1.000 |

**Dyad-level bootstrap** (difference in mean log(QoM), 95% percentile CI)

| Segment | Contrast | Estimate | 95% CI |
|---|---|---|---|
| Head | Standing − Seated | 0.161 | [0.064, 0.270] |
| Head | Standing − Semi-standing | 0.158 | [0.077, 0.236] |
| Head | Semi-standing − Seated | 0.002 | [−0.119, 0.136] |
| Upper limbs | Standing − Seated | 0.763 | [0.605, 0.936] |
| Upper limbs | Standing − Semi-standing | 0.744 | [0.487, 1.024] |
| Upper limbs | Semi-standing − Seated | 0.019 | [−0.197, 0.226] |

**Friedman tests**: head, χ²(2) = 14.55, p = .0007, Kendall's W = 0.182; upper limbs,
χ²(2) = 36.15, p < .001, W = 0.452. Post hoc Wilcoxon tests: standing differs from both
seated and semi-standing for both segments (all Bonferroni-corrected p ≤ .006); seated and
semi-standing do not differ.

**Order effect**: head, χ²(2) = 0.45, p = .799; upper limbs, χ²(2) = 7.35, p = .025,
W = 0.092.

## License

Code released under the MIT License (see `LICENSE`).

## Contact

Matys Precloux — EuroMov Digital Health in Motion, University of Montpellier
(precloux.matys@gmail.com).
