# Baseline results (ClinVar-only features, chromosome hold-out)

Generated 2026-10-08T09:43:57.354998. Train/val/test = 1,415,482/151,323/135,521 variants; 3 seeds; threshold tuned on validation (max F1); gene target encoding is out-of-fold.

## Test metrics (mean ± std over seeds)

| Model | ROC-AUC | PR-AUC | F1 | Precision | Recall | Bal. Acc. | Brier | ECE |
|---|---|---|---|---|---|---|---|---|
| Logistic Regression | 0.812 ± 0.000 | 0.683 ± 0.000 | 0.630 ± 0.000 | 0.705 ± 0.000 | 0.569 ± 0.000 | 0.746 ± 0.000 | 0.140 ± 0.000 | 0.119 ± 0.000 |
| Random Forest | 0.808 ± 0.001 | 0.668 ± 0.002 | 0.633 ± 0.004 | 0.727 ± 0.013 | 0.561 ± 0.014 | 0.746 ± 0.004 | 0.147 ± 0.001 | 0.142 ± 0.004 |
| XGBoost | 0.802 ± 0.000 | 0.681 ± 0.000 | 0.636 ± 0.000 | 0.742 ± 0.000 | 0.556 ± 0.000 | 0.747 ± 0.000 | 0.127 ± 0.000 | 0.044 ± 0.000 |
| MLP | 0.784 ± 0.009 | 0.674 ± 0.007 | 0.624 ± 0.006 | 0.738 ± 0.015 | 0.540 ± 0.017 | 0.739 ± 0.005 | 0.134 ± 0.002 | 0.083 ± 0.016 |
| LightGBM | 0.808 ± 0.002 | 0.684 ± 0.001 | 0.637 ± 0.001 | 0.730 ± 0.009 | 0.566 ± 0.003 | 0.749 ± 0.000 | 0.126 ± 0.000 | 0.036 ± 0.002 |

## 95% bootstrap CIs (seed 0)

| Model | ROC-AUC | PR-AUC | F1 |
|---|---|---|---|
| Logistic Regression | [0.809, 0.815] | [0.678, 0.687] | [0.626, 0.633] |
| Random Forest | [0.806, 0.811] | [0.666, 0.674] | [0.623, 0.631] |
| XGBoost | [0.799, 0.805] | [0.676, 0.685] | [0.632, 0.640] |
| MLP | [0.794, 0.800] | [0.678, 0.686] | [0.628, 0.635] |
| LightGBM | [0.807, 0.814] | [0.680, 0.689] | [0.631, 0.639] |

## Per-variant-type breakdown (LightGBM, chosen on validation PR-AUC)

| Type | n | Pathogenic rate | ROC-AUC | PR-AUC | F1 |
|---|---|---|---|---|---|
| deletion | 10,705 | 0.781 | 0.669 | 0.880 | 0.882 |
| duplication | 4,782 | 0.709 | 0.638 | 0.825 | 0.835 |
| indel | 705 | 0.862 | 0.853 | 0.972 | 0.926 |
| insertion | 868 | 0.571 | 0.702 | 0.770 | 0.747 |
| microsatellite | 2,411 | 0.448 | 0.836 | 0.801 | 0.743 |
| single nucleotide variant | 115,950 | 0.167 | 0.714 | 0.454 | 0.381 |

## Shortcut-feature ablation (LightGBM)

| Setting | ROC-AUC | PR-AUC | F1 |
|---|---|---|---|
| all features | 0.810 | 0.685 | 0.636 |
| without gene features | 0.802 | 0.681 | 0.635 |
| without gene features and ReviewStatus | 0.796 | 0.678 | 0.634 |

## Limitations

- Features are limited to ClinVar's own coordinates/alleles/type/gene; no external annotations are used.
- ClinVar labels depend on submitter and review status; `ReviewStatus` is a proxy for label provenance and may inflate scores.
- No external cohort or later ClinVar release has been evaluated yet.
- Not a clinically validated tool.
