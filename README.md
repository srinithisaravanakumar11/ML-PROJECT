# ClinVar Variant Pathogenicity Baseline

Binary classification (Pathogenic vs Benign) of ClinVar variants using **only
features available in the raw ClinVar export**, evaluated with a chromosome
hold-out, validation-tuned thresholds, multiple seeds and bootstrap confidence
intervals.

> **Scope and honesty note.** An earlier version of this project reported
> near-perfect scores (PR-AUC ≈ 0.9998) using gnomAD/CADD/SIFT/PolyPhen/phastCons
> values that were *simulated from the label*, plus random "sequence" features.
> Those results were invalid; see [docs/phase8_leakage_audit.md](docs/phase8_leakage_audit.md)
>.
> The current results use no simulated data. No real external annotations or
> reference sequence have been added yet, so this is a ClinVar-only baseline,
> not a clinically validated tool.

## Current results (test = chromosomes 21, 22, X, Y, MT)

1,415,482 train / 151,323 validation / 135,521 test variants. Mean over 3
seeds; threshold chosen on validation to maximise F1. Full tables (std, CIs,
calibration, per-type breakdown) are in [output/reports/baseline_report.md](output/reports/baseline_report.md).

| Model | ROC-AUC | PR-AUC | F1 | Precision | Recall |
|---|---|---|---|---|---|
| Logistic Regression | 0.812 | 0.683 | 0.630 | 0.705 | 0.569 |
| Random Forest | 0.808 | 0.668 | 0.633 | 0.727 | 0.561 |
| XGBoost | 0.802 | 0.681 | 0.636 | 0.742 | 0.556 |
| MLP | 0.784 | 0.674 | 0.624 | 0.738 | 0.540 |
| LightGBM | 0.808 | 0.684 | 0.637 | 0.730 | 0.566 |

Bootstrap 95% CIs are roughly ±0.003–0.005 wide, so the models are close and no
model is clearly best.

Findings worth knowing:

- **Performance is driven by variant type.** The test set is 86% SNVs (16.7%
  pathogenic), where PR-AUC is only 0.454 and F1 0.381; deletions/duplications are
  70–80% pathogenic, which lifts the aggregate numbers.
- **Gene features add little** (LightGBM ROC-AUC 0.810 → 0.802 without them),
  because only 0.1% of test rows have a gene seen in training. Also dropping
  `ReviewStatus` gives 0.796.
- Probabilities are poorly calibrated for Logistic Regression / Random Forest
  (ECE 0.12–0.14); XGBoost and LightGBM are better (ECE ≈ 0.04).

## Data-integrity checks

[output/reports/split_audit.md](output/reports/split_audit.md): no duplicate or
conflicting variant keys, zero variant-key overlap across train/val/test, no
unassigned chromosomes. Gene target encoding is out-of-fold for training rows
and training-only for validation/test (`src/gene_encoding.py`, unit-tested).

## Reproducing

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt   # Linux/macOS: .venv/bin/pip
# Download variant_summary.txt.gz from
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ into the project root
# and record the release date you used.
python src/preprocess_clinvar.py        # -> output/data/clinvar_ml_dataset.csv
python src/audit_splits.py              # -> output/reports/split_audit.md
python src/train_models.py              # ~2 h on CPU (3 seeds); smoke test: --seeds 1 --max-train 200000
python src/plot_model_evaluation_graphs.py
pytest
```

On Windows set `PYTHONUTF8=1` if the console encoding causes errors.

## Layout

| Path | Purpose |
|---|---|
| `src/paths.py` | All file locations |
| `src/preprocess_clinvar.py` | Clean ClinVar, binary label |
| `src/features.py` | Label-free feature engineering and chromosome split |
| `src/gene_encoding.py` | Leakage-safe gene features |
| `src/metrics.py` | Threshold selection, bootstrap CIs, calibration |
| `src/train_models.py` | Benchmark, ablation, report |
| `src/audit_splits.py` | Split-integrity audit |
| `src/join_annotations.py` | Join *real* external annotations you supply (never simulates) |
| `src/plot_model_evaluation_graphs.py` | Figures from `baseline_results.json` |
| `tests/` | pytest suite (CI in `.github/workflows`) |
| `docs/` | Historical leakage audit of the earlier, removed pipeline |

## Known gaps

1. **Real annotations** (gnomAD, CADD/SIFT/PolyPhen via dbNSFP, phastCons): put
   files in `data/annotations/`, run `src/join_annotations.py`, then extend the
   models to use them.
2. **Real sequence features / CNN**: requires the GRCh38 reference.
3. **External validation** on a later ClinVar release or independent cohort.
4. ClinVar labels depend on submitters and review status (label noise/provenance).
