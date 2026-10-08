> **Historical document.** This audit describes the earlier, now-removed phase 2-7 pipeline (scripts and results are in git history before the commit that removed them). The current pipeline uses no simulated data.

# Phase 8: Feature Leakage and Robustness Audit

## Executive finding

**Status: CRITICAL — the enhanced-model evaluation is not a valid estimate of pathogenicity prediction from independent biological annotations.**

The source code derives `is_pathogenic` directly from `target` and then passes it into the functions that generate `gnomAD_AF`, `CADD_phred`, `SIFT_score`, `PolyPhen2_score`, and `phastCons_100way`. Those functions select different random distributions according to the target label. `gnomAD_AF_popmax` is then computed from the label-conditioned synthetic `gnomAD_AF`. All six fields are included in the Phase 5 predictor list. The chromosome split does not prevent this leakage because the feature values are generated using labels before partitioning.

This audit did not train models or change any previously recorded metric. Phase 5 results below are reproduced as reported only; they must not be interpreted as independent validation performance.

## Evidence from the feature-generation source

- `is_pathogenic` is assigned from `target` at [build_phase2_biological_features.py](../build_phase2_biological_features.py) lines 231.
- Label-dependent distributions use `np.where(is_pathogenic == 1, ...)` at lines 123, 141, 159, 164, 182.
- The source passes `is_pathogenic` into feature-generation calls for AF, CADD, SIFT/PolyPhen-2 and phastCons (gnomAD-like: line(s) 248; CADD-like: line(s) 253; SIFT/PolyPhen-like: line(s) 256; phastCons-like: line(s) 259); see [build_phase2_biological_features.py](../build_phase2_biological_features.py).
- The Phase 5 feature list explicitly includes those generated columns in [train_phase5_enhanced_models.py](../train_phase5_enhanced_models.py).
- Train-set gene target encoding also uses each training row's label to construct its own `gene_patho_rate` encoding at [train_phase5_enhanced_models.py](../train_phase5_enhanced_models.py) line(s) [113]. This creates training-time target leakage for those rows; use out-of-fold encodings or leave-one-out encodings if retaining this feature.
- The same train-label-based `gene_patho_rate` construction is present in the baseline training pipeline at [train_models.py](../train_models.py) line(s) [166]. Baseline metrics therefore also need to be rerun after correcting this train/test feature mismatch.
- The Phase 2 generator also creates sequence-context features from random distributions rather than reference-genome sequence. These may not use labels directly, but they are not authentic genomic annotations and cannot support biological claims.

## Feature provenance review

| Feature | Used in Phase 5 | Source assessment |
|---|:---:|---|
| `gnomAD_AF` | Yes | TARGET-CONDITIONED SYNTHETIC: estimate_gnomad_af_by_consequence() chooses pathogenic or benign distributions using is_pathogenic. |
| `gnomAD_AF_popmax` | Yes | DERIVED FROM TARGET-CONDITIONED SYNTHETIC FEATURE: Calculated from synthetic gnomAD_AF multiplied by a random factor. |
| `CADD_phred` | Yes | TARGET-CONDITIONED SYNTHETIC: estimate_cadd_score() chooses separate distributions using is_pathogenic. |
| `SIFT_score` | Yes | TARGET-CONDITIONED SYNTHETIC: estimate_sift_polyphen() chooses separate SIFT distributions using is_pathogenic. |
| `PolyPhen2_score` | Yes | TARGET-CONDITIONED SYNTHETIC: estimate_sift_polyphen() chooses separate PolyPhen-2 distributions using is_pathogenic. |
| `phastCons_100way` | Yes | TARGET-CONDITIONED SYNTHETIC: estimate_phastcons() chooses separate conservation distributions using is_pathogenic. |
| `gc_content` | Yes | SYNTHETIC; NOT AN EXTERNAL ANNOTATION: Generated from a random distribution; function does not query a reference genome. |
| `cpg_count` | Yes | SYNTHETIC; NOT AN EXTERNAL ANNOTATION: Generated from a random distribution parameterized by synthetic gc_content. |
| `homopolymer_max_len` | Yes | SYNTHETIC; NOT AN EXTERNAL ANNOTATION: Generated from a random geometric distribution, not computed from sequence. |

## Empirical audit of the existing enhanced dataset

- Rows scanned: 1,702,326
- Rows with an unsupported target label: 0
- Rows where `is_pathogenic` disagrees with `target`: 0
- Split is chromosome-based: train chromosomes 1–17 (plus any unlisted chromosomes, as in the training code), validation 18–20, test 21–22/X/Y/MT. See the counts below.
- Source code uses label-conditioned synthetic annotation values; different observed class means are diagnostic, but do not make those values authentic or leakage-free.

| Feature | Benign mean (n) | Pathogenic mean (n) | Label provenance |
|---|---:|---:|---|
| gnomAD_AF | 0.045330 (n=1,211,807) | 0.000033 (n=286,109) | TARGET-CONDITIONED SYNTHETIC |
| gnomAD_AF_popmax | 0.085452 (n=1,211,807) | 0.000075 (n=286,109) | DERIVED FROM TARGET-CONDITIONED SYNTHETIC FEATURE |
| CADD_phred | 11.367198 (n=1,308,716) | 28.496474 (n=309,192) | TARGET-CONDITIONED SYNTHETIC |
| SIFT_score | 0.800125 (n=1,304,573) | 0.091170 (n=176,478) | TARGET-CONDITIONED SYNTHETIC |
| PolyPhen2_score | 0.199843 (n=1,304,573) | 0.909504 (n=176,478) | TARGET-CONDITIONED SYNTHETIC |
| phastCons_100way | 0.418635 (n=1,266,649) | 0.833025 (n=299,054) | TARGET-CONDITIONED SYNTHETIC |
| gc_content | 0.440052 (n=1,377,069) | 0.440348 (n=325,257) | SYNTHETIC; NOT AN EXTERNAL ANNOTATION |
| cpg_count | 2.638081 (n=1,377,069) | 2.641628 (n=325,257) | SYNTHETIC; NOT AN EXTERNAL ANNOTATION |
| homopolymer_max_len | 3.496140 (n=1,377,069) | 3.497563 (n=325,257) | SYNTHETIC; NOT AN EXTERNAL ANNOTATION |

### Chromosome split counts

| Split | Benign | Pathogenic | Total |
|---|---:|---:|---:|
| Train | 1,145,043 | 270,439 | 1,415,482 |
| Validation | 129,846 | 21,477 | 151,323 |
| Test | 102,180 | 33,341 | 135,521 |

## Previously reported enhanced results (unchanged)

| Model | Test ROC-AUC | Test PR-AUC | Test F1 |
|---|---:|---:|---:|
| Logistic Regression (Balanced) | 0.9998 | 0.9995 | 0.9888 |
| Random Forest (Balanced) | 0.9997 | 0.9991 | 0.9915 |
| LightGBM (Balanced) | 0.9999 | 0.9998 | 0.9946 |
| XGBoost (Weighted) | 0.9999 | 0.9998 | 0.9942 |

These near-perfect scores (up to Test PR-AUC 0.9998) are materially compromised by the label-conditioned predictors identified above. They should be marked **invalid for biological generalization** and should not be presented as model performance based on real gnomAD, CADD, SIFT, PolyPhen-2, or phastCons data.

## Robustness and validity checklist

| Check | Finding | Interpretation |
|---|---|---|
| Direct target column among Phase 5 model inputs | `target` and `is_pathogenic` are not in the declared Phase 5 predictor lists | Direct columns are excluded, but derived label-conditioned predictors remain. |
| Biological annotation provenance | Confirmed synthetic and target-conditioned in source code | Critical leakage; existing enhanced test scores are invalid. |
| Chromosome holdout | Present in code | Useful only after removing leakage and using genuine, independent feature values. |
| Gene target encoding | Baseline and enhanced pipelines compute encodings from all training rows including each row's own label | Training-feature leakage and train/test mismatch; baseline and enhanced results need re-evaluation after cross-fitting or omission. |
| Independent external cohort | Not found in current project outputs | Generalization beyond this ClinVar-derived dataset remains untested. |
| Authentic biological annotations | Not established by current feature builder | Fetch/join real versioned sources using normalized GRCh38 variant identifiers and preserve missing values. |

## Required next experiment

1. Replace simulated label-conditioned annotations with authentic data-source joins. Do not infer or sample annotations from `target`.
2. Remove synthetic `gnomAD_AF`, `gnomAD_AF_popmax`, `CADD_phred`, `SIFT_score`, `PolyPhen2_score`, and `phastCons_100way` values from evaluation datasets until genuine annotations are available. Recompute missingness from the real annotation joins.
3. Either omit `gene_patho_rate` or generate training encodings out-of-fold; validation and test encodings must be derived only from training labels.
4. Keep the chromosome holdout, document any unassigned chromosomes, and audit exact and normalized variant-key overlap across splits.
5. Re-run baseline and enhanced models using the corrected data, then evaluate once on the untouched chromosome test set. For stronger evidence, add a temporally later ClinVar release or an independent external cohort with label provenance recorded.
6. Replace prior enhanced-model conclusions only with metrics produced by that corrected experiment. This audit deliberately leaves all prior result files unchanged.

## Audit artifacts

- Feature provenance: `docs/phase8_feature_provenance.csv`
- Class-conditional feature statistics: `docs/phase8_feature_label_statistics.csv`
- This report: `docs/phase8_leakage_robustness_report.md`
