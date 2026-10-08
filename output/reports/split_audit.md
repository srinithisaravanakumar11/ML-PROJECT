# Split audit

Total variants: 1,702,326. Train/val/test: 1,415,482/151,323/135,521.

| Check | Result |
|---|---|
| Rows on chromosomes in no split (excluded) | 0 {} |
| Duplicate variant keys in dataset | 0 |
| Variant keys with conflicting labels | 0 |
| Train∩val variant keys | 0 |
| Train∩test variant keys | 0 |
| Val∩test variant keys | 0 |
| Test genes also present in train | 5 of 1,557 |
| Test rows whose gene appears in train | 0.1% |

Gene overlap matters: genes spanning/being annotated on several chromosomes, or gene-level features (frequency, pathogenic rate), let a model memorise gene identity rather than learn variant-level signal. See the shortcut-feature ablation in `baseline_report.md`.
