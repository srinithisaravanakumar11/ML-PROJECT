"""Audit the chromosome split for leakage-style problems.

Checks: variant-key overlap across splits, duplicate keys within the data,
conflicting labels for identical keys, gene overlap between train and test,
and rows on chromosomes that belong to no split.

Output: output/reports/split_audit.md
Usage:  python src/audit_splits.py
"""
import pandas as pd

from features import engineer, split
from paths import ML_DATASET, REPORTS_DIR


def main():
    df = engineer(pd.read_csv(ML_DATASET, low_memory=False))
    train, val, test, unassigned = split(df)
    keys = {n: set(d["variant_key"]) for n, d in (("train", train), ("val", val), ("test", test))}
    dup = int(df["variant_key"].duplicated().sum())
    conflicts = int((df.groupby("variant_key")["y"].nunique() > 1).sum())
    genes = {n: set(d["GeneSymbol"].dropna()) for n, d in (("train", train), ("val", val), ("test", test))}
    test_rows_known_gene = float(test["GeneSymbol"].isin(genes["train"]).mean())

    L = ["# Split audit", "",
         f"Total variants: {len(df):,}. Train/val/test: {len(train):,}/{len(val):,}/{len(test):,}.", "",
         "| Check | Result |", "|---|---|",
         f"| Rows on chromosomes in no split (excluded) | {int(unassigned.sum()):,} {dict(unassigned)} |",
         f"| Duplicate variant keys in dataset | {dup:,} |",
         f"| Variant keys with conflicting labels | {conflicts:,} |",
         f"| Train∩val variant keys | {len(keys['train'] & keys['val']):,} |",
         f"| Train∩test variant keys | {len(keys['train'] & keys['test']):,} |",
         f"| Val∩test variant keys | {len(keys['val'] & keys['test']):,} |",
         f"| Test genes also present in train | {len(genes['test'] & genes['train']):,} of {len(genes['test']):,} |",
         f"| Test rows whose gene appears in train | {test_rows_known_gene:.1%} |", "",
         "Gene overlap matters: genes spanning/being annotated on several chromosomes, or "
         "gene-level features (frequency, pathogenic rate), let a model memorise gene identity "
         "rather than learn variant-level signal. See the shortcut-feature ablation in `baseline_report.md`."]
    out = REPORTS_DIR / "split_audit.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
