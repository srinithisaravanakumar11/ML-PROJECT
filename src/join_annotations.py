"""Join REAL external annotations onto the ClinVar dataset.

This project previously "simulated" gnomAD/CADD/SIFT/PolyPhen/phastCons values
from the label, which leaked the target (see docs/phase8_leakage_audit.md).
This script never generates values: it left-joins user-supplied annotation
tables on a normalised (chrom, pos, ref, alt) key and leaves unmatched rows
missing (NaN).

Annotation files go in data/annotations/*.tsv(.gz|.csv) and must contain the
columns: chrom, pos, ref, alt, plus any numeric annotation columns. Record the
source name/version of each file in data/annotations/SOURCES.md.

Output: output/data/clinvar_annotated_dataset.csv  (gitignored)
Usage:  python src/join_annotations.py
"""
import sys

import pandas as pd

from paths import DATA_DIR, ML_DATASET, ROOT

ANNOT_DIR = ROOT / "data" / "annotations"
OUT = DATA_DIR / "clinvar_annotated_dataset.csv"
KEY = ["chrom", "pos", "ref", "alt"]


def norm_key(chrom, pos, ref, alt):
    chrom = chrom.astype(str).str.replace("chr", "", case=False).str.upper().replace({"M": "MT"})
    return chrom + ":" + pos.astype(str) + ":" + ref.astype(str).str.upper() + ":" + alt.astype(str).str.upper()


def main():
    files = sorted(p for p in ANNOT_DIR.glob("*") if p.suffix in {".tsv", ".csv", ".gz"})
    if not files:
        sys.exit(f"No annotation files in {ANNOT_DIR}. Nothing joined; refusing to fabricate values.")
    df = pd.read_csv(ML_DATASET, low_memory=False)
    df["variant_key"] = norm_key(df["Chromosome"], df["PositionVCF"], df["ReferenceAlleleVCF"], df["AlternateAlleleVCF"])
    for f in files:
        a = pd.read_csv(f, sep=None, engine="python")
        missing = [c for c in KEY if c not in a.columns]
        if missing:
            sys.exit(f"{f.name}: missing key columns {missing}")
        a["variant_key"] = norm_key(a["chrom"], a["pos"], a["ref"], a["alt"])
        a = a.drop(columns=KEY).drop_duplicates("variant_key")
        df = df.merge(a, on="variant_key", how="left")
        cols = [c for c in a.columns if c != "variant_key"]
        print(f"{f.name}: columns={cols} matched={df[cols[0]].notna().mean():.1%}")
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
