"""Label-free feature engineering and the chromosome split, shared by all
training / auditing scripts. Only columns that exist in the raw ClinVar
export are used; nothing here is simulated or derived from the target."""
import numpy as np
import pandas as pd

from gene_encoding import add_gene_features

TRAIN_CHRS = [str(i) for i in range(1, 18)]
VAL_CHRS = ["18", "19", "20"]
TEST_CHRS = ["21", "22", "X", "Y", "MT", "M"]

NUM_BASE = ["PositionVCF", "Stop", "ref_len", "alt_len", "var_len", "log_var_len", "len_diff"]
NUM_GENE = ["gene_freq", "gene_patho_rate"]
CAT_COLS = ["Type", "Origin", "ReviewStatus", "is_transition", "is_indel"]

_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def _is_transition(ref, alt):
    if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
        return 1 if (ref, alt) in _TRANSITIONS else 0
    return -1


def engineer(df):
    ref = df["ReferenceAlleleVCF"].astype(str).str.upper()
    alt = df["AlternateAlleleVCF"].astype(str).str.upper()
    df["y"] = (df["target"].astype(str).str.strip().str.lower() == "pathogenic").astype(int)
    df["ref_len"] = ref.str.len()
    df["alt_len"] = alt.str.len()
    df["var_len"] = (df["Stop"] - df["PositionVCF"]).abs()
    df["log_var_len"] = np.log1p(df["var_len"])
    df["len_diff"] = df["alt_len"] - df["ref_len"]
    df["is_indel"] = ((df["ref_len"] != 1) | (df["alt_len"] != 1)).astype(int)
    df["is_transition"] = [_is_transition(r, a) for r, a in zip(ref, alt)]
    df["Chr_Str"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False).str.upper()
    df["variant_key"] = (df["Chr_Str"] + ":" + df["PositionVCF"].astype(str) + ":" + ref + ":" + alt)
    return df


def split(df):
    """Chromosome hold-out. Chromosomes outside every list are *excluded*
    (not silently added to train) and reported by the caller."""
    tr = df["Chr_Str"].isin(TRAIN_CHRS)
    va = df["Chr_Str"].isin(VAL_CHRS)
    te = df["Chr_Str"].isin(TEST_CHRS)
    unassigned = df.loc[~(tr | va | te), "Chr_Str"].value_counts()
    return df[tr].copy(), df[va].copy(), df[te].copy(), unassigned


def impute(train, *others, num_cols, cat_cols=CAT_COLS):
    """Median/mode imputation fitted on train only."""
    for c in num_cols:
        med = train[c].median()
        for d in (train, *others):
            d[c] = d[c].fillna(med)
    for c in cat_cols:
        mode = train[c].mode()[0]
        for d in (train, *others):
            d[c] = d[c].fillna(mode).astype(str)


def prepare(df, use_gene_features=True):
    """Return train/val/test frames with features added + the numeric column list."""
    df = engineer(df)
    train, val, test, unassigned = split(df)
    add_gene_features(train, val, test)
    num_cols = NUM_BASE + (NUM_GENE if use_gene_features else [])
    impute(train, val, test, num_cols=num_cols)
    return train, val, test, num_cols, unassigned
