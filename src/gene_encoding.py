"""Leakage-safe gene features shared by all training pipelines.

The previous implementation encoded every training row with a gene pathogenicity
rate that included the row's own label (target leakage, and a train/test
mismatch). Here training rows receive an out-of-fold (cross-fitted) encoding,
while validation/test rows are encoded from the full training set only.
"""
import numpy as np

SMOOTHING = 10
N_SPLITS = 5
SEED = 42


def _fit_encoding(genes, y, smoothing):
    """Smoothed per-gene positive rate; returns (mapping Series, global mean)."""
    global_mean = float(np.mean(y))
    counts = genes.value_counts()
    pos = genes[np.asarray(y) == 1].value_counts().reindex(counts.index, fill_value=0)
    return (pos + smoothing * global_mean) / (counts + smoothing), global_mean


def add_gene_features(df_train, df_val, df_test, label_col="y",
                      smoothing=SMOOTHING, n_splits=N_SPLITS, seed=SEED):
    """Add `gene_freq` and leakage-safe `gene_patho_rate` in place."""
    counts = df_train["GeneSymbol"].value_counts()
    for d in (df_train, df_val, df_test):
        d["gene_freq"] = d["GeneSymbol"].map(counts).fillna(0)

    # Training rows: each fold is encoded using labels from the other folds only.
    oof = np.empty(len(df_train), dtype=float)
    genes = df_train["GeneSymbol"].reset_index(drop=True)
    y = df_train[label_col].reset_index(drop=True)
    order = np.random.default_rng(seed).permutation(len(genes))
    for enc_idx in np.array_split(order, n_splits):
        fit_idx = np.setdiff1d(order, enc_idx, assume_unique=True)
        enc, gmean = _fit_encoding(genes.iloc[fit_idx], y.iloc[fit_idx], smoothing)
        oof[enc_idx] = genes.iloc[enc_idx].map(enc).fillna(gmean).to_numpy()
    df_train["gene_patho_rate"] = oof

    # Validation/test rows: encoded from the whole training set.
    enc, gmean = _fit_encoding(df_train["GeneSymbol"], df_train[label_col], smoothing)
    for d in (df_val, df_test):
        d["gene_patho_rate"] = d["GeneSymbol"].map(enc).fillna(gmean)
