import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from features import engineer, split  # noqa: E402
from gene_encoding import add_gene_features  # noqa: E402
from metrics import best_f1_threshold, expected_calibration_error  # noqa: E402


def _auc(y, s):
    r = pd.Series(s).rank().values
    p = y.sum()
    q = len(y) - p
    return (r[y == 1].sum() - p * (p + 1) / 2) / (p * q)


def test_gene_encoding_has_no_label_leak():
    rng = np.random.default_rng(0)
    n = 20000
    df = pd.DataFrame({"GeneSymbol": rng.integers(0, 4000, n).astype(str), "y": rng.integers(0, 2, n)})
    tr, va, te = df.iloc[:14000].copy(), df.iloc[14000:17000].copy(), df.iloc[17000:].copy()
    add_gene_features(tr, va, te)
    # labels are pure noise, so a leak-free encoding must not predict them
    assert abs(_auc(tr.y.values, tr.gene_patho_rate.values) - 0.5) < 0.05
    assert not tr.gene_patho_rate.isna().any() and not te.gene_patho_rate.isna().any()


def _toy():
    return pd.DataFrame({
        "Chromosome": ["1", "chr21", "X", "19", "5", "HLA-A"],
        "PositionVCF": [10, 20, 30, 40, 50, 60], "Stop": [10, 21, 30, 40, 50, 60],
        "ReferenceAlleleVCF": ["A", "AT", "C", "G", "T", "A"],
        "AlternateAlleleVCF": ["G", "A", "A", "C", "C", "T"],
        "target": ["Pathogenic", "Benign", "Pathogenic", "Benign", "Benign", "Benign"],
    })


def test_split_is_disjoint_and_excludes_unassigned():
    tr, va, te, unassigned = split(engineer(_toy()))
    assert set(tr.Chr_Str) == {"1", "5"} and set(va.Chr_Str) == {"19"} and set(te.Chr_Str) == {"21", "X"}
    assert unassigned.to_dict() == {"HLA-A": 1}


def test_engineer_flags_and_target():
    d = engineer(_toy())
    assert d.loc[1, "is_indel"] == 1 and d.loc[0, "is_transition"] == 1 and d.loc[3, "is_transition"] == 0
    assert d["y"].tolist() == [1, 0, 1, 0, 0, 0]


def test_threshold_chosen_from_given_labels_only():
    y = np.array([0, 0, 0, 1, 1, 1])
    p = np.array([0.1, 0.2, 0.3, 0.6, 0.7, 0.9])
    assert 0.3 < best_f1_threshold(y, p) <= 0.6


def test_ece_zero_for_perfectly_calibrated_extremes():
    y = np.array([0, 0, 1, 1])
    assert expected_calibration_error(y, np.array([0.0, 0.0, 1.0, 1.0])) == 0.0
