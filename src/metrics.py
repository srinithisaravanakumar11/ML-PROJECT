"""Evaluation helpers: threshold selection, bootstrap CIs, calibration."""
import numpy as np
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             brier_score_loss, f1_score, precision_recall_curve,
                             precision_score, recall_score, roc_auc_score)


def best_f1_threshold(y_val, p_val):
    """Threshold maximising F1 on the *validation* set (never the test set)."""
    prec, rec, thr = precision_recall_curve(y_val, p_val)
    f1 = 2 * prec[:-1] * rec[:-1] / np.clip(prec[:-1] + rec[:-1], 1e-12, None)
    return float(thr[int(np.argmax(f1))])


def point_metrics(y, p, thr):
    pred = (p >= thr).astype(int)
    return {
        "ROC-AUC": roc_auc_score(y, p),
        "PR-AUC": average_precision_score(y, p),
        "F1": f1_score(y, pred),
        "Precision": precision_score(y, pred, zero_division=0),
        "Recall": recall_score(y, pred),
        "Balanced Accuracy": balanced_accuracy_score(y, pred),
        "Brier": brier_score_loss(y, np.clip(p, 0, 1)),
        "ECE": expected_calibration_error(y, p),
    }


def expected_calibration_error(y, p, bins=10):
    p = np.clip(p, 0, 1)
    idx = np.minimum((p * bins).astype(int), bins - 1)
    ece = 0.0
    for b in range(bins):
        m = idx == b
        if m.any():
            ece += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(ece)


def bootstrap_ci(y, p, thr, n_boot=200, seed=0):
    """95% percentile CI for each metric by resampling test rows."""
    rng = np.random.default_rng(seed)
    n = len(y)
    keys, rows = None, []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if y[i].min() == y[i].max():
            continue
        m = point_metrics(y[i], p[i], thr)
        keys = keys or list(m)
        rows.append([m[k] for k in keys])
    arr = np.array(rows)
    lo, hi = np.percentile(arr, [2.5, 97.5], axis=0)
    return {k: [float(a), float(b)] for k, a, b in zip(keys, lo, hi)}
