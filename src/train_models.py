"""Baseline ClinVar pathogenicity benchmark (leakage-audited).

Features come only from the raw ClinVar export (see features.py). Gene target
encoding is out-of-fold for training rows. For every model this script
reports, on the chromosome hold-out test set:

  * metrics at a threshold tuned on the VALIDATION set (not 0.5, not test),
  * mean/std across several random seeds,
  * 95% bootstrap confidence intervals,
  * calibration (Brier, ECE),
  * an ablation without gene-identity features,
  * a per-variant-type breakdown for the model chosen on validation PR-AUC.

Outputs: output/data/baseline_results.json, output/reports/baseline_report.md
Usage:   python src/train_models.py [--seeds 3] [--boot 200] [--max-train N]
"""
import argparse
import json
import time
import warnings
from datetime import datetime

import joblib
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from features import CAT_COLS, prepare
from metrics import best_f1_threshold, bootstrap_ci, point_metrics
from paths import DATA_DIR, ML_DATASET, MODEL_DIR, REPORTS_DIR

warnings.filterwarnings("ignore")

try:
    import lightgbm as lgb
except Exception:  # pragma: no cover - optional native dependency
    lgb = None

RESULTS_PATH = DATA_DIR / "baseline_results.json"
REPORT_PATH = REPORTS_DIR / "baseline_report.md"


def build_models(seed, pos_weight):
    models = {
        "Logistic Regression": LogisticRegression(class_weight="balanced", max_iter=500, random_state=seed),
        "Random Forest": RandomForestClassifier(
            n_estimators=200, max_depth=16, min_samples_leaf=5, class_weight="balanced",
            n_jobs=-1, random_state=seed),
        "XGBoost": xgb.XGBClassifier(
            n_estimators=200, max_depth=8, learning_rate=0.05, scale_pos_weight=pos_weight,
            n_jobs=-1, eval_metric="logloss", random_state=seed),
        "MLP": MLPClassifier(hidden_layer_sizes=(128, 64), max_iter=20, early_stopping=True, random_state=seed),
    }
    if lgb is not None:
        models["LightGBM"] = lgb.LGBMClassifier(
            n_estimators=300, learning_rate=0.05, num_leaves=63, class_weight="balanced",
            n_jobs=-1, random_state=seed, verbose=-1)
    else:
        models["LightGBM"] = HistGradientBoostingClassifier(
            max_iter=200, learning_rate=0.05, class_weight="balanced", random_state=seed)
    return models


def design_matrices(train, val, test, num_cols, cat_cols=CAT_COLS):
    pre = ColumnTransformer([
        ("num", StandardScaler(), num_cols),
        ("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False), cat_cols),
    ])
    return pre.fit_transform(train), pre.transform(val), pre.transform(test), pre


def evaluate(model, Xtr, ytr, Xva, yva, Xte, yte):
    t0 = time.time()
    model.fit(Xtr, ytr)
    fit_s = time.time() - t0
    p_val = model.predict_proba(Xva)[:, 1]
    p_test = model.predict_proba(Xte)[:, 1]
    thr = best_f1_threshold(yva, p_val)
    return {"fit_seconds": fit_s, "threshold": thr,
            "val": point_metrics(yva, p_val, thr), "test": point_metrics(yte, p_test, thr)}, p_test, thr


def summarise(runs):
    keys = runs[0]["test"].keys()
    return {k: {"mean": float(np.mean([r["test"][k] for r in runs])),
                "std": float(np.std([r["test"][k] for r in runs]))} for k in keys}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--max-train", type=int, default=0, help="subsample training rows (smoke test)")
    ap.add_argument("--save-models", action="store_true")
    args = ap.parse_args()

    t_start = time.time()
    if not ML_DATASET.exists():
        raise SystemExit(f"{ML_DATASET} not found - run src/preprocess_clinvar.py first.")
    df = pd.read_csv(ML_DATASET, low_memory=False)
    train, val, test, num_cols, unassigned = prepare(df, use_gene_features=True)
    if args.max_train:
        train = train.sample(min(args.max_train, len(train)), random_state=0)
    ytr, yva, yte = train["y"].values, val["y"].values, test["y"].values
    pos_weight = (ytr == 0).sum() / max(1, (ytr == 1).sum())
    print(f"train={len(train):,} val={len(val):,} test={len(test):,} unassigned-chrom rows excluded={int(unassigned.sum()):,}")

    Xtr, Xva, Xte, pre = design_matrices(train, val, test, num_cols)
    results = {"generated": datetime.now().isoformat(), "n": {"train": len(train), "val": len(val), "test": len(test)},
               "excluded_unassigned_chromosome_rows": int(unassigned.sum()),
               "seeds": args.seeds, "bootstrap_resamples": args.boot, "models": {}}

    best_name, best_val_pr, best_probs, best_thr = None, -1, None, None
    for name in build_models(0, pos_weight):
        runs, seed0_probs, seed0_thr = [], None, None
        for seed in range(args.seeds):
            model = build_models(seed, pos_weight)[name]
            r, p_test, thr = evaluate(model, Xtr, ytr, Xva, yva, Xte, yte)
            runs.append(r)
            if seed == 0:
                seed0_probs, seed0_thr = p_test, thr
                if args.save_models:
                    MODEL_DIR.mkdir(parents=True, exist_ok=True)
                    joblib.dump(model, MODEL_DIR / f"{name.lower().replace(' ', '_')}.joblib")
            print(f"{name:20s} seed={seed} test ROC={r['test']['ROC-AUC']:.4f} PR={r['test']['PR-AUC']:.4f} "
                  f"F1={r['test']['F1']:.4f} thr={thr:.3f} ({r['fit_seconds']:.0f}s)")
        results["models"][name] = {
            "test_mean_std": summarise(runs),
            "test_ci95_seed0": bootstrap_ci(yte, seed0_probs, seed0_thr, n_boot=args.boot),
            "val_pr_auc_seed0": runs[0]["val"]["PR-AUC"],
            "threshold_seed0": seed0_thr,
            "fit_seconds_mean": float(np.mean([r["fit_seconds"] for r in runs])),
            "per_seed": runs,
        }
        if runs[0]["val"]["PR-AUC"] > best_val_pr:
            best_name, best_val_pr, best_probs, best_thr = name, runs[0]["val"]["PR-AUC"], seed0_probs, seed0_thr

    # Per-variant-type breakdown for the model selected on VALIDATION PR-AUC.
    by_type = {}
    for t, idx in test.reset_index(drop=True).groupby("Type").groups.items():
        idx = np.asarray(list(idx))
        if len(idx) >= 200 and 0 < yte[idx].sum() < len(idx):
            by_type[str(t)] = {"n": int(len(idx)), "pathogenic_rate": float(yte[idx].mean()),
                               **{k: v for k, v in point_metrics(yte[idx], best_probs[idx], best_thr).items()
                                  if k in ("ROC-AUC", "PR-AUC", "F1")}}
    results["selected_model_by_val_pr_auc"] = best_name
    results["per_type_breakdown"] = by_type

    # Ablation: how much performance comes from shortcut features
    # (gene identity, ReviewStatus as a label-provenance proxy)?
    tr2, va2, te2, num2, _ = prepare(df, use_gene_features=False)
    if args.max_train:
        tr2 = tr2.sample(min(args.max_train, len(tr2)), random_state=0)
    no_review = [c for c in CAT_COLS if c != "ReviewStatus"]
    settings = {
        "all features": (train, val, test, num_cols, CAT_COLS),
        "without gene features": (tr2, va2, te2, num2, CAT_COLS),
        "without gene features and ReviewStatus": (tr2, va2, te2, num2, no_review),
    }
    abl = {}
    for label, (a_, b_, c_, nc, cc) in settings.items():
        Xa, Xb, Xc, _ = design_matrices(a_, b_, c_, nc, cc)
        r, _, _ = evaluate(build_models(0, pos_weight)["LightGBM"], Xa, a_["y"].values,
                           Xb, b_["y"].values, Xc, c_["y"].values)
        abl[label] = r["test"]
    results["shortcut_feature_ablation_lightgbm"] = abl

    results["runtime_minutes"] = (time.time() - t_start) / 60
    RESULTS_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_report(results)
    print(f"wrote {RESULTS_PATH} and {REPORT_PATH}")


def fmt(m, k):
    return f"{m['test_mean_std'][k]['mean']:.3f} ± {m['test_mean_std'][k]['std']:.3f}"


def write_report(res):
    L = ["# Baseline results (ClinVar-only features, chromosome hold-out)", "",
         f"Generated {res['generated']}. Train/val/test = {res['n']['train']:,}/{res['n']['val']:,}/{res['n']['test']:,} variants; "
         f"{res['seeds']} seeds; threshold tuned on validation (max F1); gene target encoding is out-of-fold.", "",
         "## Test metrics (mean ± std over seeds)", "",
         "| Model | ROC-AUC | PR-AUC | F1 | Precision | Recall | Bal. Acc. | Brier | ECE |", "|---|---|---|---|---|---|---|---|---|"]
    for n, m in res["models"].items():
        L.append(f"| {n} | " + " | ".join(fmt(m, k) for k in
                 ("ROC-AUC", "PR-AUC", "F1", "Precision", "Recall", "Balanced Accuracy", "Brier", "ECE")) + " |")
    L += ["", "## 95% bootstrap CIs (seed 0)", "", "| Model | ROC-AUC | PR-AUC | F1 |", "|---|---|---|---|"]
    for n, m in res["models"].items():
        c = m["test_ci95_seed0"]
        L.append(f"| {n} | " + " | ".join(f"[{c[k][0]:.3f}, {c[k][1]:.3f}]" for k in ("ROC-AUC", "PR-AUC", "F1")) + " |")
    L += ["", f"## Per-variant-type breakdown ({res['selected_model_by_val_pr_auc']}, chosen on validation PR-AUC)", "",
          "| Type | n | Pathogenic rate | ROC-AUC | PR-AUC | F1 |", "|---|---|---|---|---|---|"]
    for t, v in res["per_type_breakdown"].items():
        L.append(f"| {t} | {v['n']:,} | {v['pathogenic_rate']:.3f} | {v['ROC-AUC']:.3f} | {v['PR-AUC']:.3f} | {v['F1']:.3f} |")
    L += ["", "## Shortcut-feature ablation (LightGBM)", "", "| Setting | ROC-AUC | PR-AUC | F1 |", "|---|---|---|---|"]
    for k, v in res["shortcut_feature_ablation_lightgbm"].items():
        L.append(f"| {k} | {v['ROC-AUC']:.3f} | {v['PR-AUC']:.3f} | {v['F1']:.3f} |")
    L += ["", "## Limitations", "",
          "- Features are limited to ClinVar's own coordinates/alleles/type/gene; no external annotations are used.",
          "- ClinVar labels depend on submitter and review status; `ReviewStatus` is a proxy for label provenance and may inflate scores.",
          "- No external cohort or later ClinVar release has been evaluated yet.",
          "- Not a clinically validated tool."]
    REPORT_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
