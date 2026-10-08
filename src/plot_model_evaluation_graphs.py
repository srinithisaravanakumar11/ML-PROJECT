#!/usr/bin/env python3
"""Create publication-ready ClinVar model comparison visualizations."""

import json
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"
FIGURES_DIR = OUTPUT_DIR / "figures"
DATA_DIR = OUTPUT_DIR / "data"
MODELS = [
    "Logistic Regression",
    "Random Forest",
    "XGBoost",
    "MLP",
    "LightGBM",
]
COLORS = {
    "Logistic Regression": "#315f8c",
    "Random Forest": "#36836f",
    "XGBoost": "#c77b28",
    "MLP": "#8366a5",
    "LightGBM": "#b95151",
}
def load_results() -> Dict[str, Dict[str, float]]:
    """Read mean test metrics from output/data/baseline_results.json (written by train_models.py)."""
    path = DATA_DIR / "baseline_results.json"
    if not path.exists():
        raise SystemExit(f"{path} not found - run src/train_models.py first.")
    saved = json.loads(path.read_text(encoding="utf-8"))["models"]
    out = {}
    for name in MODELS:
        m = saved[name]
        out[name] = {k: m["test_mean_std"][k]["mean"] for k in
                     ("ROC-AUC", "PR-AUC", "F1", "Precision", "Recall", "Balanced Accuracy")}
        out[name]["Training Time (seconds)"] = m["fit_seconds_mean"]
    return out


RESULTS: Dict[str, Dict[str, float]] = load_results()
METRICS = [
    "ROC-AUC",
    "PR-AUC",
    "F1",
    "Precision",
    "Recall",
    "Balanced Accuracy",
]


def save_figure(
    fig: plt.Figure,
    filename: str,
    aliases: Optional[List[str]] = None,
) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for output_name in [filename] + (aliases or []):
        fig.savefig(FIGURES_DIR / output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def ranked(metric: str, descending: bool = True) -> List[str]:
    return sorted(
        MODELS,
        key=lambda model: RESULTS[model][metric],
        reverse=descending,
    )


def pareto_efficient(x_metric: str, y_metric: str) -> List[str]:
    """Return models not dominated by a faster, higher-scoring model."""
    efficient = []
    for model in MODELS:
        x_value = RESULTS[model][x_metric]
        y_value = RESULTS[model][y_metric]
        dominated = any(
            RESULTS[other][x_metric] <= x_value
            and RESULTS[other][y_metric] >= y_value
            and (
                RESULTS[other][x_metric] < x_value
                or RESULTS[other][y_metric] > y_value
            )
            for other in MODELS
            if other != model
        )
        if not dominated:
            efficient.append(model)
    return sorted(efficient, key=lambda model: RESULTS[model][x_metric])


def make_roc_radar() -> None:
    angles = np.linspace(0, 2 * np.pi, len(MODELS), endpoint=False)
    values = [RESULTS[model]["ROC-AUC"] for model in MODELS]
    closed_angles = np.append(angles, angles[0])
    closed_values = values + [values[0]]

    fig, ax = plt.subplots(figsize=(8.0, 7.2), subplot_kw={"polar": True})
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_xticks(angles, MODELS)
    ax.set_rlabel_position(20)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], color="#5a6570")
    ax.plot(closed_angles, closed_values, color="#315f8c", linewidth=2.2)
    ax.fill(closed_angles, closed_values, color="#4d83ad", alpha=0.14)
    ax.scatter(angles, values, color="#315f8c", s=52, zorder=3)
    best_index = int(np.argmax(values))
    ax.scatter(
        [angles[best_index]],
        [values[best_index]],
        color="#d04a3a",
        s=105,
        marker="*",
        zorder=4,
        label=f"Highest ROC-AUC: {MODELS[best_index]}",
    )
    for angle, value in zip(angles, values):
        ax.annotate(
            f"{value:.6f}",
            (angle, value),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )
    ax.set_title("ROC-AUC Model Profiles", pad=28, fontweight="semibold")
    ax.grid(color="#cdd5dd", linewidth=0.8)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.17), frameon=False)
    fig.text(
        0.5,
        0.015,
        "Radial scale: 0 to 1 (larger radius indicates higher ROC-AUC).",
        ha="center",
        fontsize=9,
        color="#55616c",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.98))
    save_figure(fig, "roc_auc_radar.png", ["roc_auc_comparison.png"])


def make_pr_auc_dotplot() -> None:
    models = ranked("PR-AUC")
    values = [RESULTS[model]["PR-AUC"] for model in models]
    best_model = models[0]
    fig, ax = plt.subplots(figsize=(8.6, 5.3))
    positions = np.arange(len(models))
    for position, model, value in zip(positions, models, values):
        ax.hlines(position, 0.65, value, color="#d8e0e7", linewidth=1.4, zorder=1)
        ax.scatter(
            value,
            position,
            color=COLORS[model],
            edgecolor="white",
            linewidth=0.8,
            s=90,
            zorder=3,
        )
        ax.annotate(
            f"{value:.6f}",
            (value, position),
            xytext=(9, 0),
            textcoords="offset points",
            va="center",
            fontsize=9,
        )
    best_value = RESULTS[best_model]["PR-AUC"]
    ax.axvline(best_value, color="#bb493c", linestyle=(0, (4, 3)), linewidth=1.2)
    ax.set_yticks(positions, models)
    ax.invert_yaxis()
    ax.set_xlim(0.65, 0.71)
    ax.set_xlabel("PR-AUC")
    ax.set_ylabel("Model (ranked)")
    ax.set_title("PR-AUC Ranked Dot Plot", loc="left", pad=14, fontweight="semibold")
    ax.text(
        0.99,
        0.02,
        f"Best reference: {best_model} ({best_value:.6f})",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color="#9c3d34",
        fontsize=9,
    )
    ax.grid(axis="x", color="#dce2e8", linewidth=0.8)
    ax.set_axisbelow(True)
    fig.text(
        0.125,
        0.025,
        "PR-AUC is particularly useful when evaluating performance under class imbalance.",
        fontsize=9,
        color="#55616c",
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save_figure(fig, "pr_auc_cleveland_dotplot.png", ["pr_auc_comparison.png"])


def make_f1_dumbbell() -> None:
    models = ranked("F1")
    positions = np.arange(len(models))
    fig, ax = plt.subplots(figsize=(10.2, 5.7))
    for position, model in zip(positions, models):
        precision = RESULTS[model]["Precision"]
        f1 = RESULTS[model]["F1"]
        recall = RESULTS[model]["Recall"]
        left, right = sorted([precision, recall])
        ax.hlines(position, left, right, color="#b9c4ce", linewidth=2.0, zorder=1)
        ax.scatter(
            precision,
            position,
            marker="o",
            s=65,
            color=COLORS[model],
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        ax.scatter(
            f1,
            position,
            marker="D",
            s=78,
            color="#273746",
            edgecolor="white",
            linewidth=0.8,
            zorder=4,
        )
        ax.scatter(
            recall,
            position,
            marker="s",
            s=58,
            color=COLORS[model],
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        ax.text(
            1.02,
            position,
            f"P {precision:.6f}   F1 {f1:.6f}   R {recall:.6f}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=8.5,
            color="#27313a",
        )
    ax.set_yticks(positions, models)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Score")
    ax.set_ylabel("Model (ranked by F1)")
    ax.set_title(
        "Precision–F1–Recall Balance by Model",
        loc="left",
        pad=14,
        fontweight="semibold",
    )
    ax.grid(axis="x", color="#dce2e8", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.scatter([], [], marker="o", color="#718096", label="Precision")
    ax.scatter([], [], marker="D", color="#273746", label="F1")
    ax.scatter([], [], marker="s", color="#718096", label="Recall")
    ax.legend(loc="lower left", ncol=3, frameon=False, bbox_to_anchor=(0, -0.23))
    fig.text(
        0.125,
        0.015,
        "F1 is the harmonic mean of precision and recall; its marker lies between the two endpoint scores.",
        fontsize=9,
        color="#55616c",
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save_figure(
        fig,
        "f1_precision_recall_dumbbell.png",
        ["f1_score_comparison.png"],
    )


def make_precision_beeswarm() -> None:
    models = ranked("Precision")
    values = [RESULTS[model]["Precision"] for model in models]
    offsets = np.array([0.10, -0.10, 0.10, -0.10, 0.0])
    positions = np.arange(len(models), dtype=float) + offsets
    fig, ax = plt.subplots(figsize=(8.4, 5.3))
    for model, value, position in zip(models, values, positions):
        ax.scatter(
            value,
            position,
            s=105,
            color=COLORS[model],
            edgecolor="white",
            linewidth=0.9,
            zorder=3,
        )
        ax.annotate(
            f"{value:.6f}",
            (value, position),
            xytext=(9, 0),
            textcoords="offset points",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(np.arange(len(models)), models)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Precision")
    ax.set_ylabel("Model (ranked)")
    ax.set_title("Precision by Model — Jittered Dot Plot", loc="left", pad=14, fontweight="semibold")
    ax.grid(axis="x", color="#dce2e8", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.text(
        0.99,
        0.02,
        "Vertical jitter is visual only; score coordinates are unchanged.",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        color="#55616c",
    )
    fig.tight_layout()
    save_figure(
        fig,
        "precision_jittered_dotplot.png",
        ["precision_comparison.png"],
    )


def make_recall_lollipop() -> None:
    models = ranked("Recall")
    positions = np.arange(len(models))
    fig, ax = plt.subplots(figsize=(8.5, 5.3))
    for position, model in zip(positions, models):
        recall = RESULTS[model]["Recall"]
        ax.hlines(position, 0, recall, color="#cbd5df", linewidth=2)
        ax.scatter(
            recall,
            position,
            s=95,
            color=COLORS[model],
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        ax.annotate(
            f"{recall:.6f}",
            (recall, position),
            xytext=(8, 0),
            textcoords="offset points",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(positions, models)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.7)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Model (ranked)")
    ax.set_title("Pathogenic Variant Detection Rate", loc="left", pad=14, fontweight="semibold")
    ax.grid(axis="x", color="#dce2e8", linewidth=0.8)
    ax.set_axisbelow(True)
    fig.text(
        0.125,
        0.025,
        "Recall is the proportion of actual pathogenic variants correctly identified by the model.",
        fontsize=9,
        color="#55616c",
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save_figure(
        fig,
        "recall_detection_lollipop.png",
        ["recall_comparison.png"],
    )


def make_balanced_accuracy_ranking() -> None:
    models = ranked("Balanced Accuracy")
    positions = np.arange(len(models))
    fig, ax = plt.subplots(figsize=(8.6, 5.3))
    for rank, (position, model) in enumerate(zip(positions, models), start=1):
        value = RESULTS[model]["Balanced Accuracy"]
        ax.scatter(
            rank,
            position,
            s=135,
            color=COLORS[model],
            edgecolor="white",
            linewidth=0.9,
            zorder=3,
        )
        ax.text(
            5.5,
            position,
            f"{model}   {value:.6f}",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(positions, [f"Rank {i}" for i in range(1, len(models) + 1)])
    ax.invert_yaxis()
    ax.set_xlim(0.5, 8.9)
    ax.set_xticks(range(1, 6), ["1st", "2nd", "3rd", "4th", "5th"])
    ax.set_xlabel("Balanced Accuracy rank (1st = highest)")
    ax.set_ylabel("Rank position")
    ax.set_title("Balanced Accuracy Model Ranking", loc="left", pad=14, fontweight="semibold")
    ax.grid(axis="x", color="#dce2e8", linewidth=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    save_figure(
        fig,
        "balanced_accuracy_rank_plot.png",
        ["balanced_accuracy_comparison.png"],
    )


def make_training_balanced_accuracy_scatter() -> None:
    efficient = pareto_efficient("Training Time (seconds)", "Balanced Accuracy")
    fig, ax = plt.subplots(figsize=(8.7, 5.7))
    for model in MODELS:
        x_value = RESULTS[model]["Training Time (seconds)"]
        y_value = RESULTS[model]["Balanced Accuracy"]
        ax.scatter(
            x_value,
            y_value,
            s=105,
            color=COLORS[model],
            edgecolor="white",
            linewidth=0.9,
            zorder=3,
        )
        offset = (8, 9) if model not in {"XGBoost", "LightGBM"} else (8, -16)
        ax.annotate(
            f"{model}\n({x_value:.3f}s, {y_value:.6f})",
            (x_value, y_value),
            xytext=offset,
            textcoords="offset points",
            fontsize=8.5,
        )
    frontier_x = [RESULTS[model]["Training Time (seconds)"] for model in efficient]
    frontier_y = [RESULTS[model]["Balanced Accuracy"] for model in efficient]
    ax.plot(
        frontier_x,
        frontier_y,
        color="#596a78",
        linestyle=(0, (4, 3)),
        linewidth=1.1,
        zorder=2,
        label="Pareto frontier",
    )
    ax.set_xlim(0, 62)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Training Time (seconds)")
    ax.set_ylabel("Balanced Accuracy")
    ax.set_title(
        "Performance vs Computation",
        loc="left",
        pad=14,
        fontweight="semibold",
    )
    ax.grid(color="#e0e5ea", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    ax.text(
        0.02,
        0.97,
        "Efficient frontier: " + " and ".join(efficient),
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        color="#3f4d59",
    )
    fig.tight_layout()
    save_figure(
        fig,
        "performance_vs_computation.png",
        ["training_time_comparison.png"],
    )


def make_model_fingerprint() -> None:
    angles = np.linspace(0, 2 * np.pi, len(METRICS), endpoint=False)
    closed_angles = np.append(angles, angles[0])
    fig, ax = plt.subplots(figsize=(9.3, 8.1), subplot_kw={"polar": True})
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_xticks(angles, METRICS)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], color="#66717b")
    ax.set_rlabel_position(15)
    for model in MODELS:
        values = [RESULTS[model][metric] for metric in METRICS]
        closed_values = values + [values[0]]
        ax.plot(
            closed_angles,
            closed_values,
            color=COLORS[model],
            linewidth=1.9,
            marker="o",
            markersize=4,
            label=model,
        )
        ax.fill(closed_angles, closed_values, color=COLORS[model], alpha=0.035)
    ax.grid(color="#cdd5dd", linewidth=0.8)
    ax.set_title("Model Performance Fingerprint", pad=28, fontweight="semibold")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=3,
        frameon=False,
    )
    fig.text(
        0.5,
        0.015,
        "All six metrics use the same 0–1 radial scale; larger values extend farther from the center.",
        ha="center",
        fontsize=9,
        color="#55616c",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.98))
    save_figure(
        fig,
        "model_performance_fingerprint.png",
        ["model_performance_summary.png"],
    )


def make_model_efficiency_scatter() -> None:
    efficient = pareto_efficient("Training Time (seconds)", "ROC-AUC")
    fig, ax = plt.subplots(figsize=(8.8, 5.7))
    for model in MODELS:
        x_value = RESULTS[model]["Training Time (seconds)"]
        y_value = RESULTS[model]["ROC-AUC"]
        ax.scatter(
            x_value,
            y_value,
            s=105,
            color=COLORS[model],
            edgecolor="white",
            linewidth=0.9,
            zorder=3,
        )
        offset = (8, 8)
        if model == "LightGBM":
            offset = (10, -19)
        elif model == "XGBoost":
            offset = (8, -18)
        ax.annotate(
            f"{model}\n({x_value:.3f}s, {y_value:.6f})",
            (x_value, y_value),
            xytext=offset,
            textcoords="offset points",
            fontsize=8.5,
        )
    frontier_x = [RESULTS[model]["Training Time (seconds)"] for model in efficient]
    frontier_y = [RESULTS[model]["ROC-AUC"] for model in efficient]
    ax.plot(
        frontier_x,
        frontier_y,
        color="#596a78",
        linestyle=(0, (4, 3)),
        linewidth=1.1,
        label="Pareto frontier",
        zorder=2,
    )
    ax.set_xlim(0, 62)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Training Time (seconds)")
    ax.set_ylabel("ROC-AUC")
    ax.set_title(
        "Model Efficiency: ROC-AUC vs Training Time",
        loc="left",
        pad=14,
        fontweight="semibold",
    )
    ax.grid(color="#e0e5ea", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    ax.annotate(
        "Fast, competitive ROC-AUC",
        xy=(
            RESULTS["LightGBM"]["Training Time (seconds)"],
            RESULTS["LightGBM"]["ROC-AUC"],
        ),
        xytext=(18, 24),
        textcoords="offset points",
        fontsize=8.5,
        color="#42515d",
        arrowprops={"arrowstyle": "->", "color": "#71808c", "lw": 0.9},
    )
    ax.text(
        0.02,
        0.97,
        "Pareto-efficient: " + " and ".join(efficient),
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        color="#3f4d59",
    )
    fig.tight_layout()
    save_figure(fig, "model_efficiency_roc_auc_vs_time.png")


def make_comparison_table(data: pd.DataFrame) -> None:
    metric_columns = [
        "ROC-AUC",
        "PR-AUC",
        "F1",
        "Precision",
        "Recall",
        "Balanced Accuracy",
        "Training Time (seconds)",
    ]
    winners = {
        metric: (
            min(MODELS, key=lambda model: data.loc[model, metric])
            if metric == "Training Time (seconds)"
            else max(MODELS, key=lambda model: data.loc[model, metric])
        )
        for metric in metric_columns
    }
    display_columns = [
        "Model",
        "ROC-AUC",
        "PR-AUC",
        "F1",
        "Precision",
        "Recall",
        "Balanced Accuracy",
        "Training Time",
    ]
    table_rows = []
    for model in MODELS:
        table_rows.append(
            [
                model,
                *[f"{data.loc[model, metric]:.6f}" for metric in metric_columns[:-1]],
                f"{data.loc[model, 'Training Time (seconds)']:.3f}",
            ]
        )

    fig, ax = plt.subplots(figsize=(15.5, 4.5))
    ax.axis("off")
    table = ax.table(
        cellText=table_rows,
        colLabels=display_columns,
        cellLoc="center",
        colLoc="center",
        loc="center",
        colWidths=[0.19, 0.105, 0.105, 0.09, 0.10, 0.09, 0.15, 0.13],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.8)
    for (row, column), cell in table.get_celld().items():
        cell.set_edgecolor("#d5dce2")
        if row == 0:
            cell.set_facecolor("#263b4d")
            cell.set_text_props(color="white", weight="bold")
        else:
            model = MODELS[row - 1]
            cell.set_facecolor("#f3f6f8" if row % 2 == 0 else "white")
            if column == 0:
                cell.set_text_props(ha="left", weight="medium")
            if column > 0 and winners[metric_columns[column - 1]] == model:
                cell.set_facecolor("#dcefe8")
                cell.set_text_props(weight="bold", color="#174d3b")
    ax.set_title(
        "Baseline Model Evaluation — Test Set",
        loc="left",
        pad=20,
        fontsize=14,
        fontweight="semibold",
    )
    fig.text(
        0.02,
        0.04,
        "Bold green cells identify the best value in each metric column; lowest training time is best.",
        fontsize=9,
        color="#55616c",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save_figure(fig, "model_performance_comparison_table.png")


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    data = pd.DataFrame.from_dict(RESULTS, orient="index")
    data.index.name = "Model"
    data.to_csv(DATA_DIR / "model_performance_comparison.csv", index=True, float_format="%.6f")

    make_roc_radar()
    make_pr_auc_dotplot()
    make_f1_dumbbell()
    make_precision_beeswarm()
    make_recall_lollipop()
    make_balanced_accuracy_ranking()
    make_training_balanced_accuracy_scatter()
    make_model_fingerprint()
    make_model_efficiency_scatter()
    make_comparison_table(data)

    print("Generated ten 300-DPI figures and model_performance_comparison.csv.")
    print("ROC-AUC Pareto-efficient:", ", ".join(pareto_efficient("Training Time (seconds)", "ROC-AUC")))
    print(
        "Balanced-accuracy Pareto-efficient:",
        ", ".join(pareto_efficient("Training Time (seconds)", "Balanced Accuracy")),
    )


if __name__ == "__main__":
    main()
