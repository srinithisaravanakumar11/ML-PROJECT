"""Central project paths so every script resolves files the same way,
regardless of the directory it is launched from."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"
DATA_DIR = OUTPUT_DIR / "data"
REPORTS_DIR = OUTPUT_DIR / "reports"
FIGURES_DIR = OUTPUT_DIR / "figures"
MODEL_DIR = OUTPUT_DIR / "saved_models"

# Raw ClinVar download (gitignored); see README "Reproducing the pipeline".
RAW_CLINVAR = ROOT / "variant_summary.txt.gz"

# Large intermediate datasets (gitignored).
ML_DATASET = DATA_DIR / "clinvar_ml_dataset.csv"
ENHANCED_DATASET = DATA_DIR / "clinvar_enhanced_dataset.csv"

for _d in (DATA_DIR, REPORTS_DIR, FIGURES_DIR):
    _d.mkdir(parents=True, exist_ok=True)
