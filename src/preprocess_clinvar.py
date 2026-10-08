#!/usr/bin/env python3
"""
ClinVar Variant Summary Preprocessing Pipeline
================================================
Purpose: Clean and preprocess variant_summary.txt.gz for binary
         pathogenicity classification (Pathogenic vs Benign).

Author: Genome Variant Research Project
Date:   2026-10-05

This script:
  1.  Inspects the raw dataset (schema, types, missing values)
  2.  Maps clinical significance → binary target
  3.  Removes unreliable records (conservative filtering)
  4.  Resolves duplicate / conflicting variants
  5.  Cleans categorical fields
  6.  Handles missing values
  7.  Selects features, flags leakage columns
  8.  Produces Dataset A (ClinVar baseline) & Dataset B (extended stub)
  9.  Recommends splitting / imbalance strategies
  10. Generates all output artefacts

Outputs (written to ./output/):
  - clinvar_cleaned.csv
  - clinvar_ml_dataset.csv
  - clinvar_preprocessing_report.txt
  - removed_records.csv
  - feature_dictionary.csv
  - (this script itself: preprocess_clinvar.py)
"""

import argparse
import gzip
import os
import sys
import warnings
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # Windows consoles default to cp1252

warnings.filterwarnings("ignore")

# ──────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent.parent  # project root
DEFAULT_INPUT_FILE = "variant_summary.txt.gz"
DEFAULT_OUTPUT_DIR = "output/data"
INPUT_FILE = str((SCRIPT_DIR / DEFAULT_INPUT_FILE).resolve())
OUTPUT_DIR = str((SCRIPT_DIR / DEFAULT_OUTPUT_DIR).resolve())
REPORT_LINES = []          # accumulates the preprocessing report


def parse_args():
    """Parse command-line arguments with sensible defaults relative to the project root."""
    parser = argparse.ArgumentParser(description="Preprocess ClinVar variant_summary.txt.gz for binary pathogenicity classification.")
    parser.add_argument("--input", default=DEFAULT_INPUT_FILE, help="Path to the raw ClinVar TSV/GZ file.")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Directory for all generated outputs.")
    return parser.parse_args()


def resolve_input_output_paths(input_path: str, output_dir: str):
    """Resolve input and output paths relative to the project root."""
    input_path_obj = Path(input_path)
    if not input_path_obj.is_absolute():
        input_path_obj = (SCRIPT_DIR / input_path_obj).resolve()

    output_dir_obj = Path(output_dir)
    if not output_dir_obj.is_absolute():
        output_dir_obj = (SCRIPT_DIR / output_dir_obj).resolve()

    return str(input_path_obj), str(output_dir_obj)


def report(msg: str, also_print: bool = True):
    """Append a line to the report and optionally print it."""
    REPORT_LINES.append(msg)
    if also_print:
        print(msg)


def section(title: str):
    """Print and record a section header."""
    bar = "=" * 80
    report("")
    report(bar)
    report(f"  {title}")
    report(bar)


def filter_log(before: int, removed: int, reason: str):
    """Standard filtering log entry."""
    after = before - removed
    report(f"  Rows before: {before:>10,}")
    report(f"  Rows removed:{removed:>10,}  — {reason}")
    report(f"  Rows after:  {after:>10,}")
    report("")
    return after


# ──────────────────────────────────────────────
# 0. Setup
# ──────────────────────────────────────────────
def main():
    global INPUT_FILE, OUTPUT_DIR, REPORT_LINES
    args = parse_args()
    INPUT_FILE, OUTPUT_DIR = resolve_input_output_paths(args.input, args.output_dir)
    REPORT_LINES = []

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    start_time = datetime.now()
    section("0. SETUP")
    report(f"  Script started at: {start_time.isoformat()}")
    report(f"  Input file: {INPUT_FILE}")
    report(f"  Output directory: {OUTPUT_DIR}/")

    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    # ──────────────────────────────────────────────
    # 1. INSPECT THE DATASET
    # ──────────────────────────────────────────────
    section("1. DATASET INSPECTION")

    report("  Loading variant_summary.txt.gz …")
    df_raw = pd.read_csv(INPUT_FILE, sep="\t", compression="gzip", low_memory=False)
    original_row_count = len(df_raw)
    original_col_count = len(df_raw.columns)

    report(f"  Rows:    {original_row_count:,}")
    report(f"  Columns: {original_col_count}")

    # Column names
    report("\n  --- Column Names ---")
    for i, col in enumerate(df_raw.columns):
        report(f"    [{i:>2}] {col}")

    # Data types
    report("\n  --- Data Types ---")
    for col in df_raw.columns:
        report(f"    {col:<40s}  {str(df_raw[col].dtype)}")

    # Missing values
    report("\n  --- Missing Values ---")
    miss = df_raw.isnull().sum()
    miss_pct = (miss / len(df_raw) * 100).round(2)
    for col in df_raw.columns:
        report(f"    {col:<40s}  {miss[col]:>10,}  ({miss_pct[col]:>6.2f}%)")

    # Unique-value counts for likely-categorical columns
    report("\n  --- Unique Values (categorical candidates) ---")
    cat_candidates = df_raw.select_dtypes(include=["object"]).columns.tolist()
    for col in cat_candidates:
        n = df_raw[col].nunique()
        report(f"    {col:<40s}  {n:>10,} unique")

    # Sample rows
    report("\n  --- Sample Rows (first 5) ---")
    sample = df_raw.head(5).to_string(max_colwidth=40)
    for line in sample.split("\n"):
        report(f"    {line}")

    # ──────────────────────────────────────────────
    # 1b. Map known ClinVar columns to semantic roles
    # ──────────────────────────────────────────────
    report("\n  --- Column Role Mapping ---")
    # Build a mapping by inspecting actual column names
    col_lower_map = {c.lower().replace(" ", "").replace("#", ""): c for c in df_raw.columns}

    # Helper: find a column whose lower-cased name contains ANY of the given substrings
    def find_col(*substrings):
        """Prefer exact column names, then fall back to substring matches."""
        for sub in substrings:
            for raw_name in df_raw.columns:
                if sub.lower() == raw_name.lower():
                    return raw_name
        for sub in substrings:
            for raw_name in df_raw.columns:
                if sub.lower() in raw_name.lower():
                    return raw_name
        return None

    COL_MAP = OrderedDict()
    COL_MAP["variant_id"]       = find_col("VariationID", "AlleleID", "#AlleleID")
    COL_MAP["chromosome"]       = find_col("Chromosome")
    COL_MAP["position_start"]   = find_col("PositionVCF", "Start")
    COL_MAP["position_stop"]    = find_col("Stop")
    COL_MAP["ref_allele"]       = find_col("ReferenceAlleleVCF", "ReferenceAllele", "RefAllele")
    COL_MAP["alt_allele"]       = find_col("AlternateAlleleVCF", "AlternateAllele", "AltAllele")
    COL_MAP["gene_symbol"]      = find_col("GeneSymbol", "Gene")
    COL_MAP["gene_id"]          = find_col("GeneID")
    COL_MAP["variant_type"]     = find_col("Type")
    COL_MAP["consequence"]      = find_col("MolecularConsequence", "Consequence")
    COL_MAP["clinical_sig"]     = find_col("ClinicalSignificance", "ClinSig")
    COL_MAP["review_status"]    = find_col("ReviewStatus")
    COL_MAP["phenotype"]        = find_col("PhenotypeIDS", "Phenotype")
    COL_MAP["phenotype_list"]   = find_col("PhenotypeList")
    COL_MAP["assembly"]         = find_col("Assembly")
    COL_MAP["rs_id"]            = find_col("RS#", "rsid", "RS (dbSNP)")
    COL_MAP["nsv_esv"]          = find_col("nsv/esv", "nsv")
    COL_MAP["origin"]           = find_col("Origin")
    COL_MAP["rcv_accession"]    = find_col("RCVaccession", "RCV")
    COL_MAP["variation_id"]     = find_col("VariationID")
    COL_MAP["name"]             = find_col("Name")
    COL_MAP["submitted_assembly"] = find_col("SubmittedAssembly")
    COL_MAP["number_submitters"]  = find_col("NumberOfSubmitters")
    COL_MAP["guidelines"]       = find_col("Guidelines")
    COL_MAP["other_ids"]        = find_col("OtherIDs")

    for role, actual in COL_MAP.items():
        status = actual if actual else "** NOT FOUND **"
        report(f"    {role:<25s} → {status}")

    # Convenience accessors (will be None if column not found)
    COL_CHROM     = COL_MAP["chromosome"]
    COL_START     = COL_MAP["position_start"]
    COL_STOP      = COL_MAP["position_stop"]
    COL_REF       = COL_MAP["ref_allele"]
    COL_ALT       = COL_MAP["alt_allele"]
    COL_GENE      = COL_MAP["gene_symbol"]
    COL_GENEID    = COL_MAP["gene_id"]
    COL_TYPE      = COL_MAP["variant_type"]
    COL_CONSEQ    = COL_MAP["consequence"]
    COL_CLINSIG   = COL_MAP["clinical_sig"]
    COL_REVIEW    = COL_MAP["review_status"]
    COL_ASSEMBLY  = COL_MAP["assembly"]
    COL_RSID      = COL_MAP["rs_id"]
    COL_ALLELE_ID = COL_MAP["variant_id"]
    COL_VAR_ID    = COL_MAP["variation_id"]
    COL_ORIGIN    = COL_MAP["origin"]
    COL_PHENO     = COL_MAP["phenotype_list"]
    COL_NAME      = COL_MAP["name"]
    COL_NSUBS     = COL_MAP["number_submitters"]

    # ──────────────────────────────────────────────
    # 2. DETERMINE TARGET LABEL
    # ──────────────────────────────────────────────
    section("2. TARGET LABEL DEFINITION")

    if COL_CLINSIG is None:
        report("  ERROR: Cannot find ClinicalSignificance column. Aborting.")
        sys.exit(1)

    # Show raw distribution
    report(f"  Column used: {COL_CLINSIG}")
    sig_counts = df_raw[COL_CLINSIG].value_counts(dropna=False)
    report("\n  --- Raw ClinicalSignificance Distribution (top 30) ---")
    for val, cnt in sig_counts.head(30).items():
        pct = cnt / len(df_raw) * 100
        report(f"    {str(val):<60s}  {cnt:>10,}  ({pct:>5.2f}%)")

    # Define mapping
    PATHOGENIC_TERMS = [
        "pathogenic",
        "pathogenic/likely pathogenic",
        "likely pathogenic",
    ]
    BENIGN_TERMS = [
        "benign",
        "benign/likely benign",
        "likely benign",
    ]
    UNCERTAIN_TERMS = [
        "uncertain significance",
        "conflicting interpretations of pathogenicity",
        "conflicting classifications of pathogenicity",
        "not provided",
        "drug response",
        "risk factor",
        "association",
        "protective",
        "affects",
        "other",
        "confers sensitivity",
        "uncertain risk allele",
        "likely risk allele",
        "established risk allele",
    ]

    def map_clinsig(value):
        """Map a ClinicalSignificance string to Pathogenic / Benign / Excluded."""
        if pd.isna(value):
            return "Excluded"
        val = str(value).strip().lower()
        if val in PATHOGENIC_TERMS:
            return "Pathogenic"
        if val in BENIGN_TERMS:
            return "Benign"
        return "Excluded"

    def map_clinsig_strict(value):
        """Strict mapping: only pure Pathogenic / Benign (no Likely)."""
        if pd.isna(value):
            return "Excluded"
        val = str(value).strip().lower()
        if val == "pathogenic":
            return "Pathogenic"
        if val == "benign":
            return "Benign"
        return "Excluded"

    # We use the broader mapping that includes Likely pathogenic / Likely benign
    # This is standard in ClinVar-based ML studies and explicitly documented.
    df_raw["target"] = df_raw[COL_CLINSIG].apply(map_clinsig)
    df_raw["target_strict"] = df_raw[COL_CLINSIG].apply(map_clinsig_strict)
    df_raw["label_definition"] = df_raw[COL_CLINSIG]  # preserve original mapping

    report("\n  --- Target Distribution (broad: includes Likely) ---")
    td = df_raw["target"].value_counts()
    for val, cnt in td.items():
        report(f"    {val:<20s}  {cnt:>10,}  ({cnt/len(df_raw)*100:.2f}%)")

    report("\n  --- Target Distribution (strict: only pure P/B) ---")
    td2 = df_raw["target_strict"].value_counts()
    for val, cnt in td2.items():
        report(f"    {val:<20s}  {cnt:>10,}  ({cnt/len(df_raw)*100:.2f}%)")

    report("\n  --- Excluded Category Breakdown ---")
    excluded_mask = df_raw["target"] == "Excluded"
    excluded_dist = df_raw.loc[excluded_mask, COL_CLINSIG].value_counts(dropna=False)
    for val, cnt in excluded_dist.head(20).items():
        report(f"    {str(val):<60s}  {cnt:>10,}")

    report("""
      DECISION: We use the BROAD mapping for binary classification:
        Pathogenic class = 'Pathogenic' + 'Pathogenic/Likely pathogenic' + 'Likely pathogenic'
        Benign class     = 'Benign' + 'Benign/Likely benign' + 'Likely benign'

      RATIONALE:
        - 'Likely pathogenic' and 'Likely benign' carry strong clinical evidence.
        - Combining them is standard practice in published ClinVar ML studies.
        - The 'label_definition' column preserves the original ClinVar label
          so users can re-stratify later.

      EXCLUDED categories (NOT used in binary classification):
        - Uncertain significance — ambiguous; would add noise
        - Conflicting interpretations — no clear ground truth
        - Not provided, Drug response, Risk factor, etc. — not binary P/B
      These are saved in removed_records.csv for transparency.
    """)

    # ──────────────────────────────────────────────
    # 3. REMOVE UNRELIABLE RECORDS
    # ──────────────────────────────────────────────
    section("3. FILTERING UNRELIABLE RECORDS")

    # Start with only Pathogenic / Benign rows
    removed_records = []  # list of (index, reason) pairs

    # 3a. Keep only target = Pathogenic or Benign
    before = len(df_raw)
    mask_target = df_raw["target"].isin(["Pathogenic", "Benign"])
    excluded_target = df_raw[~mask_target].copy()
    excluded_target["exclusion_reason"] = "Non-binary clinical significance"
    removed_records.append(excluded_target)

    df = df_raw[mask_target].copy()
    filter_log(before, before - len(df), "Non-binary clinical significance (Excluded category)")

    # 3b. Filter to a single genome assembly (prefer GRCh38, fall back to GRCh37)
    if COL_ASSEMBLY:
        report("  --- Assembly Distribution ---")
        asm_counts = df[COL_ASSEMBLY].value_counts(dropna=False)
        for val, cnt in asm_counts.items():
            report(f"    {str(val):<30s}  {cnt:>10,}")

        preferred_assembly = "GRCh38"
        if preferred_assembly not in df[COL_ASSEMBLY].values:
            preferred_assembly = "GRCh37"
        report(f"\n  Selected assembly: {preferred_assembly}")

        before = len(df)
        mask_asm = df[COL_ASSEMBLY] == preferred_assembly
        exc = df[~mask_asm].copy()
        exc["exclusion_reason"] = f"Assembly != {preferred_assembly}"
        removed_records.append(exc)
        df = df[mask_asm].copy()
        filter_log(before, before - len(df), f"Assembly != {preferred_assembly}")

    # 3c. Missing target labels (should be zero after step 3a, but safety check)
    before = len(df)
    mask_miss_target = df["target"].isna()
    exc = df[mask_miss_target].copy()
    exc["exclusion_reason"] = "Missing target label"
    removed_records.append(exc)
    df = df[~mask_miss_target].copy()
    filter_log(before, before - len(df), "Missing target label")

    # 3d. Missing essential variant information (Chromosome, Start position)
    essential_cols = [c for c in [COL_CHROM, COL_START] if c is not None]
    if essential_cols:
        before = len(df)
        mask_essential = df[essential_cols].isna().any(axis=1)
        exc = df[mask_essential].copy()
        exc["exclusion_reason"] = "Missing chromosome or position"
        removed_records.append(exc)
        df = df[~mask_essential].copy()
        filter_log(before, before - len(df), "Missing chromosome or position")

    # 3e. Invalid chromosome values
    if COL_CHROM:
        valid_chroms = set([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "x", "y", "mt"])
        before = len(df)
        df[COL_CHROM] = df[COL_CHROM].astype(str).str.strip()
        mask_valid_chrom = df[COL_CHROM].isin(valid_chroms)
        exc = df[~mask_valid_chrom].copy()
        exc["exclusion_reason"] = "Invalid chromosome value"
        removed_records.append(exc)
        df = df[mask_valid_chrom].copy()
        filter_log(before, before - len(df), "Invalid chromosome value")

    # 3f. Invalid position (non-numeric or negative)
    if COL_START:
        before = len(df)
        df[COL_START] = pd.to_numeric(df[COL_START], errors="coerce")
        mask_pos = df[COL_START].isna() | (df[COL_START] < 0)
        exc = df[mask_pos].copy()
        exc["exclusion_reason"] = "Invalid position value"
        removed_records.append(exc)
        df = df[~mask_pos].copy()
        df[COL_START] = df[COL_START].astype(int)
        filter_log(before, before - len(df), "Invalid position value")

    # 3g. Missing or invalid reference / alternate alleles
    # ClinVar uses 'na', '-', '' for missing alleles — keep '-' as it can mean deletion
    if COL_REF and COL_ALT:
        before = len(df)
        df[COL_REF] = df[COL_REF].astype(str).str.strip()
        df[COL_ALT] = df[COL_ALT].astype(str).str.strip()
        invalid_allele_vals = {"nan", "", "na"}
        mask_allele = (
            df[COL_REF].str.lower().isin(invalid_allele_vals) |
            df[COL_ALT].str.lower().isin(invalid_allele_vals)
        )
        exc = df[mask_allele].copy()
        exc["exclusion_reason"] = "Missing or invalid ref/alt allele"
        removed_records.append(exc)
        df = df[~mask_allele].copy()
        filter_log(before, before - len(df), "Missing or invalid ref/alt allele")

    # 3h. Records with very low review evidence
    # ClinVar review statuses (ascending confidence):
    #   no assertion criteria provided
    #   criteria provided, single submitter
    #   criteria provided, conflicting interpretations
    #   criteria provided, multiple submitters, no conflicts
    #   reviewed by expert panel
    #   practice guideline
    if COL_REVIEW:
        report("  --- Review Status Distribution (after prior filters) ---")
        rs_counts = df[COL_REVIEW].value_counts(dropna=False)
        for val, cnt in rs_counts.items():
            report(f"    {str(val):<65s}  {cnt:>10,}")

        # Remove only "no assertion" records — very low reliability
        before = len(df)
        low_evidence_patterns = [
            "no assertion",
            "no classification",
        ]
        mask_low = df[COL_REVIEW].astype(str).str.lower().apply(
            lambda x: any(pat in x for pat in low_evidence_patterns)
        )
        exc = df[mask_low].copy()
        exc["exclusion_reason"] = "Insufficient review evidence (no assertion criteria)"
        removed_records.append(exc)
        df = df[~mask_low].copy()
        filter_log(before, before - len(df), "Insufficient review evidence")

    # ──────────────────────────────────────────────
    # 4. HANDLE DUPLICATE VARIANTS
    # ──────────────────────────────────────────────
    section("4. DUPLICATE VARIANT HANDLING")

    # Define variant key
    VARIANT_KEY_COLS = []
    if COL_CHROM:
        VARIANT_KEY_COLS.append(COL_CHROM)
    if COL_START:
        VARIANT_KEY_COLS.append(COL_START)
    if COL_REF:
        VARIANT_KEY_COLS.append(COL_REF)
    if COL_ALT:
        VARIANT_KEY_COLS.append(COL_ALT)

    report(f"  Variant key columns: {VARIANT_KEY_COLS}")

    if len(VARIANT_KEY_COLS) >= 2:
        df["_variant_key"] = df[VARIANT_KEY_COLS].astype(str).agg("|".join, axis=1)
        dup_counts = df["_variant_key"].value_counts()
        n_dup_keys = (dup_counts > 1).sum()
        n_dup_rows = dup_counts[dup_counts > 1].sum()
        report(f"  Unique variant keys: {dup_counts.shape[0]:,}")
        report(f"  Keys with >1 row:    {n_dup_keys:,}")
        report(f"  Total duplicate rows: {n_dup_rows:,}")

        # Analyse conflict among duplicates
        dup_keys = dup_counts[dup_counts > 1].index
        dup_df = df[df["_variant_key"].isin(dup_keys)]

        # For each duplicate key, check if all targets agree
        dup_target = dup_df.groupby("_variant_key")["target"].nunique()
        n_conflict = (dup_target > 1).sum()
        n_agree = (dup_target == 1).sum()
        report(f"  Duplicate keys with SAME target:      {n_agree:,}")
        report(f"  Duplicate keys with CONFLICTING target: {n_conflict:,}")

        # Strategy:
        # - Conflicting targets → REMOVE (unreliable ground truth)
        # - Agreeing targets → keep the row with the highest review status
        report("""
      DUPLICATE STRATEGY:
        1. Conflicting targets (same variant, different P/B label) → Remove entirely.
           Reason: no reliable ground truth for that variant.
        2. Agreeing targets → Keep the row with the highest review confidence.
           Review status priority (high → low):
             practice guideline > reviewed by expert panel >
             criteria provided, multiple submitters, no conflicts >
             criteria provided, conflicting interpretations >
             criteria provided, single submitter
        3. Build a unique variant key for downstream train/test splitting
           to prevent data leakage.
    """)

        # Remove conflicting-target duplicates
        conflicting_keys = dup_target[dup_target > 1].index
        before = len(df)
        mask_conflict = df["_variant_key"].isin(conflicting_keys)
        exc = df[mask_conflict].copy()
        exc["exclusion_reason"] = "Conflicting clinical significance for same variant"
        removed_records.append(exc)
        df = df[~mask_conflict].copy()
        filter_log(before, before - len(df), "Conflicting clinical significance duplicates")

        # For agreeing duplicates, keep highest-review-status row
        REVIEW_PRIORITY = {
            "practice guideline": 6,
            "reviewed by expert panel": 5,
            "criteria provided, multiple submitters, no conflicts": 4,
            "criteria provided, conflicting interpretations": 3,
            "criteria provided, conflicting classifications": 3,
            "criteria provided, single submitter": 2,
        }

        if COL_REVIEW:
            df["_review_rank"] = df[COL_REVIEW].astype(str).str.lower().str.strip().map(REVIEW_PRIORITY).fillna(0)
        else:
            df["_review_rank"] = 0

        before = len(df)
        # Sort by review rank descending, keep first per variant key
        df = df.sort_values("_review_rank", ascending=False)
        dup_mask_keep = df.duplicated(subset=["_variant_key"], keep="first")
        exc = df[dup_mask_keep].copy()
        exc["exclusion_reason"] = "Duplicate variant (lower review evidence)"
        removed_records.append(exc)
        df = df[~dup_mask_keep].copy()
        filter_log(before, before - len(df), "Duplicate variants (kept highest review)")

        df.drop(columns=["_review_rank"], inplace=True)
    else:
        report("  WARNING: Insufficient columns to build variant key. Skipping dedup.")
        df["_variant_key"] = range(len(df))

    # ──────────────────────────────────────────────
    # 5. CLEAN CATEGORICAL FIELDS
    # ──────────────────────────────────────────────
    section("5. CATEGORICAL FIELD CLEANING")

    # 5a. Chromosome — standardize
    if COL_CHROM:
        df[COL_CHROM] = df[COL_CHROM].astype(str).str.strip().str.upper()
        # Standardize MT
        df[COL_CHROM] = df[COL_CHROM].replace({"MT": "MT", "M": "MT"})
        report(f"  Chromosome unique values: {sorted(df[COL_CHROM].unique())}")

    # 5b. Clinical significance — already mapped; clean original column too
    if COL_CLINSIG:
        df[COL_CLINSIG] = df[COL_CLINSIG].astype(str).str.strip()

    # 5c. Variant type
    if COL_TYPE:
        df[COL_TYPE] = df[COL_TYPE].astype(str).str.strip().str.lower()
        report(f"  Variant type unique values ({df[COL_TYPE].nunique()}):")
        for val, cnt in df[COL_TYPE].value_counts().head(15).items():
            report(f"    {val:<40s}  {cnt:>10,}")

    # 5d. Molecular consequence
    if COL_CONSEQ:
        df[COL_CONSEQ] = df[COL_CONSEQ].astype(str).str.strip()
        report(f"\n  Molecular consequence unique values ({df[COL_CONSEQ].nunique()}):")
        for val, cnt in df[COL_CONSEQ].value_counts().head(15).items():
            report(f"    {val:<50s}  {cnt:>10,}")

    # 5e. Review status
    if COL_REVIEW:
        df[COL_REVIEW] = df[COL_REVIEW].astype(str).str.strip().str.lower()
        report(f"\n  Review status unique values ({df[COL_REVIEW].nunique()}):")
        for val, cnt in df[COL_REVIEW].value_counts().items():
            report(f"    {val:<65s}  {cnt:>10,}")

    # 5f. Gene symbol
    if COL_GENE:
        df[COL_GENE] = df[COL_GENE].astype(str).str.strip()

    # 5g. Origin
    if COL_ORIGIN:
        df[COL_ORIGIN] = df[COL_ORIGIN].astype(str).str.strip().str.lower()

    # ──────────────────────────────────────────────
    # 6. HANDLE MISSING VALUES
    # ──────────────────────────────────────────────
    section("6. MISSING VALUE ANALYSIS")

    important_features = [c for c in [
        COL_CHROM, COL_START, COL_STOP, COL_REF, COL_ALT,
        COL_GENE, COL_TYPE, COL_CONSEQ, COL_REVIEW, COL_ORIGIN,
        COL_RSID, COL_GENEID, COL_NAME, COL_NSUBS
    ] if c is not None]

    report("  --- Missing Value Summary for Important Features ---")
    report(f"  {'Feature':<40s}  {'Missing':>10s}  {'Pct':>8s}  {'Decision'}")
    report(f"  {'-'*40}  {'-'*10}  {'-'*8}  {'-'*35}")

    missing_decisions = {}
    for col in important_features:
        # Treat 'nan', 'na', '-1' strings as missing for reporting
        is_null = df[col].isna()
        is_str_null = df[col].astype(str).str.lower().isin(["nan", ""])
        total_miss = (is_null | is_str_null).sum()
        pct = total_miss / len(df) * 100

        if col in [COL_CHROM, COL_START]:
            decision = "REQUIRED — rows already removed"
        elif pct > 80:
            decision = "DROP COLUMN — >80% missing"
        elif pct > 50:
            decision = "KEEP — mark as 'unknown'"
        elif col == COL_CONSEQ:
            decision = "KEEP — fill missing with 'unknown'"
        elif col == COL_GENE:
            decision = "KEEP — fill missing with 'unknown'"
        elif col == COL_RSID:
            decision = "KEEP — optional identifier"
        else:
            decision = "KEEP — acceptable missingness"

        missing_decisions[col] = decision
        report(f"  {col:<40s}  {total_miss:>10,}  {pct:>7.2f}%  {decision}")

    # Apply decisions: fill selected columns
    for col, dec in missing_decisions.items():
        if "fill missing with" in dec:
            fill_val = "unknown"
            mask = df[col].isna() | (df[col].astype(str).str.lower().isin(["nan", ""]))
            df.loc[mask, col] = fill_val
            report(f"\n  Filled {mask.sum():,} missing values in '{col}' with '{fill_val}'")
        if "DROP COLUMN" in dec:
            report(f"\n  Dropping column '{col}' (>80% missing)")
            df.drop(columns=[col], inplace=True, errors="ignore")

    report("""
      MISSING VALUE POLICY:
        - Essential fields (chromosome, position): rows already removed if missing
        - Biological values are NOT filled with 0 or arbitrary numbers
        - 'unknown' is used to distinguish genuinely missing values from 'not applicable'
        - No imputation of biological features (alleles, sequence, scores)
    """)

    # ──────────────────────────────────────────────
    # 7. FEATURE SELECTION
    # ──────────────────────────────────────────────
    section("7. FEATURE SELECTION")

    # Identify potential leakage columns
    leakage_columns = []
    leakage_candidates = {
        COL_CLINSIG: "Directly encodes the target",
        COL_REVIEW: "Correlated with interpretation confidence, but not the target itself. SAFE with caution.",
        COL_PHENO: "Disease associations may correlate with pathogenicity",
    }
    report("  --- Potential Leakage Columns ---")
    for col, reason in leakage_candidates.items():
        if col and col in df.columns:
            report(f"    {col:<40s}  {reason}")
            if "Directly" in reason:
                leakage_columns.append(col)

    # Select features for ML dataset
    ml_features = []
    feature_descriptions = []

    def add_feature(col, desc, role):
        if col and col in df.columns:
            ml_features.append(col)
            feature_descriptions.append({"feature": col, "description": desc, "role": role})

    add_feature(COL_CHROM, "Chromosome (1-22, X, Y, MT)", "categorical")
    add_feature(COL_START, "Genomic start position (GRCh38)", "numeric")
    if COL_STOP and COL_STOP in df.columns:
        add_feature(COL_STOP, "Genomic stop position (GRCh38)", "numeric")
    add_feature(COL_REF, "Reference allele", "categorical/sequence")
    add_feature(COL_ALT, "Alternate allele", "categorical/sequence")
    add_feature(COL_GENE, "Gene symbol (HGNC)", "categorical")
    add_feature(COL_TYPE, "Variant type (e.g., single nucleotide variant)", "categorical")
    add_feature(COL_CONSEQ, "Molecular consequence (SO terms)", "categorical")
    add_feature(COL_REVIEW, "ClinVar review status (evidence level)", "ordinal")
    add_feature(COL_ORIGIN, "Allele origin (germline, somatic, etc.)", "categorical")

    # Derived features
    # Variant length
    if COL_START and COL_STOP and COL_STOP in df.columns:
        df["variant_length"] = (pd.to_numeric(df[COL_STOP], errors="coerce") -
                                pd.to_numeric(df[COL_START], errors="coerce")).fillna(0).astype(int)
        ml_features.append("variant_length")
        feature_descriptions.append({"feature": "variant_length",
                                     "description": "Length of variant (stop - start)", "role": "numeric"})

    # Ref/Alt allele length
    if COL_REF and COL_REF in df.columns:
        df["ref_length"] = df[COL_REF].astype(str).apply(lambda x: len(x) if x != "-" else 0)
        ml_features.append("ref_length")
        feature_descriptions.append({"feature": "ref_length",
                                     "description": "Length of reference allele", "role": "numeric"})

    if COL_ALT and COL_ALT in df.columns:
        df["alt_length"] = df[COL_ALT].astype(str).apply(lambda x: len(x) if x != "-" else 0)
        ml_features.append("alt_length")
        feature_descriptions.append({"feature": "alt_length",
                                     "description": "Length of alternate allele", "role": "numeric"})

    # Is indel flag
    if COL_REF and COL_ALT:
        df["is_indel"] = (df["ref_length"] != df["alt_length"]).astype(int)
        ml_features.append("is_indel")
        feature_descriptions.append({"feature": "is_indel",
                                     "description": "1 if insertion/deletion, 0 if substitution", "role": "binary"})

    report("\n  --- Selected ML Features ---")
    for fd in feature_descriptions:
        report(f"    {fd['feature']:<30s}  [{fd['role']:<20s}]  {fd['description']}")

    report(f"\n  Total features: {len(ml_features)}")
    report(f"  Target column: 'target'")
    report(f"  Leakage columns excluded: {leakage_columns}")

    # ──────────────────────────────────────────────
    # 8. EXTERNAL FEATURES (not fabricated)
    # ──────────────────────────────────────────────
    section("8. ADDITIONAL FEATURES (not in ClinVar — DO NOT FABRICATE)")

    report("""
      The following features are NOT present in variant_summary.txt.gz.
      DO NOT create fake values. They should be obtained from external databases.

      Feature                   Source                  Join Key
      ─────────────────────────────────────────────────────────────────────
      Population frequency      gnomAD                  chr:pos:ref:alt
      Conservation score        UCSC / phastCons        chr:pos (liftover if needed)
      CADD score                CADD database           chr:pos:ref:alt
      SIFT score                Ensembl VEP             chr:pos:ref:alt or rs_id
      PolyPhen2 score           Ensembl VEP             chr:pos:ref:alt or rs_id
      DNA sequence context      UCSC Genome Browser     chr:pos (±50bp window)
      Protein domain            InterPro / UniProt      Gene + protein position
      GC content                Reference genome        chr:pos (local window)
      Exon / intron location    Ensembl / GENCODE       chr:pos + gene

      How to join:
        1. Build a variant key: chromosome + ":" + position + ":" + ref + ":" + alt
        2. Use tabix-indexed VCFs from gnomAD for efficient lookup
        3. Use Ensembl REST API or VEP for consequence annotations
        4. For sequence context, query the reference genome FASTA (hg38)
    """)

    # ──────────────────────────────────────────────
    # 9. CREATE DATASETS A & B
    # ──────────────────────────────────────────────
    section("9. DATASET CREATION")

    # Variant ID for the dataset
    if "_variant_key" in df.columns:
        df["variant_id"] = df["_variant_key"]
    else:
        df["variant_id"] = range(len(df))

    # Dataset A: ClinVar baseline
    dataset_a_cols = ["variant_id"] + ml_features + ["target", "label_definition"]
    dataset_a_cols = [c for c in dataset_a_cols if c in df.columns]
    dataset_a = df[dataset_a_cols].copy()

    report(f"  Dataset A (ClinVar baseline):")
    report(f"    Rows:     {len(dataset_a):,}")
    report(f"    Columns:  {len(dataset_a.columns)}")
    report(f"    Features: {ml_features}")

    # Dataset B: Extended stub (same as A, plus placeholder columns)
    dataset_b = dataset_a.copy()
    # Add empty columns for future external features — but do NOT fill them
    external_feature_stubs = [
        "gnomad_af", "cadd_score", "sift_score", "polyphen2_score",
        "phastcons_score", "gc_content", "seq_context_50bp"
    ]
    for stub in external_feature_stubs:
        dataset_b[stub] = np.nan

    report(f"\n  Dataset B (Extended research stub):")
    report(f"    Rows:     {len(dataset_b):,}")
    report(f"    Columns:  {len(dataset_b.columns)}")
    report(f"    Stub columns (to fill from external sources): {external_feature_stubs}")

    # ──────────────────────────────────────────────
    # 10. DATA LEAKAGE PREVENTION
    # ──────────────────────────────────────────────
    section("10. DATA LEAKAGE PREVENTION")

    report("""
      CRITICAL: Why random row-level splitting is DANGEROUS for genomic data
      ─────────────────────────────────────────────────────────────────────
      Even after deduplication, related variants can leak information:

      1. Same gene: Variants in the same gene share functional context.
         A model might learn "gene X → pathogenic" rather than true biology.

      2. Linkage disequilibrium: Nearby variants on the same chromosome
         are co-inherited and correlated.

      3. Same ClinVar submission: Multiple variants submitted together
         may share clinical interpretation methodology.

      RECOMMENDED SPLITTING STRATEGIES (safest first):
      ─────────────────────────────────────────────────
      1. CHROMOSOME-BASED SPLIT (safest, recommended for first experiment)
         - Train: chromosomes 1-17
         - Validation: chromosomes 18-20
         - Test: chromosomes 21-22, X
         - Ensures zero overlap of genomic regions

      2. GENE-BASED SPLIT (second choice)
         - Split by gene: all variants in a gene go to the same set
         - Prevents gene-level leakage
         - Use GroupKFold with gene as the group

      3. STRATIFIED RANDOM SPLIT (least safe, but common)
         - 70/15/15 or 80/10/10 stratified by target
         - Only acceptable if combined with variant-key deduplication
         - Risk: related variants may appear in both train and test
    """)

    # ──────────────────────────────────────────────
    # 11. CLASS IMBALANCE
    # ──────────────────────────────────────────────
    section("11. CLASS IMBALANCE ANALYSIS")

    n_path = (df["target"] == "Pathogenic").sum()
    n_ben = (df["target"] == "Benign").sum()
    total = n_path + n_ben
    ratio = max(n_path, n_ben) / min(n_path, n_ben) if min(n_path, n_ben) > 0 else float("inf")

    report(f"  Pathogenic:  {n_path:>10,}  ({n_path/total*100:.2f}%)")
    report(f"  Benign:      {n_ben:>10,}  ({n_ben/total*100:.2f}%)")
    report(f"  Total:       {total:>10,}")
    report(f"  Imbalance ratio: {ratio:.2f}:1  ({'Pathogenic' if n_path > n_ben else 'Benign'} majority)")

    report("""
      RECOMMENDED IMBALANCE HANDLING:
      ───────────────────────────────
      1. CLASS WEIGHTS (safest, recommended)
         - Assign inverse-frequency weights during model training
         - Supported by sklearn, PyTorch, TensorFlow
         - No synthetic data created → no risk of artefacts

      2. STRATIFIED SPLITTING
         - Ensure train/val/test maintain the same class proportions
         - Always use regardless of other balancing methods

      3. FOCAL LOSS
         - For deep learning models, focal loss down-weights easy examples
         - Effective for moderate imbalance

      4. SMOTE (use with EXTREME caution)
         - Only apply AFTER train/test split (never on test data)
         - For genomic features, SMOTE interpolation may create
           biologically impossible variants
         - Only acceptable on numeric feature representations,
           NOT on raw allele sequences

      RECOMMENDATION: Start with class weights + stratified splitting.
      These are safe, simple, and biologically sound.
    """)

    # ──────────────────────────────────────────────
    # 12. SAVE ALL OUTPUTS
    # ──────────────────────────────────────────────
    section("12. SAVING OUTPUT FILES")

    # 12a. clinvar_cleaned.csv — full cleaned dataframe
    cleaned_path = os.path.join(OUTPUT_DIR, "clinvar_cleaned.csv")
    # Drop internal helper columns
    drop_internal = [c for c in ["_variant_key", "target_strict"] if c in df.columns]
    df_save = df.drop(columns=drop_internal, errors="ignore")
    df_save.to_csv(cleaned_path, index=False)
    report(f"  Saved: {cleaned_path} ({len(df_save):,} rows, {len(df_save.columns)} cols)")

    # 12b. clinvar_ml_dataset.csv — Dataset A
    ml_path = os.path.join(OUTPUT_DIR, "clinvar_ml_dataset.csv")
    dataset_a.to_csv(ml_path, index=False)
    report(f"  Saved: {ml_path} ({len(dataset_a):,} rows, {len(dataset_a.columns)} cols)")

    # 12c. removed_records.csv
    removed_path = os.path.join(OUTPUT_DIR, "removed_records.csv")
    if removed_records:
        df_removed = pd.concat(removed_records, ignore_index=True)
        # Ensure exclusion_reason column is present
        if "exclusion_reason" not in df_removed.columns:
            df_removed["exclusion_reason"] = "Unknown"
        df_removed.to_csv(removed_path, index=False)
        report(f"  Saved: {removed_path} ({len(df_removed):,} rows)")
    else:
        pd.DataFrame().to_csv(removed_path, index=False)
        report(f"  Saved: {removed_path} (0 rows)")

    # 12d. feature_dictionary.csv
    feat_dict_path = os.path.join(OUTPUT_DIR, "feature_dictionary.csv")
    feat_df = pd.DataFrame(feature_descriptions)
    feat_df.to_csv(feat_dict_path, index=False)
    report(f"  Saved: {feat_dict_path} ({len(feat_df)} features)")

    # 12e. clinvar_preprocessing_report.txt
    report_path = os.path.join(os.path.dirname(OUTPUT_DIR), "reports", "clinvar_preprocessing_report.txt")
    os.makedirs(os.path.dirname(report_path), exist_ok=True)

    # ──────────────────────────────────────────────
    # 13. FINAL QUALITY REPORT
    # ──────────────────────────────────────────────
    section("13. FINAL QUALITY REPORT")

    report(f"  Original row count:        {original_row_count:>12,}")
    report(f"  Final row count:           {len(df):>12,}")
    report(f"  Removed rows:              {original_row_count - len(df):>12,}")
    report(f"  Pathogenic variants:       {n_path:>12,}")
    report(f"  Benign variants:           {n_ben:>12,}")

    # Count excluded uncertain/other
    n_excluded_uncertain = 0
    if removed_records:
        df_rem_all = pd.concat(removed_records, ignore_index=True)
        n_excluded_uncertain = len(df_rem_all[
            df_rem_all.get("exclusion_reason", pd.Series()) == "Non-binary clinical significance"
        ]) if "exclusion_reason" in df_rem_all.columns else 0
    report(f"  Excluded/uncertain:        {n_excluded_uncertain:>12,}")

    n_dup_removed = 0
    if removed_records:
        for r in removed_records:
            if "exclusion_reason" in r.columns:
                n_dup_removed += r["exclusion_reason"].str.contains("uplicate|onflict", case=False, na=False).sum()
    report(f"  Duplicate-related removed:  {n_dup_removed:>12,}")

    report(f"\n  Final feature list: {ml_features}")
    report(f"  Potential leakage columns: {leakage_columns}")
    report(f"  Class distribution: Pathogenic={n_path:,} ({n_path/total*100:.1f}%), Benign={n_ben:,} ({n_ben/total*100:.1f}%)")
    report(f"  Recommended strategy: Chromosome-based split + class weights")

    end_time = datetime.now()
    report(f"\n  Script completed at: {end_time.isoformat()}")
    report(f"  Total runtime: {end_time - start_time}")

    # Write the report file
    with open(report_path, "w") as f:
        f.write("\n".join(REPORT_LINES))
    report(f"\n  Report saved to: {report_path}")

    # 12f. Save Dataset B stub
    dataset_b_path = os.path.join(OUTPUT_DIR, "clinvar_extended_dataset.csv")
    dataset_b.to_csv(dataset_b_path, index=False)
    report(f"  Saved: {dataset_b_path} ({len(dataset_b):,} rows, {len(dataset_b.columns)} cols)")

    print("\n" + "=" * 80)
    print("  PREPROCESSING COMPLETE")
    print("=" * 80)
    print(f"\n  All outputs saved to: ./{OUTPUT_DIR}/")
    print(f"  See {report_path} for the full preprocessing report.")

if __name__ == "__main__":
    main()
