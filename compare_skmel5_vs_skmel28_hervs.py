#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ------------------ INPUTS ------------------
FILE_5  = "differential_expression_ttest_fdr_5.tsv"    # SK-MEL-5 results (1% vs 10%)
FILE_28 = "differential_expression_ttest_fdr_28.tsv"   # SK-MEL-28 results (1% vs 10%)

OUTDIR = "./SKMEL5_vs_SKMEL28_comparison"
FDR_CUTOFF = 0.05

# Column names expected in your result files:
COL_TX   = "transcript"
COL_L2FC = "log2_fold_change_01_over_10"
COL_FDR  = "fdr_adjusted_p"
COL_P    = "p_value"  # optional, but nice to keep
# --------------------------------------------

os.makedirs(OUTDIR, exist_ok=True)

def load_de_table(path: str, label: str) -> pd.DataFrame:
    """Load DE table, validate columns, add significance flag."""
    df = pd.read_csv(path, sep="\t")

    # Basic validation
    needed = {COL_TX, COL_L2FC, COL_FDR}
    missing = needed - set(df.columns)
    if missing:
        raise RuntimeError(f"{label}: missing columns {missing} in {path}\nFound: {list(df.columns)}")

    # Clean
    df = df.copy()
    df[COL_TX] = df[COL_TX].astype(str).str.strip()
    df[COL_L2FC] = pd.to_numeric(df[COL_L2FC], errors="coerce")
    df[COL_FDR] = pd.to_numeric(df[COL_FDR], errors="coerce")

    df["significant"] = df[COL_FDR] < FDR_CUTOFF
    return df

# ---------- Load ----------
de5  = load_de_table(FILE_5,  "SK-MEL-5")
de28 = load_de_table(FILE_28, "SK-MEL-28")

# ---------- Sets: significant loci ----------
sig5_set  = set(de5.loc[de5["significant"], COL_TX])
sig28_set = set(de28.loc[de28["significant"], COL_TX])

shared_sig = sig5_set & sig28_set
only5_sig  = sig5_set - sig28_set
only28_sig = sig28_set - sig5_set

# ---------- Summary ----------
summary = pd.DataFrame({
    "category": ["SK-MEL-5 significant", "SK-MEL-28 significant", "Shared significant", "SK-MEL-5 only", "SK-MEL-28 only"],
    "count":    [len(sig5_set),           len(sig28_set),          len(shared_sig),       len(only5_sig),  len(only28_sig)]
})

summary_path = os.path.join(OUTDIR, "summary_counts.tsv")
summary.to_csv(summary_path, sep="\t", index=False)

print("\n=== Summary (FDR < {:.2g}) ===".format(FDR_CUTOFF))
print(summary.to_string(index=False))
print("\nSaved:", summary_path)

# ---------- Merge tables for shared significant loci ----------
de5_small = de5[[COL_TX, COL_L2FC, COL_FDR] + ([COL_P] if COL_P in de5.columns else [])].copy()
de28_small = de28[[COL_TX, COL_L2FC, COL_FDR] + ([COL_P] if COL_P in de28.columns else [])].copy()

de5_small = de5_small.rename(columns={
    COL_L2FC: "log2FC_skmel5",
    COL_FDR: "FDR_skmel5",
    **({COL_P: "p_skmel5"} if COL_P in de5_small.columns else {})
})
de28_small = de28_small.rename(columns={
    COL_L2FC: "log2FC_skmel28",
    COL_FDR: "FDR_skmel28",
    **({COL_P: "p_skmel28"} if COL_P in de28_small.columns else {})
})

shared_df = pd.merge(de5_small, de28_small, on=COL_TX, how="inner")
shared_df_sig = shared_df[shared_df[COL_TX].isin(shared_sig)].copy()

# Direction category for shared significant loci
def direction(row):
    a = row["log2FC_skmel5"]
    b = row["log2FC_skmel28"]
    if pd.isna(a) or pd.isna(b):
        return "unknown"
    if a > 0 and b > 0:
        return "concordant_up"
    if a < 0 and b < 0:
        return "concordant_down"
    if a > 0 and b < 0:
        return "opposite"
    if a < 0 and b > 0:
        return "opposite"
    return "zero_or_unknown"

shared_df_sig["direction_class"] = shared_df_sig.apply(direction, axis=1)
shared_df_sig["abs_log2FC_skmel5"] = shared_df_sig["log2FC_skmel5"].abs()
shared_df_sig["abs_log2FC_skmel28"] = shared_df_sig["log2FC_skmel28"].abs()

shared_out = os.path.join(OUTDIR, "shared_significant_hervs.tsv")
shared_df_sig.sort_values(["direction_class", "FDR_skmel5", "FDR_skmel28"]).to_csv(shared_out, sep="\t", index=False)
print("Saved shared significant table:", shared_out)

# ---------- SK-MEL-5-only and SK-MEL-28-only tables ----------
only5_df = de5[de5[COL_TX].isin(only5_sig)].copy().sort_values(COL_FDR)
only28_df = de28[de28[COL_TX].isin(only28_sig)].copy().sort_values(COL_FDR)

only5_out = os.path.join(OUTDIR, "skmel5_only_significant_hervs.tsv")
only28_out = os.path.join(OUTDIR, "skmel28_only_significant_hervs.tsv")

only5_df.to_csv(only5_out, sep="\t", index=False)
only28_df.to_csv(only28_out, sep="\t", index=False)

print("Saved SK-MEL-5 only table:", only5_out)
print("Saved SK-MEL-28 only table:", only28_out)

# ---------- Direction counts ----------
dir_counts = shared_df_sig["direction_class"].value_counts().reset_index()
dir_counts.columns = ["direction_class", "count"]
dir_out = os.path.join(OUTDIR, "shared_direction_counts.tsv")
dir_counts.to_csv(dir_out, sep="\t", index=False)
print("Saved direction counts:", dir_out)

# ------------------ PLOTS ------------------

# 1) Bar plot of summary counts
plt.figure(figsize=(8, 5))
plt.bar(summary["category"], summary["count"])
plt.xticks(rotation=30, ha="right")
plt.ylabel("Number of HERV loci")
plt.title("Stress-responsive HERV loci (FDR < {:.2g})".format(FDR_CUTOFF))
plt.tight_layout()
p1 = os.path.join(OUTDIR, "summary_counts_barplot.png")
plt.savefig(p1, dpi=150)
print("Saved plot:", p1)
plt.close()

# 2) Scatter plot log2FC SK-MEL-5 vs SK-MEL-28 for shared significant
if len(shared_df_sig) > 0:
    plt.figure(figsize=(6, 6))
    x = shared_df_sig["log2FC_skmel5"].values
    y = shared_df_sig["log2FC_skmel28"].values
    plt.scatter(x, y, s=18, alpha=0.8)
    plt.axhline(0, linestyle="--")
    plt.axvline(0, linestyle="--")
    plt.xlabel("log2FC (SK-MEL-5: 1% / 10%)")
    plt.ylabel("log2FC (SK-MEL-28: 1% / 10%)")
    plt.title("Shared significant HERVs: log2FC comparison")
    plt.tight_layout()
    p2 = os.path.join(OUTDIR, "shared_log2fc_scatter.png")
    plt.savefig(p2, dpi=150)
    print("Saved plot:", p2)
    plt.close()

# 3) Direction class bar plot
plt.figure(figsize=(6, 4))
plt.bar(dir_counts["direction_class"], dir_counts["count"])
plt.xticks(rotation=20, ha="right")
plt.ylabel("Number of shared significant loci")
plt.title("Direction of regulation (shared significant HERVs)")
plt.tight_layout()
p3 = os.path.join(OUTDIR, "shared_direction_barplot.png")
plt.savefig(p3, dpi=150)
print("Saved plot:", p3)
plt.close()

print("\nAll done ")
