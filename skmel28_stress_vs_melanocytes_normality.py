#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

# ------------------ INPUTS ------------------
SK28_01_MATRIX = "merged_final_count_matrix_common_finalcount1.tsv"     # SK-MEL-28 1% FBS
SK28_10_MATRIX = "merged_final_count_matrix_common_finalcount10.tsv"    # SK-MEL-28 10% FBS

# Normal melanocytes merged matrix (all normal samples columns)
MEL_MATRIX = "melanocytes_merged_finalcount.tsv"

# SK-MEL-28 DE results (from your Welch t-test + BH-FDR)
SK28_DE_TABLE = "differential_expression_ttest_fdr_28.tsv"

OUTDIR = "./SKMEL28_vs_Melanocytes_normality"
FDR_CUTOFF = 0.05
EPS = 1e-12
# -------------------------------------------
os.makedirs(OUTDIR, exist_ok=True)

def load_matrix(path: str) -> pd.DataFrame:
    """
    Load a merged final_count matrix:
      - rows: transcripts (index)
      - columns: samples
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return df

def counts_to_prop(df: pd.DataFrame) -> pd.DataFrame:
    """Convert counts to within-sample proportions (each column sums to 1)."""
    col_sums = df.sum(axis=0)
    return df.div(col_sums.replace(0, np.nan), axis=1).fillna(0.0)

# ---------- Load matrices ----------
sk01 = load_matrix(SK28_01_MATRIX)
sk10 = load_matrix(SK28_10_MATRIX)
mel  = load_matrix(MEL_MATRIX)

# ---------- Load DE table ----------
de28 = pd.read_csv(SK28_DE_TABLE, sep="\t")
de28["fdr_adjusted_p"] = pd.to_numeric(de28["fdr_adjusted_p"], errors="coerce")
de28["transcript"] = de28["transcript"].astype(str).str.strip()

# Stress-responsive set in SK-MEL-28
stress_tx = set(de28.loc[de28["fdr_adjusted_p"] < FDR_CUTOFF, "transcript"])
print(f"Stress-responsive (SK-MEL-28, FDR<{FDR_CUTOFF}): {len(stress_tx)}")

# ---------- Align transcript universe ----------
common = sk01.index.intersection(sk10.index).intersection(mel.index)
if len(common) == 0:
    raise RuntimeError("No overlapping transcripts among SK-MEL-28 (1%,10%) and melanocytes matrices.")

sk01 = sk01.loc[common]
sk10 = sk10.loc[common]
mel  = mel.loc[common]

# Keep only stress-responsive that exist in all 3
stress_common = sorted(list(set(common) & stress_tx))
print("Stress-responsive present in all three datasets:", len(stress_common))

# ---------- Convert to proportions ----------
sk01p = counts_to_prop(sk01)
sk10p = counts_to_prop(sk10)
melp  = counts_to_prop(mel)

# Group means (mean proportion per transcript)
sk01_mean = sk01p.mean(axis=1)
sk10_mean = sk10p.mean(axis=1)
mel_mean  = melp.mean(axis=1)

# ---------- Per-transcript distance to normal ----------
dist_10 = (sk10_mean - mel_mean).abs()
dist_01 = (sk01_mean - mel_mean).abs()
delta_dist = dist_01 - dist_10   # <0 means moved closer to normal under stress

# Build results table (only stress-responsive loci)
res = pd.DataFrame({
    "transcript": stress_common,
    "mel_mean_prop": mel_mean.loc[stress_common].values,
    "skmel28_10_mean_prop": sk10_mean.loc[stress_common].values,
    "skmel28_01_mean_prop": sk01_mean.loc[stress_common].values,
    "dist_to_mel_10": dist_10.loc[stress_common].values,
    "dist_to_mel_01": dist_01.loc[stress_common].values,
    "delta_dist_(01-10)": delta_dist.loc[stress_common].values
})

# Classification
res["stress_effect_vs_normal"] = np.where(
    res["delta_dist_(01-10)"] < 0, "more_normal_under_stress",
    np.where(res["delta_dist_(01-10)"] > 0, "more_abnormal_under_stress", "no_change")
)

# ---------- Save outputs ----------
out_table = os.path.join(OUTDIR, "skmel28_stress_responsive_vs_melanocytes.tsv")
res.sort_values("delta_dist_(01-10)").to_csv(out_table, sep="\t", index=False)
print("Saved table:", out_table)

# ---------- Summaries ----------
n_total = res.shape[0]
n_more_normal = (res["stress_effect_vs_normal"] == "more_normal_under_stress").sum()
n_more_abnormal = (res["stress_effect_vs_normal"] == "more_abnormal_under_stress").sum()
n_no = (res["stress_effect_vs_normal"] == "no_change").sum()

summary = pd.DataFrame({
    "category": ["more_normal_under_stress", "more_abnormal_under_stress", "no_change", "total"],
    "count": [n_more_normal, n_more_abnormal, n_no, n_total],
    "fraction": [
        n_more_normal / n_total if n_total else np.nan,
        n_more_abnormal / n_total if n_total else np.nan,
        n_no / n_total if n_total else np.nan,
        1.0
    ]
})

summary_path = os.path.join(OUTDIR, "summary_normality_shift.tsv")
summary.to_csv(summary_path, sep="\t", index=False)
print("Saved summary:", summary_path)
print(summary.to_string(index=False))

# ---------- Simple paired test (global tendency) ----------
# Wilcoxon signed-rank: compare dist_01 vs dist_10 across loci
if n_total >= 10:
    stat, p = wilcoxon(res["dist_to_mel_01"], res["dist_to_mel_10"], alternative="two-sided")
    with open(os.path.join(OUTDIR, "wilcoxon_test.txt"), "w") as fh:
        fh.write("Wilcoxon signed-rank test (dist_01 vs dist_10) on stress-responsive loci\n")
        fh.write(f"n = {n_total}\nstat = {stat}\np = {p}\n")
    print(f"Wilcoxon test saved (p={p:.4g})")
else:
    print("Not enough loci for Wilcoxon test (need ~10+).")

# ---------- Plot: dist to melanocytes (10% vs 1%) ----------
plt.figure(figsize=(6, 6))
plt.scatter(res["dist_to_mel_10"], res["dist_to_mel_01"], s=18, alpha=0.8)
mx = max(res["dist_to_mel_10"].max(), res["dist_to_mel_01"].max()) if n_total else 1
plt.plot([0, mx], [0, mx], linestyle="--")  # diagonal: no change

plt.xlabel("Distance to melanocytes (SK-MEL-28 10%)")
plt.ylabel("Distance to melanocytes (SK-MEL-28 1%)")
plt.title("Do stress-responsive HERVs move toward normal? (SK-MEL-28)")
plt.tight_layout()

plot_path = os.path.join(OUTDIR, "distance_to_melanocytes_scatter.png")
plt.savefig(plot_path, dpi=150)
print("Saved plot:", plot_path)

print("Done ")
