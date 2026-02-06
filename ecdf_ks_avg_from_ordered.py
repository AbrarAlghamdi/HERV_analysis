#!/usr/bin/env python3

"""
ECDF + KS test on *_filtered.tsv, skipping broken/missing files.
- Searches current dir + telescope_ordered_Alph/
- Accepts 'transcript' or 'name'; accepts 'final_prop' or computes from 'final_count'
- Skips files that are missing, empty, fail to parse, or have too few rows
- Union-align transcripts (outer join), fill missing with 0
- Averages technical replicates per group (B1/B2/B3)
- Saves ECDF plot and pairwise KS results
"""

import os
import glob
import re
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

# ---------------- Config ----------------
search_roots = ["telescope_ordered_Alph", "."]
patterns = ["*_filtered.tsv"]
output_dir = "ecdf_ks_avg_results"
os.makedirs(output_dir, exist_ok=True)

# Skip criteria
min_rows_threshold = 10       # skip file if < 10 transcripts (adjust as needed)
drop_version_suffix = False   # normalize transcript IDs like "X.1" -> "X" if True

# Keep only transcripts present (non-zero) in >= K samples (set False to disable)
use_presence_filter = False
presence_K = 3

# ---------------- Helpers ----------------
def normalize_id(x: str) -> str:
    s = str(x).strip().replace('"', '')
    if drop_version_suffix and "." in s:
        s = s.split(".", 1)[0]
    return s

def read_filtered_table(path: str) -> pd.DataFrame:
    """Read a *_filtered.tsv and return DF indexed by transcript with one column 'final_prop'."""
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()

    # ID column
    if "transcript" in df.columns:
        id_col = "transcript"
    elif "name" in df.columns:
        id_col = "name"
        df.rename(columns={"name": "transcript"}, inplace=True)
    else:
        raise ValueError("Missing identifier column ('transcript' or 'name').")

    # Use final_prop if present; otherwise compute from final_count
    if "final_prop" in df.columns:
        out = df[["transcript", "final_prop"]].copy()
    else:
        if "final_count" not in df.columns:
            raise ValueError("Missing 'final_prop' and 'final_count'; need at least one.")
        tmp = df[["transcript", "final_count"]].copy()
        tot = tmp["final_count"].sum()
        tmp["final_prop"] = 0.0 if tot == 0 else tmp["final_count"] / tot
        out = tmp[["transcript", "final_prop"]].copy()

    out["transcript"] = out["transcript"].map(normalize_id)
    out = out.dropna(subset=["transcript"]).drop_duplicates(subset=["transcript"])
    out.set_index("transcript", inplace=True)
    out["final_prop"] = pd.to_numeric(out["final_prop"], errors="coerce").fillna(0.0)
    return out

def in_group(colname: str, group_num: int) -> bool:
    # Matches B1 / B_1 / B-1 / b1 / b_1 / b-1 at token boundaries
    return re.search(fr'(^|[^A-Za-z0-9])B[_-]?{group_num}([^A-Za-z0-9]|$)',
                     colname, flags=re.IGNORECASE) is not None

def ecdf(data: np.ndarray) -> np.ndarray:
    if len(data) == 0:
        return np.array([])
    return np.arange(1, len(data) + 1) / len(data)

# ---------------- Find files ----------------
candidates = []
for root in search_roots:
    for pat in patterns:
        candidates.extend(glob.glob(os.path.join(root, pat)))
candidates = sorted(set(candidates))

print(f"📄 Found {len(candidates)} candidate files")
for f in candidates:
    print("  -", f)
if not candidates:
    raise RuntimeError(" No *_filtered.tsv files found.")

# ---------------- Load & skip broken ----------------
dfs = {}
skipped = []   # (file, reason)

for fp in candidates:
    sample = os.path.splitext(os.path.basename(fp))[0].replace("_filtered", "")
    try:
        if not os.path.exists(fp):
            skipped.append((fp, "missing (path does not exist)"))
            continue
        if os.path.getsize(fp) == 0:
            skipped.append((fp, "empty file (0 bytes)"))
            continue

        df = read_filtered_table(fp)

        if df.shape[0] < min_rows_threshold:
            skipped.append((fp, f"too few rows (<{min_rows_threshold})"))
            continue

        dfs[sample] = df
        print(f" Loaded {sample}: {df.shape[0]} transcripts")

    except Exception as e:
        skipped.append((fp, f"read/parse error: {e}"))

if not dfs:
    raise RuntimeError(" All files were skipped as broken/empty/invalid. Nothing to analyze.")

# ---------------- Diagnostics: overlap ----------------
samples = list(dfs.keys())
for s in samples:
    print(f"• {s}: {dfs[s].shape[0]} unique transcripts")
for s1, s2 in itertools.combinations(samples, 2):
    a = set(dfs[s1].index); b = set(dfs[s2].index)
    inter, union = len(a & b), len(a | b)
    jac = inter / union if union else 0
    print(f"   Overlap {s1} vs {s2}: intersect={inter}, union={union}, Jaccard={jac:.3f}")

# ---------------- UNION alignment ----------------
combined_df = pd.concat([dfs[s] for s in samples], axis=1, join="outer")
combined_df.columns = samples
combined_df = combined_df.fillna(0.0).sort_index()
print("Combined matrix shape (rows x samples):", combined_df.shape)

# Optional ≥K presence filter
if use_presence_filter:
    present_counts = (combined_df > 0).astype(int).sum(axis=1)
    keep = present_counts >= presence_K
    kept_rows = int(keep.sum())
    print(f"🔎 Presence filter: keep rows present in ≥ {presence_K} samples -> {kept_rows} transcripts")
    if kept_rows == 0:
        raise RuntimeError(" No transcripts meet the ≥K presence criterion.")
    combined_df = combined_df.loc[keep]

# ---------------- Group detection ----------------
candidate_groups = ["B1", "B2", "B3"]
grouped = {g: [c for c in combined_df.columns if in_group(c, int(g[1:]))] for g in candidate_groups}
present = {g: cols for g, cols in grouped.items() if cols}

print("Groups detected (column counts):", {g: len(v) for g, v in present.items()})
if len(present) < 2:
    raise RuntimeError(f" Need at least two groups among {candidate_groups} after skipping broken files. "
                       f"Found: { {g: len(v) for g,v in present.items()} }")

# ---------------- Average tech replicates per group ----------------
group_avgs_sorted = {}
for g, cols in present.items():
    series = combined_df[cols].mean(axis=1)
    arr = np.sort(series.values)
    group_avgs_sorted[g] = arr

# ---------------- Plot ECDFs ----------------
plt.figure(figsize=(8, 6))
for g, arr in group_avgs_sorted.items():
    plt.plot(arr, ecdf(arr), lw=2, label=g)
plt.xlabel("final_prop (HERV Expression Proportion)")
plt.ylabel("Cumulative Probability")
plt.title("Empirical CDF - Kolmogorov–Smirnov Test (group means, union-aligned)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plot_path = os.path.join(output_dir, "ecdf_ks_avg_plot.png")
plt.savefig(plot_path, dpi=150)
print(f" Saved ECDF plot: {plot_path}")

# ---------------- Pairwise KS tests ----------------
rows = []
for g1, g2 in itertools.combinations(group_avgs_sorted.keys(), 2):
    res = ks_2samp(group_avgs_sorted[g1], group_avgs_sorted[g2])
    rows.append({"group1": g1, "group2": g2, "ks_stat": res.statistic, "p_value": res.pvalue})
    print(f"KS {g1} vs {g2}: stat={res.statistic:.4f}, p={res.pvalue:.4g}")

ks_df = pd.DataFrame(rows)
ks_out = os.path.join(output_dir, "ks_results.tsv")
ks_df.to_csv(ks_out, sep="\t", index=False)
print(f"📄 Saved KS table: {ks_out}")

# ---------------- Report skipped files ----------------
if skipped:
    print("\n Skipped files:")
    for fp, reason in skipped:
        print(f"  - {fp}: {reason}")
