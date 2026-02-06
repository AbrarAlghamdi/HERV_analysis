#!/usr/bin/env python3

"""
Spearman correlation (group means) on *_filtered.tsv, skipping broken/missing files.
- Searches current dir + telescope_ordered_Alph/
- Accepts 'transcript' or 'name'; accepts 'final_prop' or computes from 'final_count'
- Skips files that are missing, empty, fail to parse, or have too few rows
- Union-align transcripts (outer join), fill missing with 0
- Requires >= 2 groups (B1/B2/B3) after skipping; averages tech replicates
- Heatmap shows numbers
"""

import os
import glob
import re
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# ---------------- Config ----------------
search_roots = ["telescope_ordered_Alph", "."]
patterns = ["*_filtered.tsv"]
output_dir = "spearman_grouped_results"
os.makedirs(output_dir, exist_ok=True)

# Skip criteria
min_rows_threshold = 10       # skip file if it has < 10 transcripts (tune as you like)
drop_version_suffix = False   # normalize transcript IDs like "X.1" -> "X" if True
use_presence_filter = False   # keep only transcripts present in >= K samples
presence_K = 3

# ---------------- Helpers ----------------
def normalize_id(x: str) -> str:
    s = str(x).strip().replace('"', '')
    if drop_version_suffix and "." in s:
        s = s.split(".", 1)[0]
    return s

def read_filtered_table(path: str) -> pd.DataFrame:
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
    # final_prop or compute from final_count
    if "final_prop" in df.columns:
        out = df[["transcript", "final_prop"]].copy()
    else:
        if "final_count" not in df.columns:
            raise ValueError("Missing 'final_prop' and 'final_count'; need at least one.")
        tmp = df[["transcript", "final_count"]].copy()
        tot = tmp["final_count"].sum()
        tmp["final_prop"] = 0.0 if tot == 0 else tmp["final_count"] / tot
        out = tmp[["transcript", "final_prop"]].copy()
    # clean
    out["transcript"] = out["transcript"].map(normalize_id)
    out = out.dropna(subset=["transcript"]).drop_duplicates(subset=["transcript"])
    out.set_index("transcript", inplace=True)
    out["final_prop"] = pd.to_numeric(out["final_prop"], errors="coerce").fillna(0.0)
    return out

def in_group(colname: str, group_num: int) -> bool:
    return re.search(fr'(^|[^A-Za-z0-9])B[_-]?{group_num}([^A-Za-z0-9]|$)',
                     colname, flags=re.IGNORECASE) is not None

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
    raise RuntimeError("❌ No *_filtered.tsv files found.")

# ---------------- Load & skip broken ----------------
dfs = {}
skipped = []   # (file, reason)

for fp in candidates:
    sample = os.path.splitext(os.path.basename(fp))[0].replace("_filtered", "")
    try:
        # missing or empty file
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
        print(f"Loaded {sample}: {df.shape[0]} transcripts")

    except Exception as e:
        skipped.append((fp, f"read/parse error: {e}"))

if not dfs:
    raise RuntimeError("❌ All files were skipped as broken/empty/invalid. Nothing to analyze.")

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
combined = pd.concat([dfs[s] for s in samples], axis=1, join="outer")
combined.columns = samples
combined = combined.fillna(0.0).sort_index()
print(" Combined matrix shape (rows x samples):", combined.shape)

# Optional presence filter
if use_presence_filter:
    present_counts = (combined > 0).astype(int).sum(axis=1)
    keep = present_counts >= presence_K
    kept_rows = int(keep.sum())
    print(f"🔎 Presence filter: keep rows present in ≥ {presence_K} samples -> {kept_rows} transcripts")
    if kept_rows == 0:
        raise RuntimeError(" No transcripts meet the ≥K presence criterion.")
    combined = combined.loc[keep]

# ---------------- Group detection ----------------
candidate_groups = ["B1", "B2", "B3"]
grouped = {g: [c for c in combined.columns if in_group(c, int(g[1:]))] for g in candidate_groups}
present = {g: cols for g, cols in grouped.items() if cols}

print(" Groups detected (column counts):", {g: len(v) for g, v in present.items()})
if len(present) < 2:
    raise RuntimeError(f"❌ Need at least two groups among {candidate_groups} after skipping broken files. "
                       f"Found: { {g: len(v) for g,v in present.items()} }")

# ---------------- Average tech replicates per group ----------------
group_means = {g: combined[cols].mean(axis=1) for g, cols in present.items()}
groups = list(group_means.keys())

# ---------------- Spearman correlation ----------------
corr_mat = pd.DataFrame(index=groups, columns=groups, dtype=float)
pval_mat = pd.DataFrame(index=groups, columns=groups, dtype=float)

for g1, g2 in itertools.product(groups, groups):
    rho, p = spearmanr(group_means[g1], group_means[g2])
    corr_mat.loc[g1, g2] = rho
    pval_mat.loc[g1, g2] = p

# ---------------- Save tables ----------------
os.makedirs(output_dir, exist_ok=True)
corr_path = os.path.join(output_dir, "spearman_corr_matrix.tsv")
pval_path = os.path.join(output_dir, "spearman_pvalues_matrix.tsv")
corr_mat.to_csv(corr_path, sep="\t")
pval_mat.to_csv(pval_path, sep="\t")
print(f" Saved Spearman correlation matrix: {corr_path}")
print(f" Saved p-value matrix: {pval_path}")

# ---------------- Heatmap with numbers ----------------
fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(corr_mat.values.astype(float), vmin=-1, vmax=1)
cbar = plt.colorbar(im, ax=ax, label="Spearman ρ")

ax.set_xticks(range(len(groups)))
ax.set_xticklabels(groups, rotation=45, ha="right")
ax.set_yticks(range(len(groups)))
ax.set_yticklabels(groups)
ax.set_title("Spearman correlation (group means of final_prop)")

vals = corr_mat.values.astype(float)
for i in range(len(groups)):
    for j in range(len(groups)):
        v = vals[i, j]
        txt_color = "white" if abs(v) > 0.5 else "black"
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", color=txt_color, fontsize=11)

plt.tight_layout()
plot_path = os.path.join(output_dir, "spearman_corr_heatmap.png")
plt.savefig(plot_path, dpi=150)
print(f" Saved heatmap: {plot_path}")

# ---------------- Report skipped files ----------------
if skipped:
    print("\n Skipped files:")
    for fp, reason in skipped:
        print(f"  - {fp}: {reason}")
