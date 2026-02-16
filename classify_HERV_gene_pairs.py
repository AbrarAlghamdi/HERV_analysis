#!/usr/bin/env python3
import pandas as pd

# =========================
# Inputs
# =========================
HERV_DE = "DESeq2_HERV_Carboplatin_vs_Untreated.csv"
GENE_DE = "DESeq2_gene_Carboplatin_Treated_vs_Untreated_full.csv"
MAPPING = "herv_nearest_gene_10kb.tsv"  # bedtools closest -d output filtered to <=10kb

HERV_PADJ_CUTOFF = 0.05
GENE_PADJ_CUTOFF = 0.05

# =========================
# Load DE results
# =========================
herv = pd.read_csv(HERV_DE)
gene = pd.read_csv(GENE_DE)

# Validate HERV cols
required_herv_cols = {"transcript", "log2FoldChange", "padj"}
missing_herv = required_herv_cols - set(herv.columns)
if missing_herv:
    raise ValueError(f"HERV DE file missing columns: {missing_herv}")

# Identify gene ID column
gene_id_col = None
if "gene_id" in gene.columns:
    gene_id_col = "gene_id"
elif "ENSEMBL" in gene.columns:
    gene_id_col = "ENSEMBL"
else:
    raise ValueError("Gene DE file must contain 'gene_id' or 'ENSEMBL' column.")

required_gene_cols = {gene_id_col, "log2FoldChange", "padj"}
missing_gene = required_gene_cols - set(gene.columns)
if missing_gene:
    raise ValueError(f"Gene DE file missing columns: {missing_gene}")

# Clean Ensembl IDs: remove version suffix (e.g., ENSG... .14)
gene["gene_id_clean"] = gene[gene_id_col].astype(str).str.split(".").str[0]

# =========================
# Load mapping (whitespace-separated)
# =========================
map_df = pd.read_csv(MAPPING, sep=r"\s+", header=None)

# Expected columns from bedtools closest -d
expected_ncols = 13
if map_df.shape[1] != expected_ncols:
    raise ValueError(
        f"Unexpected number of columns in mapping file: {map_df.shape[1]} (expected {expected_ncols}).\n"
        f"Tip: run `head {MAPPING} | cat -A` to check delimiters."
    )

map_df.columns = [
    "h_chr", "h_start", "h_end", "transcript", "score", "h_strand",
    "g_chr", "g_start", "g_end", "gene_name", "gene_id", "g_strand",
    "distance_bp"
]

# Clean gene IDs in mapping too
map_df["gene_id_clean"] = map_df["gene_id"].astype(str).str.split(".").str[0]

# =========================
# Merge DE onto mapping
# =========================
map_df = map_df.merge(
    herv[["transcript", "log2FoldChange", "padj"]],
    on="transcript",
    how="left"
).rename(columns={
    "log2FoldChange": "herv_log2FC",
    "padj": "herv_padj"
})

map_df = map_df.merge(
    gene[["gene_id_clean", "log2FoldChange", "padj"]],
    on="gene_id_clean",
    how="left"
).rename(columns={
    "log2FoldChange": "gene_log2FC",
    "padj": "gene_padj"
})

# Ensure numeric
map_df["herv_log2FC"] = pd.to_numeric(map_df["herv_log2FC"], errors="coerce")
map_df["gene_log2FC"] = pd.to_numeric(map_df["gene_log2FC"], errors="coerce")
map_df["herv_padj"] = pd.to_numeric(map_df["herv_padj"], errors="coerce")
map_df["gene_padj"] = pd.to_numeric(map_df["gene_padj"], errors="coerce")

# =========================
# Classification
# =========================
def classify(row) -> str:
    # Gene not significant or missing -> Independent
    if pd.isna(row["gene_padj"]) or row["gene_padj"] >= GENE_PADJ_CUTOFF:
        return "Independent"

    # HERV should be significant (because we filtered), but keep safe
    if pd.isna(row["herv_padj"]) or row["herv_padj"] >= HERV_PADJ_CUTOFF:
        return "HERV_not_significant"

    h = row["herv_log2FC"]
    g = row["gene_log2FC"]
    if pd.isna(h) or pd.isna(g):
        return "Missing_log2FC"

    if h > 0 and g > 0:
        return "Up–Up"
    if h < 0 and g < 0:
        return "Down–Down"
    if h > 0 and g < 0:
        return "Up–Down"
    if h < 0 and g > 0:
        return "Down–Up"
    return "Other"

map_df["Regulatory_Pattern"] = map_df.apply(classify, axis=1)

# =========================
# Summary table
# =========================
summary = (map_df["Regulatory_Pattern"]
           .value_counts()
           .rename_axis("Regulatory Pattern")
           .reset_index(name="Number of Pairs"))

total = int(summary["Number of Pairs"].sum())
summary["Percentage"] = (summary["Number of Pairs"] / total * 100).round(1)

print(f"\nTotal pairs within ±10kb: {total}\n")
print("Summary:")
print(summary.to_string(index=False))

# Save outputs
map_df.to_csv("HERV_gene_pairs_classified.tsv", sep="\t", index=False)
summary.to_csv("HERV_gene_regulatory_summary.tsv", sep="\t", index=False)

print("\nFiles written:")
print(" - HERV_gene_pairs_classified.tsv")
print(" - HERV_gene_regulatory_summary.tsv")
