#!/usr/bin/env python3
import pandas as pd

# =========================
# INPUT FILES
# =========================
CIS_FILE = "DESeq2_HERV_Cisplatin_vs_Untreated.csv"
CARB_FILE = "DESeq2_HERV_Carboplatin_vs_Untreated.csv"

FDR_THRESHOLD = 0.05

# =========================
# LOAD DATA
# =========================
cis = pd.read_csv(CIS_FILE)
carb = pd.read_csv(CARB_FILE)

# Keep necessary columns
cis = cis[["transcript", "log2FoldChange", "padj"]].rename(columns={
    "log2FoldChange": "log2FC_Cis",
    "padj": "padj_Cis"
})

carb = carb[["transcript", "log2FoldChange", "padj"]].rename(columns={
    "log2FoldChange": "log2FC_Carb",
    "padj": "padj_Carb"
})

# Merge by transcript
df = pd.merge(cis, carb, on="transcript", how="outer")

# Ensure numeric
df["log2FC_Cis"] = pd.to_numeric(df["log2FC_Cis"], errors="coerce")
df["padj_Cis"] = pd.to_numeric(df["padj_Cis"], errors="coerce")
df["log2FC_Carb"] = pd.to_numeric(df["log2FC_Carb"], errors="coerce")
df["padj_Carb"] = pd.to_numeric(df["padj_Carb"], errors="coerce")

# =========================
# SIGNIFICANCE CLASSIFICATION
# =========================
def classify_significance(row):
    cis_sig = row["padj_Cis"] < FDR_THRESHOLD if pd.notna(row["padj_Cis"]) else False
    carb_sig = row["padj_Carb"] < FDR_THRESHOLD if pd.notna(row["padj_Carb"]) else False

    if cis_sig and not carb_sig:
        return "Cisplatin-only significant"
    elif carb_sig and not cis_sig:
        return "Carboplatin-only significant"
    elif cis_sig and carb_sig:
        return "Shared significant loci"
    else:
        return "Non-significant loci"

df["Significance_Category"] = df.apply(classify_significance, axis=1)

# =========================
# DIRECTIONAL CLASSIFICATION
# =========================
def direction_label(log2fc):
    if pd.isna(log2fc):
        return "NA"
    elif log2fc > 0:
        return "Up"
    elif log2fc < 0:
        return "Down"
    else:
        return "Zero"

df["Direction_Cis"] = df["log2FC_Cis"].apply(direction_label)
df["Direction_Carb"] = df["log2FC_Carb"].apply(direction_label)

# Combined direction (only if both directions exist)
def combined_direction(row):
    if row["Direction_Cis"] in ["Up", "Down"] and row["Direction_Carb"] in ["Up", "Down"]:
        return f"{row['Direction_Cis']}_{row['Direction_Carb']}"
    else:
        return "NA"

df["Combined_Direction"] = df.apply(combined_direction, axis=1)

# =========================
# SUMMARY TABLES
# =========================

print("\n===== Significance Categories =====")
print(df["Significance_Category"].value_counts())

print("\n===== Combined Direction Categories =====")
print(df["Combined_Direction"].value_counts())

# Save output
df.to_csv("HERV_Cis_Carb_Directional_Classification.tsv", sep="\t", index=False)

print("\nSaved file:")
print(" - HERV_Cis_Carb_Directional_Classification.tsv")
