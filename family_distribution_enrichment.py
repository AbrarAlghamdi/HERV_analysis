#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

INFILE = "HERV_Cis_Carb_Directional_Classification.tsv"
TOP_N = 15               # show top families in plots
DPI = 300

# -----------------------------
# Helper: extract "family"
# -----------------------------
def extract_family(transcript: str) -> str:
    """
    Extract family label from transcript name.
    Default: take text before first underscore.
    Examples:
      HERVH_12p12.2b -> HERVH
      ERVLE_1p36.13b -> ERVLE
      HARLEQUIN_1p36.33 -> HARLEQUIN
    """
    if pd.isna(transcript):
        return "NA"
    t = str(transcript)
    return t.split("_")[0] if "_" in t else t

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(INFILE, sep="\t")

# Expect columns produced by your earlier script:
# transcript, padj_Cis, padj_Carb, Significance_Category
required_cols = {"transcript", "Significance_Category"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in {INFILE}: {missing}")

df["Family"] = df["transcript"].apply(extract_family)

# Categories
cis_only = df[df["Significance_Category"] == "Cisplatin-only significant"].copy()
carb_only = df[df["Significance_Category"] == "Carboplatin-only significant"].copy()
shared = df[df["Significance_Category"] == "Shared significant loci"].copy()

print("Counts:")
print("  Cisplatin-only:", len(cis_only))
print("  Carboplatin-only:", len(carb_only))
print("  Shared:", len(shared))
print("  Total:", len(df))

# -----------------------------
# Family distribution tables
# -----------------------------
def family_distribution(subset: pd.DataFrame, label: str) -> pd.DataFrame:
    counts = subset["Family"].value_counts().reset_index()
    counts.columns = ["Family", "Count"]
    total = counts["Count"].sum()
    counts["Percentage"] = (counts["Count"] / total * 100).round(1)
    counts["Group"] = label
    return counts

cis_dist = family_distribution(cis_only, "Cisplatin-only")
carb_dist = family_distribution(carb_only, "Carboplatin-only")
shared_dist = family_distribution(shared, "Shared")

cis_dist.to_csv("Cisplatin_only_family_distribution.tsv", sep="\t", index=False)
carb_dist.to_csv("Carboplatin_only_family_distribution.tsv", sep="\t", index=False)
shared_dist.to_csv("Shared_family_distribution.tsv", sep="\t", index=False)

print("\nSaved family distribution tables:")
print(" - Cisplatin_only_family_distribution.tsv")
print(" - Carboplatin_only_family_distribution.tsv")
print(" - Shared_family_distribution.tsv")

# -----------------------------
# Plot function (single plot)
# -----------------------------
def plot_top_families(dist_df: pd.DataFrame, title: str, outpng: str):
    top = dist_df.sort_values("Count", ascending=False).head(TOP_N).copy()
    # reverse so largest appears at top if horizontal
    top = top.sort_values("Count", ascending=True)

    plt.figure(figsize=(8, 5))
    plt.barh(top["Family"], top["Count"])
    plt.xlabel("Number of significant loci")
    plt.ylabel("HERV family")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng, dpi=DPI)
    plt.show()
    print("Saved:", outpng)

# Cis-only plot
plot_top_families(
    cis_dist,
    "HERV Family Distribution Among Cisplatin-only Significant Loci",
    "Figure_Cisplatin_only_family_distribution.png"
)

# Carb-only plot
plot_top_families(
    carb_dist,
    "HERV Family Distribution Among Carboplatin-only Significant Loci",
    "Figure_Carboplatin_only_family_distribution.png"
)

# Optional: Shared plot
plot_top_families(
    shared_dist,
    "HERV Family Distribution Among Shared Significant Loci",
    "Figure_Shared_family_distribution.png"
)

# -----------------------------
# Enrichment testing (Fisher)
# -----------------------------
def fisher_enrichment(target_mask, label: str) -> pd.DataFrame:
    """
    For each family, test enrichment in target group vs all other loci.
    2x2 table per family:
        in_target_family     in_target_not_family
        out_target_family    out_target_not_family
    """
    results = []
    families = sorted(df["Family"].dropna().unique())

    for fam in families:
        in_target_family = ((df["Family"] == fam) & target_mask).sum()
        in_target_notfam = ((df["Family"] != fam) & target_mask).sum()
        out_target_family = ((df["Family"] == fam) & (~target_mask)).sum()
        out_target_notfam = ((df["Family"] != fam) & (~target_mask)).sum()

        table = [
            [in_target_family, in_target_notfam],
            [out_target_family, out_target_notfam],
        ]

        # Fisher's exact test (enrichment): alternative='greater'
        odds, p = fisher_exact(table, alternative="greater")

        results.append({
            "Family": fam,
            "InTarget_Family": in_target_family,
            "InTarget_Total": int(target_mask.sum()),
            "OutTarget_Family": out_target_family,
            "OutTarget_Total": int((~target_mask).sum()),
            "OddsRatio": odds,
            "PValue": p
        })

    res = pd.DataFrame(results)

    # Multiple testing correction (BH FDR)
    res["FDR"] = multipletests(res["PValue"], method="fdr_bh")[1]
    res["GroupTested"] = label

    # Sort by FDR then p
    res = res.sort_values(["FDR", "PValue"], ascending=True).reset_index(drop=True)
    return res

cis_mask = df["Significance_Category"] == "Cisplatin-only significant"
carb_mask = df["Significance_Category"] == "Carboplatin-only significant"

cis_enrich = fisher_enrichment(cis_mask, "Cisplatin-only vs Others")
carb_enrich = fisher_enrichment(carb_mask, "Carboplatin-only vs Others")

cis_enrich.to_csv("Cisplatin_only_family_enrichment_fisher.tsv", sep="\t", index=False)
carb_enrich.to_csv("Carboplatin_only_family_enrichment_fisher.tsv", sep="\t", index=False)

print("\nSaved enrichment tables:")
print(" - Cisplatin_only_family_enrichment_fisher.tsv")
print(" - Carboplatin_only_family_enrichment_fisher.tsv")

# Print top enriched families (FDR < 0.05) if any
print("\nTop Cisplatin-only enriched families (FDR < 0.05):")
print(cis_enrich[cis_enrich["FDR"] < 0.05].head(10).to_string(index=False))

print("\nTop Carboplatin-only enriched families (FDR < 0.05):")
print(carb_enrich[carb_enrich["FDR"] < 0.05].head(10).to_string(index=False))
