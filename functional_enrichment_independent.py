#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Inputs
# -----------------------------
CARB_FILE = "/home/a1321/Desktop/Bladder/carp/indepent/HERV_gene_pairs_classified.tsv"
CIS_FILE  = "/home/a1321/Desktop/Bladder/Cis/conf/HERV_gene_pairs_classified.tsv"

# Enrichment libraries
GENE_SETS = [
    "GO_Biological_Process_2021",
    "GO_Molecular_Function_2021",
    "KEGG_2021_Human"
]

# Significance / reporting
FDR_CUTOFF = 0.05         # strict significance
RELAX_FDR_FOR_PLOT = 0.20 # optional for plotting when strict has none
TOP_TERMS = 15
DPI = 300

# -----------------------------
# gseapy
# -----------------------------
try:
    import gseapy as gp
except ImportError:
    raise SystemExit("gseapy not installed. Run: pip install gseapy")

# -----------------------------
# Helpers
# -----------------------------
def load_independent_genes(path: str):
    df = pd.read_csv(path, sep="\t")
    genes = (df[df["Regulatory_Pattern"] == "Independent"]["gene_name"]
             .dropna()
             .astype(str)
             .unique()
             .tolist())
    return genes

def run_enrichr_safe(gene_list, label):
    """
    Run enrichr but do NOT crash if nothing passes cutoff.
    Set cutoff=1.0 so gseapy returns results; we filter later ourselves.
    Also set no_plot=True to avoid gseapy internal plotting crash.
    """
    outdir = f"Independent_Enrichment_{label}"
    os.makedirs(outdir, exist_ok=True)

    enr = gp.enrichr(
        gene_list=gene_list,
        organism="Human",
        gene_sets=GENE_SETS,
        outdir=outdir,
        cutoff=1.0,     # return everything; we'll filter ourselves
        no_plot=True    # critical: avoids ValueError when nothing significant
    )
    res = enr.results.copy()
    res["Dataset"] = label
    return res, outdir

def plot_top_terms(res, label, outpng, fdr_for_plot=RELAX_FDR_FOR_PLOT):
    """
    Plot top terms using a relaxed FDR threshold if needed.
    If still empty, skip plot gracefully.
    """
    if res is None or res.empty:
        print(f"[{label}] No results returned. Skipping plot.")
        return

    # Choose terms for plotting
    res_plot = res.sort_values("Adjusted P-value").copy()

    # Try strict first
    strict = res_plot[res_plot["Adjusted P-value"] < FDR_CUTOFF]
    if not strict.empty:
        top = strict.head(TOP_TERMS)
        used = f"FDR<{FDR_CUTOFF}"
    else:
        relaxed = res_plot[res_plot["Adjusted P-value"] < fdr_for_plot]
        if not relaxed.empty:
            top = relaxed.head(TOP_TERMS)
            used = f"FDR<{fdr_for_plot} (relaxed for visualization)"
        else:
            # As last resort: take top nominal p-values (still informative)
            top = res_plot.head(TOP_TERMS)
            used = "top nominal terms (no FDR-significant terms)"

    # Plot
    plt.figure(figsize=(9, 6))
    plt.barh(top["Term"], -np.log10(top["Adjusted P-value"].astype(float)))
    plt.xlabel("-log10 Adjusted P-value")
    plt.ylabel("Enriched term")
    plt.title(f"Independent-HERV Gene Enrichment ({label})\n({used})")
    plt.tight_layout()
    plt.savefig(outpng, dpi=DPI)
    plt.show()
    print(f"Saved: {outpng}")

# -----------------------------
# Main
# -----------------------------
print("Loading Independent genes...")

carb_genes = load_independent_genes(CARB_FILE)
cis_genes  = load_independent_genes(CIS_FILE)

print(f"Carboplatin Independent genes: {len(carb_genes)}")
print(f"Cisplatin Independent genes:  {len(cis_genes)}")

pd.Series(carb_genes).to_csv("Independent_genes_Carboplatin.txt", index=False, header=False)
pd.Series(cis_genes).to_csv("Independent_genes_Cisplatin.txt", index=False, header=False)
print("Saved gene lists.")

print("\nRunning Enrichr enrichment (safe mode, no_plot=True)...")
carb_res, _ = run_enrichr_safe(carb_genes, "Carboplatin")
cis_res,  _ = run_enrichr_safe(cis_genes,  "Cisplatin")

carb_res.to_csv("Independent_enrichment_Carboplatin.tsv", sep="\t", index=False)
cis_res.to_csv("Independent_enrichment_Cisplatin.tsv", sep="\t", index=False)
print("Saved enrichment tables.")

# Strict significant hits (for reporting)
carb_sig = carb_res[carb_res["Adjusted P-value"] < FDR_CUTOFF]
cis_sig  = cis_res[cis_res["Adjusted P-value"] < FDR_CUTOFF]

print("\nSignificant enriched terms (FDR < 0.05):")
print("  Carboplatin:", len(carb_sig))
print("  Cisplatin: ", len(cis_sig))

# Plot (will gracefully relax if none)
plot_top_terms(carb_res, "Carboplatin", "Independent_enrichment_Carboplatin_barplot.png")
plot_top_terms(cis_res,  "Cisplatin",  "Independent_enrichment_Cisplatin_barplot.png")

# Combined table for convenience
combined = pd.concat([carb_res, cis_res], ignore_index=True)
combined.to_csv("Independent_enrichment_combined.tsv", sep="\t", index=False)
print("\nSaved combined table: Independent_enrichment_combined.tsv")
print("Done.")
