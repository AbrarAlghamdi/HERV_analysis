import pandas as pd

############################################
# Load Differential Expression Results
############################################

skmel5 = pd.read_csv(
    "differential_expression_significant_fdr_skmel_5.tsv",
    sep="\t"
)

skmel28 = pd.read_csv(
    "differential_expression_ttest_fdr_skmel28.tsv",
    sep="\t"
)

# Keep only significant
skmel5 = skmel5[skmel5["significant_fdr05"] == True]
skmel28 = skmel28[skmel28["significant_fdr05"] == True]

############################################
#  Parse GTF Annotation
############################################

gtf_file = "HERV_rmsk.hg38.v2.2.gtf"

locus_data = {}

with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        chrom = fields[0]
        feature = fields[2]
        attributes = fields[8]

        # Extract gene_id
        gene_id = None
        for attr in attributes.split(";"):
            if "gene_id" in attr:
                gene_id = attr.split('"')[1]

        if gene_id is None:
            continue

        if gene_id not in locus_data:
            locus_data[gene_id] = {
                "chromosome": chrom,
                "repFamily": set(),
                "has_internal": False,
                "has_ltr": False
            }

        if "repFamily" in attributes:
            for attr in attributes.split(";"):
                if "repFamily" in attr:
                    fam = attr.split('"')[1]
                    locus_data[gene_id]["repFamily"].add(fam)

        if 'geneRegion "internal"' in attributes:
            locus_data[gene_id]["has_internal"] = True

        if 'geneRegion "ltr"' in attributes:
            locus_data[gene_id]["has_ltr"] = True

############################################
#  Convert to DataFrame
############################################

gtf_df = []

for locus, info in locus_data.items():
    if info["has_internal"] and info["has_ltr"]:
        structure = "provirus_like"
    elif info["has_ltr"] and not info["has_internal"]:
        structure = "solo_LTR"
    else:
        structure = "fragmented"

    gtf_df.append({
        "transcript": locus,
        "chromosome": info["chromosome"],
        "repFamily": ",".join(info["repFamily"]),
        "structure": structure
    })

gtf_df = pd.DataFrame(gtf_df)

############################################
#  Merge With DE Lists
############################################

skmel5_annot = skmel5.merge(gtf_df, on="transcript", how="left")
skmel28_annot = skmel28.merge(gtf_df, on="transcript", how="left")

############################################
#  Chromosome Distribution
############################################

print("\nSK-MEL-5 Chromosome Distribution:")
print(skmel5_annot["chromosome"].value_counts())

print("\nSK-MEL-28 Chromosome Distribution:")
print(skmel28_annot["chromosome"].value_counts())

############################################
# Family Distribution
############################################

print("\nSK-MEL-5 repFamily Distribution:")
print(skmel5_annot["repFamily"].value_counts())

print("\nSK-MEL-28 repFamily Distribution:")
print(skmel28_annot["repFamily"].value_counts())

############################################
# Structural Pattern
############################################

print("\nSK-MEL-5 Structure Distribution:")
print(skmel5_annot["structure"].value_counts())

print("\nSK-MEL-28 Structure Distribution:")
print(skmel28_annot["structure"].value_counts())
