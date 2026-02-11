import pandas as pd
import matplotlib.pyplot as plt

# Load files
skmel5 = pd.read_csv("differential_expression_significant_fdr_skmel_5.tsv", sep="\t")
skmel28 = pd.read_csv("differential_expression_ttest_fdr_skmel28.tsv", sep="\t")

# Keep only significant rows
skmel5 = skmel5[skmel5["significant_fdr05"] == True]
skmel28 = skmel28[skmel28["significant_fdr05"] == True]

def extract_chr(name):
    try:
        part = name.split("_")[1]
        chr_part = part.split("p")[0].split("q")[0]
        return "chr" + chr_part
    except:
        return "unknown"

skmel5["chromosome"] = skmel5["transcript"].apply(extract_chr)
skmel28["chromosome"] = skmel28["transcript"].apply(extract_chr)

# Count chromosomes
chr5_counts = skmel5["chromosome"].value_counts().sort_index()
chr28_counts = skmel28["chromosome"].value_counts().sort_index()

print("SK-MEL-5 chromosome distribution:")
print(chr5_counts)

print("\nSK-MEL-28 chromosome distribution:")
print(chr28_counts)

# Plot SK-MEL-5
plt.figure()
chr5_counts.plot(kind="bar")
plt.title("Chromosome Distribution - SK-MEL-5")
plt.ylabel("Number of significant HERV loci")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

def extract_family(name):
    return name.split("_")[0]

skmel5["family"] = skmel5["transcript"].apply(extract_family)
skmel28["family"] = skmel28["transcript"].apply(extract_family)

fam5_counts = skmel5["family"].value_counts()
fam28_counts = skmel28["family"].value_counts()

print("\nSK-MEL-5 family distribution:")
print(fam5_counts)

print("\nSK-MEL-28 family distribution:")
print(fam28_counts)

plt.figure()
fam5_counts.plot(kind="bar")
plt.title("HERV Family Distribution - SK-MEL-5")
plt.ylabel("Count")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
