A) Build BED for HERV loci from HERV_rmsk.hg38.v2.gtf

(collapses exons → one interval per transcript)

awk -F'\t' 'BEGIN{OFS="\t"} $0!~/^#/ && $3=="exon" {
  tid="";
  if (match($9, /transcript_id "([^"]+)"/, a)) tid=a[1];
  if (tid=="") next;
  chr=$1; start=$4; end=$5; strand=$7;
  key=chr"|"tid"|"strand;
  if (!(key in min) || start<min[key]) min[key]=start;
  if (!(key in max) || end>max[key])   max[key]=end;
}
END{
  for (k in min){
    split(k,b,"|");
    chr=b[1]; tid=b[2]; strand=b[3];
    print chr, min[k]-1, max[k], tid, 0, strand;   # BED 0-based start
  }
}' HERV_rmsk.hg38.v2.gtf | sort -k1,1 -k2,2n > herv_loci.bed

B) Build BED for protein-coding genes from GENCODE v44

(using gene features)

awk -F'\t' 'BEGIN{OFS="\t"} $0!~/^#/ && $3=="gene" {
  gtype=""; gname=""; gid="";
  if (match($9, /gene_type "([^"]+)"/, a)) gtype=a[1];
  if (gtype!="protein_coding") next;
  if (match($9, /gene_name "([^"]+)"/, b)) gname=b[1];
  if (match($9, /gene_id "([^"]+)"/, c)) gid=c[1];
  print $1, $4-1, $5, gname, gid, $7;
}' gencode.v44.annotation.gtf | sort -k1,1 -k2,2n > gencode_pc_genes.bed

C) Extract significant HERV IDs (padj < 0.05) and subset BED
python3 - <<'PY'
import pandas as pd
df = pd.read_csv("DESeq2_HERV_Carboplatin_vs_Untreated.csv")
sig = df[df["padj"] < 0.05]["transcript"].dropna().astype(str)
sig.to_csv("sig_herv_ids.txt", index=False, header=False)
print("Significant HERV loci:", len(sig))
PY

grep -F -f sig_herv_ids.txt herv_loci.bed > sig_herv_loci.bed

D) Map each significant HERV to nearest protein-coding gene + keep ≤10kb
bedtools closest -a sig_herv_loci.bed -b gencode_pc_genes.bed -d > herv_nearest_gene.tsv
awk 'BEGIN{OFS="\t"} $NF<=10000 {print}' herv_nearest_gene.tsv > herv_nearest_gene_10kb.tsv
