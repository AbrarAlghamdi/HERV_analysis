#!/usr/bin/env bash
set -euo pipefail

OUT="melanocyte_aug_summary.tsv"

echo -e "sample\tbam_size_GB\tmapped_reads\tproperly_paired_reads\therv_loci\ttotal_final_count\ttsv" > "$OUT"

# loop over all Results directories
find . -type d -name Results | while read R; do
    SAMPLE=$(basename "$(dirname "$R")")

    BAM="$R/telescope-updated.bam"
    TSV=$(ls "$R"/Aug_*.tsv 2>/dev/null | head -n 1)

    # skip if missing
    [[ -f "$BAM" && -f "$TSV" ]] || continue

    # BAM size in GB
    SIZE_GB=$(du -b "$BAM" | awk '{printf "%.3f", $1/1024/1024/1024}')

    # alignment stats
    MAPPED=$(samtools view -c -F 4 "$BAM")
    PPAIRED=$(samtools view -c -f 2 "$BAM")

    # HERV stats
    HERV_LOCI=$(awk 'NR>1{c++} END{print c+0}' "$TSV")
    TOTAL_FINAL=$(awk 'NR>1{s+=$3} END{printf "%.0f", s}' "$TSV")

    echo -e "${SAMPLE}\t${SIZE_GB}\t${MAPPED}\t${PPAIRED}\t${HERV_LOCI}\t${TOTAL_FINAL}\t$(basename "$TSV")" >> "$OUT"
done

echo "Saved $OUT"



A. Properly paired %
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{print $0,"\tproperly_paired_percent"}
NR>1{printf "%s\t%.1f\n",$0,($4/$3)*100}' melanocyte_aug_summary.tsv \
> melanocyte_aug_summary_with_percent.tsv



B. % of total HERV counts
TOTAL=$(awk -F'\t' 'NR>1{s+=$6} END{print s}' melanocyte_aug_summary.tsv)

awk -F'\t' -v T=$TOTAL 'BEGIN{OFS="\t"}
NR==1{print $0,"\therv_fraction_percent"}
NR>1{printf "%s\t%.1f\n",$0,($6/T)*100}' melanocyte_aug_summary_with_percent.tsv \
> melanocyte_aug_summary_final.tsv
