#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

OUT="skmel28_01_summary.tsv"
echo -e "sample\tfastq_R1_reads\tfastq_R2_reads\taligned_pairs_total\taligned_pairs_mapped\tmapped_percent\tproperly_paired_percent\tte_bam_pairs\therv_loci\ttotal_herv_counts" > "$OUT"

for d in B*_r_*_s_*; do
  [[ -d "$d" ]] || continue
  sample="$(basename "$d")"
  R="$d/Results"

  # -------------------------
  # 1) FASTQ read counts
  # -------------------------
  R1_LINES=$(zcat "$d"/*_R1_*.fastq.gz 2>/dev/null | wc -l || true)
  R2_LINES=$(zcat "$d"/*_R2_*.fastq.gz 2>/dev/null | wc -l || true)
  R1_READS=$((R1_LINES/4))
  R2_READS=$((R2_LINES/4))

  # ---------------------------------------------
  # 2) Alignment stats from per-lane sorted BAMs
  # ---------------------------------------------
  # Prefer per-lane *sorted.bam; sum across lanes.
  total_pairs=0
  mapped_pairs=0
  pp_pairs=0

  bams=("$R"/*_L???.sorted.bam "$R"/*_L???.sorted.bam 2>/dev/null || true)
  # If glob fails silently, just use ls:
  if ! ls "$R"/*_L???.sorted.bam >/dev/null 2>&1; then
    bams=()
  else
    mapfile -t bams < <(ls "$R"/*_L???.sorted.bam | sort)
  fi

  if [[ ${#bams[@]} -gt 0 ]]; then
    for bam in "${bams[@]}"; do
      # paired in sequencing = number of read pairs (counts alignments, but stable for QC)
      p=$(samtools flagstat "$bam" | awk '/paired in sequencing/ {print $1; exit}')
      m=$(samtools flagstat "$bam" | awk '/ mapped \(/ && !/primary/ {print $1; exit}')
      pp=$(samtools flagstat "$bam" | awk '/properly paired/ {print $1; exit}')

      # totals are integers
      total_pairs=$((total_pairs + p))
      mapped_pairs=$((mapped_pairs + m))
      pp_pairs=$((pp_pairs + pp))
    done

    # percentages
    if [[ $total_pairs -gt 0 ]]; then
      mapped_pct=$(awk -v a="$mapped_pairs" -v b="$total_pairs" 'BEGIN{printf "%.2f", (a/b)*100}')
      pp_pct=$(awk -v a="$pp_pairs" -v b="$total_pairs" 'BEGIN{printf "%.2f", (a/b)*100}')
    else
      mapped_pct="NA"
      pp_pct="NA"
    fi
  else
    total_pairs="NA"
    mapped_pairs="NA"
    mapped_pct="NA"
    pp_pct="NA"
  fi

  # -------------------------
  # 3) Telescope BAM stats
  # -------------------------
  TE_BAM="$R/telescope-updated.bam"
  TE_PAIRS="NA"
  if [[ -f "$TE_BAM" ]]; then
    TE_PAIRS=$(samtools flagstat "$TE_BAM" | awk '/paired in sequencing/ {print $1; exit}')
  fi

  # -------------------------
  # 4) Telescope TSV stats
  # -------------------------
  KEY=$(echo "$sample" | sed -E 's/(B[0-9]+)_r_([0-9]+)_s_([0-9]+)/\1_R\2_S\3/' | tr '[:lower:]' '[:upper:]')

  TSV=""
  if [[ -f "$R/telescope-telescope_report.tsv" ]]; then
    TSV="$R/telescope-telescope_report.tsv"
  else
    TSV=$(ls "$R"/*.tsv 2>/dev/null | grep -i "/${KEY}\.tsv$" | head -n 1 || true)
    [[ -z "${TSV:-}" ]] && TSV=$(ls "$R"/*.tsv 2>/dev/null | grep -vi "checkpoint" | head -n 1 || true)
  fi

  HERV_LOCI="NA"
  TOTAL_HERV="NA"
  if [[ -n "${TSV:-}" && -f "$TSV" ]]; then
    HERV_LOCI=$(awk 'NR>1{c++} END{print c+0}' "$TSV")
    TOTAL_HERV=$(awk 'NR>1{s+=$3} END{printf "%.0f", s+0}' "$TSV")
  fi

  echo -e "${sample}\t${R1_READS}\t${R2_READS}\t${total_pairs}\t${mapped_pairs}\t${mapped_pct}\t${pp_pct}\t${TE_PAIRS}\t${HERV_LOCI}\t${TOTAL_HERV}" >> "$OUT"
done

echo "Saved QC table: $OUT"
echo "Preview:"
column -t "$OUT" | head
