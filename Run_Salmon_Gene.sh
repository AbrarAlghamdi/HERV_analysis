#!/usr/bin/env bash
set -euo pipefail

# -------- paths --------
BASE="/DataDrive/4TB/bladder_data/Carboplatin"
TREATED_DIR="${BASE}/treated"
UNTREATED_DIR="${BASE}/Untreated"

# Salmon index (transcriptome index)
SALMON_INDEX="/path/to/salmon_index"   # 

# Output folder
OUTDIR="${BASE}/salmon_gene_quant_all"
mkdir -p "${OUTDIR}"

# -------- function to quantify 1 sample --------
run_one () {
  local fq="$1"
  local sample="$2"
  local out="${OUTDIR}/${sample}"

  echo "==> Quantifying: ${sample}"
  salmon quant \
    -i "${SALMON_INDEX}" \
    -l A \
    -r "${fq}" \
    -p 12 \
    --validateMappings \
    -o "${out}"

  # keep a convenient copy of quant.sf in one place (optional)
  cp "${out}/quant.sf" "${OUTDIR}/${sample}.quant.sf"
}

# -------- run treated --------
for fq in "${TREATED_DIR}"/*.fastq.gz; do
  sample="$(basename "$fq" .fastq.gz)"      # e.g. SRR24938919_1
  sample="${sample%_1}"                      # -> SRR24938919
  sample="Treated_${sample}"                 # -> Treated_SRR...
  run_one "$fq" "$sample"
done

# -------- run untreated --------
for fq in "${UNTREATED_DIR}"/*.fastq.gz; do
  sample="$(basename "$fq" .fastq.gz)"
  sample="${sample%_1}"
  sample="Untreated_${sample}"
  run_one "$fq" "$sample"
done

echo "Done. Per-sample outputs are in: ${OUTDIR}/<sample>/quant.sf"
echo "And copies are in: ${OUTDIR}/*.quant.sf"
