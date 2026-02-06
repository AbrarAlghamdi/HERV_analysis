#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

######################## CONFIG (overridable) ########################
ROOT="${ROOT:-/DataDrive/4TB/Skmel5_data/RNAseq}"    # expects subfolders 1/ and 10/
THREADS="${THREADS:-16}"

# Bowtie2 index **basename** (NO .fa). Use one that exists on your box:
#   /home/aasa13/Desktop/Telescope/refs/hg38_bt2/hg38
#   /home/aasa13/Desktop/Telescope/refs/hg38.fasta
#   /home/aasa13/Desktop/Telescope/GTF/hg38
BOWTIE2_INDEX="${BOWTIE2_INDEX:-/home/aasa13/Desktop/Telescope/refs/hg38_bt2/hg38}"

# HERV GTF
GTF_FILE="${GTF_FILE:-/home/aasa13/Desktop/Telescope/refs/HERV_rmsk.hg38.v2.2.gtf}"

# Conda envs
ALIGN_ENV="${ALIGN_ENV:-singleloc_align}"
TELE_ENV="${TELE_ENV:-singleloc_tele}"
######################################################################

# Sanity checks
[[ -d "$ROOT/1" && -d "$ROOT/10" ]] || { echo "ERR: expected $ROOT/{1,10}"; exit 1; }
[[ -f "$GTF_FILE" ]] || { echo "ERR: GTF not found: $GTF_FILE"; exit 1; }
if ! ls "${BOWTIE2_INDEX}".1.bt2 >/dev/null 2>&1 && ! ls "${BOWTIE2_INDEX}".1.bt2l >/dev/null 2>&1; then
  echo "ERR: Bowtie2 index not found for basename: ${BOWTIE2_INDEX}"
  echo "     Expected files like ${BOWTIE2_INDEX}.1.bt2 (or .1.bt2l)"
  exit 1
fi
command -v conda >/dev/null || { echo "ERR: conda not found in PATH"; exit 1; }
CONDA_BASE=$(conda info --base); source "$CONDA_BASE/etc/profile.d/conda.sh"

echo "== Run all (no merge) =="
echo "ROOT         : $ROOT"
echo "BOWTIE2_INDEX: $BOWTIE2_INDEX"
echo "GTF_FILE     : $GTF_FILE"
echo "THREADS      : $THREADS"
echo

process_folder () {
  local FASTQ_DIR="$1"     # e.g., /.../RNAseq/10/B_1_R_1_s_2
  [[ -d "$FASTQ_DIR" ]] || return 0

  local SAMPLE_TAG; SAMPLE_TAG="$(basename "$FASTQ_DIR")"
  local OUTDIR="${FASTQ_DIR}/Results"
  mkdir -p "$OUTDIR"

  # 1) lanes present in BOTH R1 and R2
  mapfile -t lanes_R1 < <(ls "$FASTQ_DIR"/*_L???_R1_*.fastq.gz 2>/dev/null | sed -E 's/.*_L([0-9]{3})_.*/\1/' | sort -u)
  mapfile -t lanes_R2 < <(ls "$FASTQ_DIR"/*_L???_R2_*.fastq.gz 2>/dev/null | sed -E 's/.*_L([0-9]{3})_.*/\1/' | sort -u)
  if [[ ${#lanes_R1[@]} -eq 0 || ${#lanes_R2[@]} -eq 0 ]]; then
    echo "[SKIP] $FASTQ_DIR (missing R1 or R2 lanes)"; return 0
  fi
  local common=()
  for l in "${lanes_R1[@]}"; do
    if printf '%s\n' "${lanes_R2[@]}" | grep -qx "$l"; then common+=("$l"); fi
  done
  if [[ ${#common[@]} -eq 0 ]]; then
    echo "[SKIP] $FASTQ_DIR (no common lanes)"; return 0
  fi
  [[ ${#common[@]} -ne ${#lanes_R1[@]} || ${#common[@]} -ne ${#lanes_R2[@]} ]] && \
    echo "[NOTE] $SAMPLE_TAG: using common lanes only: ${common[*]}"

  # 2) build comma-separated lists for bowtie2 (NO merging)
  local R1_LIST=() R2_LIST=()
  for L in "${common[@]}"; do
    R1_LIST+=("$FASTQ_DIR"/*_L${L}_R1_*.fastq.gz)
    R2_LIST+=("$FASTQ_DIR"/*_L${L}_R2_*.fastq.gz)
  done
  local R1_CSV; R1_CSV=$(IFS=, ; echo "${R1_LIST[*]}")
  local R2_CSV; R2_CSV=$(IFS=, ; echo "${R2_LIST[*]}")

  # 3) align → SAM → BAM → name-sort (collate)
  local SAM="$OUTDIR/${SAMPLE_TAG}.sam"
  local BAM="$OUTDIR/${SAMPLE_TAG}.bam"
  local NSORT="$OUTDIR/${SAMPLE_TAG}.namesort.bam"

  echo "[Bowtie2] $SAMPLE_TAG"
  conda activate "$ALIGN_ENV"
  bowtie2 \
    -k 50 \
    --very-sensitive-local \
    --score-min "L,0,2" \
    --rg-id "$SAMPLE_TAG" \
    -x "$BOWTIE2_INDEX" \
    -1 "$R1_CSV" \
    -2 "$R2_CSV" \
    -S "$SAM" \
    -p "$THREADS" \
    2>&1 | tee "$OUTDIR/${SAMPLE_TAG}.bowtie2.log"

  samtools view -@ "$THREADS" -Sb "$SAM" > "$BAM"
  samtools collate -@ "$THREADS" "$BAM" -o "$NSORT"
  conda deactivate

  # 4) Telescope
  echo "[Telescope] $SAMPLE_TAG"
  local TE_SAM="$OUTDIR/${SAMPLE_TAG}-telescope-updated.sam"
  conda activate "$TELE_ENV"
  telescope assign \
    --theta_prior 200000 \
    --max_iter 200 \
    --updated_sam \
    --outdir "$OUTDIR" \
    --attribute gene_id \
    --exp_tag telescope \
    "$NSORT" \
    "$GTF_FILE" \
    2>&1 | tee "$OUTDIR/${SAMPLE_TAG}.telescope.log"

  if [[ -f "$OUTDIR/telescope-updated.bam" ]]; then
    samtools view -h -o "$TE_SAM" "$OUTDIR/telescope-updated.bam"
    mv -f "$OUTDIR/telescope-updated.bam" "$OUTDIR/${SAMPLE_TAG}-telescope-updated.bam"
  else
    echo "[WARN] telescope-updated.bam not found in $OUTDIR"
  fi
  conda deactivate

  echo "[DONE] $SAMPLE_TAG → $OUTDIR"; echo
}

# Run for all folders in 1/ and 10/
for cond in 1 10; do
  for d in "$ROOT/$cond"/B_*_R_*_s_*; do
    [[ -d "$d" ]] || continue
    echo "=== Processing: $d ==="
    process_folder "$d"
  done
done

echo "All finished."
