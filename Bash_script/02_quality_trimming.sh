#!/bin/bash
#
# quality_trimming.sh
#
# Performs quality and length trimming of paired-end FASTQ reads using
# Trimmomatic (paired-end mode). Reads are trimmed for adapter-region
# artifacts (fixed-length head crop), low-quality leading/trailing bases,
# a sliding-window quality filter, a minimum length, and a minimum average
# quality.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Usage:
#   ./quality_trimming.sh -i INPUT_DIR -o OUTPUT_DIR
#
# Requirements: Trimmomatic and Java in PATH
#
# Note: trimming parameters (HEADCROP:19, SLIDINGWINDOW:4:15, MINLEN:100,
# AVGQUAL:30, etc.) are unchanged from the original script and reflect the
# settings used for this study's sequencing run/read length.

set -euo pipefail

usage() {
  echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR"
  exit 1
}

INPUT_DIR=""
OUTPUT_DIR=""

while getopts ":i:o:" opt; do
  case "$opt" in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    *) usage ;;
  esac
done

[ -z "${INPUT_DIR}" ] && usage
[ -z "${OUTPUT_DIR}" ] && usage

[ -d "${INPUT_DIR}" ] || { echo "Error: input directory not found: ${INPUT_DIR}" >&2; exit 1; }
mkdir -p "${OUTPUT_DIR}"

PAIRED_DIR="${OUTPUT_DIR%/}/paired"
UNPAIRED_DIR="${OUTPUT_DIR%/}/unpaired"
LOG_DIR="${OUTPUT_DIR%/}/logs"
mkdir -p "${PAIRED_DIR}" "${UNPAIRED_DIR}" "${LOG_DIR}"

command -v trimmomatic >/dev/null 2>&1 || { echo "Error: trimmomatic not found in PATH" >&2; exit 1; }
command -v java >/dev/null 2>&1 || { echo "Error: java not found in PATH" >&2; exit 1; }

mapfile -d '' R1_FILES < <(find "${INPUT_DIR}" -type f -name "*_R1_001.fastq.gz" -print0)

if [ "${#R1_FILES[@]}" -eq 0 ]; then
  echo "Error: no *_R1_001.fastq.gz files found in ${INPUT_DIR}" >&2
  exit 1
fi

echo "Found ${#R1_FILES[@]} sample(s) (R1)."

for R1 in "${R1_FILES[@]}"; do
  base="$(basename "$R1")"
  sample_prefix="${base%_R1_001.fastq.gz}"

  R2="${INPUT_DIR%/}/${sample_prefix}_R2_001.fastq.gz"
  if [ ! -f "$R2" ]; then
    echo "Warning: matching R2 file not found for ${sample_prefix}; sample skipped." >&2
    continue
  fi

  out_R1_paired="${PAIRED_DIR}/${sample_prefix}_R1_001.fastq.gz"
  out_R1_unpaired="${UNPAIRED_DIR}/${sample_prefix}_R1_unpaired.fastq.gz"
  out_R2_paired="${PAIRED_DIR}/${sample_prefix}_R2_001.fastq.gz"
  out_R2_unpaired="${UNPAIRED_DIR}/${sample_prefix}_R2_unpaired.fastq.gz"
  log_file="${LOG_DIR}/${sample_prefix}.trim.log"

  echo "Processing: ${sample_prefix}"
  set +e
  trimmomatic PE -threads "$(nproc)" -phred33 \
    "$R1" "$R2" \
    "$out_R1_paired" "$out_R1_unpaired" \
    "$out_R2_paired" "$out_R2_unpaired" \
    HEADCROP:19 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 AVGQUAL:30 \
    2> "$log_file"
  status=$?
  set -e

  if [ $status -ne 0 ]; then
    echo "Error: Trimmomatic failed for ${sample_prefix}. See ${log_file}" >&2
    continue
  fi

  echo "OK: ${sample_prefix}"
done

echo "Quality trimming completed. Output:"
echo "  Paired:   ${PAIRED_DIR}"
echo "  Unpaired: ${UNPAIRED_DIR}"
echo "  Logs:     ${LOG_DIR}"
