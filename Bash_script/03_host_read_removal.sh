#!/bin/bash
#
# host_read_removal.sh
#
# Removes host (human) reads from paired-end FASTQ files by aligning them
# against a human reference genome with Bowtie2 and retaining only the
# read pairs that do NOT align concordantly (i.e., non-host reads).
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Usage:
#   ./host_read_removal.sh -x HUMAN_BOWTIE2_INDEX -i INPUT_DIR -o OUTPUT_DIR [-t THREADS]
#
#   -x  Path/prefix of the Bowtie2 index built from the human reference genome
#   -i  Input directory containing paired-end FASTQ files
#       (*_R1_001.fastq.gz / *_R2_001.fastq.gz)
#   -o  Output directory for non-host (unaligned) FASTQ files
#   -t  Threads for Bowtie2 (default: 8)
#
# Requirements: Bowtie2 in PATH

set -uo pipefail

usage() {
  echo "Usage: $0 -x HUMAN_BOWTIE2_INDEX -i INPUT_DIR -o OUTPUT_DIR [-t THREADS]"
  exit 1
}

HUMAN_INDEX=""
INPUT_DIR=""
OUTPUT_DIR=""
THREADS=8

while getopts ":x:i:o:t:" opt; do
  case "$opt" in
    x) HUMAN_INDEX="$OPTARG" ;;
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    *) usage ;;
  esac
done

[ -z "$HUMAN_INDEX" ] && usage
[ -z "$INPUT_DIR" ] && usage
[ -z "$OUTPUT_DIR" ] && usage

command -v bowtie2 >/dev/null 2>&1 || { echo "Error: bowtie2 not found in PATH" >&2; exit 1; }
[ -d "$INPUT_DIR" ] || { echo "Error: input directory not found: $INPUT_DIR" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"

find "$INPUT_DIR" -maxdepth 1 -name "*_R1_001.fastq.gz" -print0 | while IFS= read -r -d $'\0' FILE1; do
  SAMPLE_NAME="${FILE1##*/}"
  SAMPLE_NAME="${SAMPLE_NAME%_R1_001.fastq.gz}"

  FILE2="${INPUT_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"

  if [ ! -f "$FILE2" ]; then
    echo "Warning: matching R2 file not found for ${SAMPLE_NAME}; sample skipped." >&2
    continue
  fi

  echo "Processing ${SAMPLE_NAME}..."

  # Align to the human reference genome. Read pairs that do not align
  # concordantly (non-host reads) are written to
  # ${OUTPUT_DIR}/${SAMPLE_NAME}.1 and .2 via --un-conc-gz.
  bowtie2 -p "$THREADS" -x "$HUMAN_INDEX" -1 "$FILE1" -2 "$FILE2" \
    --very-sensitive-local \
    --un-conc-gz "${OUTPUT_DIR}/${SAMPLE_NAME}" \
    > "${OUTPUT_DIR}/${SAMPLE_NAME}_mapped_and_unmapped.sam"

  mv "${OUTPUT_DIR}/${SAMPLE_NAME}.1" "${OUTPUT_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz"
  mv "${OUTPUT_DIR}/${SAMPLE_NAME}.2" "${OUTPUT_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"

  rm "${OUTPUT_DIR}/${SAMPLE_NAME}_mapped_and_unmapped.sam"

  echo "Non-host reads written: ${OUTPUT_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz, ${OUTPUT_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"
done

echo "Host read removal completed."
