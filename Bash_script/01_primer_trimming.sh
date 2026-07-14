#!/bin/bash
#
# primer_trimming.sh
#
# Removes Illumina sequencing adapters and locus-specific PCR primers from
# paired-end amplicon sequencing reads using Cutadapt. Supports the three
# marker genes used in this study: 16S rRNA (bacteria/archaea), 18S rRNA
# (eukaryotes), and the fungal ITS region.
#
# This script consolidates what were originally three near-identical
# scripts (one per marker) into a single, parameterized script. The
# Cutadapt command and the adapter/primer sequences for each marker are
# unchanged from the original per-marker scripts.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Usage:
#   ./primer_trimming.sh -i INPUT_DIR -o OUTPUT_DIR -m MARKER
#
#   -i  Input directory containing paired-end FASTQ files
#       (*_R1_001.fastq.gz / *_R2_001.fastq.gz)
#   -o  Output directory for primer-trimmed FASTQ files
#   -m  Marker gene: 16S, 18S, or ITS
#
# Requirements: Cutadapt in PATH

set -uo pipefail

usage() {
  echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -m {16S|18S|ITS}"
  exit 1
}

INPUT_DIR=""
OUTPUT_DIR=""
MARKER=""

while getopts ":i:o:m:" opt; do
  case "$opt" in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    m) MARKER="$OPTARG" ;;
    *) usage ;;
  esac
done

[ -z "$INPUT_DIR" ] && usage
[ -z "$OUTPUT_DIR" ] && usage
[ -z "$MARKER" ] && usage

command -v cutadapt >/dev/null 2>&1 || { echo "Error: cutadapt not found in PATH" >&2; exit 1; }
[ -d "$INPUT_DIR" ] || { echo "Error: input directory not found: $INPUT_DIR" >&2; exit 1; }

# Adapter and primer sequences per marker gene.
# These are identical to the sequences used in the original per-marker scripts.
case "$MARKER" in
  16S)
    ADAPTER_FWD="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    ADAPTER_REV="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    PRIMER_FWD="CCTACGGGNGGCWGCAG"
    PRIMER_REV="GACTACHVGGGTATCTAATCC"
    ;;
  18S)
    ADAPTER_FWD="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    ADAPTER_REV="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    PRIMER_FWD="GCCGCGGTAATTCCAGCTC"
    PRIMER_REV="CYTTCGYYCTTGATTRA"
    ;;
  ITS)
    ADAPTER_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    ADAPTER_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
    PRIMER_FWD="GTGAATTCTCGAGAATTC"
    PRIMER_REV="CTGTCTCTTATACACATCT"
    ;;
  *)
    echo "Error: unknown marker '$MARKER'. Must be one of: 16S, 18S, ITS" >&2
    exit 1
    ;;
esac

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob
R1_FILES=("$INPUT_DIR"/*_R1_001.fastq.gz)
shopt -u nullglob

if [ "${#R1_FILES[@]}" -eq 0 ]; then
  echo "Error: no *_R1_001.fastq.gz files found in $INPUT_DIR" >&2
  exit 1
fi

echo "Marker: $MARKER"
echo "Found ${#R1_FILES[@]} R1 file(s) in $INPUT_DIR"

for R1 in "${R1_FILES[@]}"; do
  R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

  if [ -f "$R2" ]; then
    R1_trimmed="$OUTPUT_DIR/$(basename "$R1")"
    R2_trimmed="$OUTPUT_DIR/$(basename "$R2")"

    cutadapt -a "$ADAPTER_FWD" -A "$ADAPTER_REV" \
             -g "$PRIMER_FWD" -G "$PRIMER_REV" \
             -o "$R1_trimmed" -p "$R2_trimmed" \
             "$R1" "$R2"

    echo "Processed: $(basename "$R1") / $(basename "$R2")"
  else
    echo "Warning: matching R2 file not found for $R1; sample skipped." >&2
  fi
done

echo "Primer trimming completed. Output: $OUTPUT_DIR"
