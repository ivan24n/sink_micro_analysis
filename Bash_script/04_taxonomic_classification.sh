#!/bin/bash
#
# taxonomic_classification.sh
#
# Runs Kraken2 taxonomic classification on paired-end FASTQ files against a
# user-provided Kraken2 database, producing per-sample classification
# output and reports. Samples can optionally be processed in parallel
# using GNU parallel.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Usage:
#   ./taxonomic_classification.sh -i INPUT_DIR -o OUTPUT_DIR -d KRAKEN2_DB [-t THREADS] [-p PARALLEL_JOBS]
#
#   -i  Input directory containing paired-end FASTQ files
#       (*_R1_001.fastq.gz / *_R2_001.fastq.gz)
#   -o  Output directory
#   -d  Path to the Kraken2 database directory (contains hash.k2d, opts.k2d,
#       taxo.k2d, etc.)
#   -t  Threads per Kraken2 process (default: 4)
#   -p  Number of samples to process in parallel (default: 1; requires
#       GNU parallel if >1)
#
# Requirements: Kraken2 in PATH and an accessible Kraken2 database
set -euo pipefail

usage() {
  echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -d KRAKEN2_DB [-t THREADS] [-p PARALLEL_JOBS]"
  echo "  -i  Input directory with FASTQ files (*.fastq.gz) ending in _R1_001/_R2_001"
  echo "  -o  Output directory"
  echo "  -d  Kraken2 database directory (contains hash.k2d, opts.k2d, taxo.k2d, etc.)"
  echo "  -t  Threads per Kraken2 process (default: 4)"
  echo "  -p  Number of samples to process in parallel (default: 1; requires GNU parallel if >1)"
  exit 1
}

INPUT_DIR=""
OUTPUT_DIR=""
DB_DIR=""
THREADS=4
PARALLEL_JOBS=1

while getopts ":i:o:d:t:p:" opt; do
  case "$opt" in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    d) DB_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    p) PARALLEL_JOBS="$OPTARG" ;;
    *) usage ;;
  esac
done

[ -z "${INPUT_DIR}" ] && usage
[ -z "${OUTPUT_DIR}" ] && usage
[ -z "${DB_DIR}" ] && usage

[ -d "${INPUT_DIR}" ] || { echo "Error: input directory not found: ${INPUT_DIR}" >&2; exit 1; }
mkdir -p "${OUTPUT_DIR}"
LOG_DIR="${OUTPUT_DIR%/}/logs"
REPORT_DIR="${OUTPUT_DIR%/}/reports"
mkdir -p "${LOG_DIR}" "${REPORT_DIR}"

command -v kraken2 >/dev/null 2>&1 || { echo "Error: kraken2 not found in PATH" >&2; exit 1; }
[ -d "${DB_DIR}" ] || { echo "Error: database directory not found: ${DB_DIR}" >&2; exit 1; }
if ! ls "${DB_DIR}"/hash.k2d "${DB_DIR}"/opts.k2d "${DB_DIR}"/taxo.k2d >/dev/null 2>&1; then
  echo "Warning: hash.k2d/opts.k2d/taxo.k2d not found in ${DB_DIR}; verify that the database has been built." >&2
fi

mapfile -d '' R1_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*_R1_001.fastq.gz" -print0)
if [ "${#R1_FILES[@]}" -eq 0 ]; then
  echo "Error: no *_R1_001.fastq.gz files found in ${INPUT_DIR}" >&2
  exit 1
fi
echo "Found ${#R1_FILES[@]} sample(s)."

run_one() {
  local R1="$1"
  local THREADS="$2"
  local DB_DIR="$3"
  local OUTPUT_DIR="$4"
  local LOG_DIR="$5"
  local REPORT_DIR="$6"

  local base="$(basename "$R1" _R1_001.fastq.gz)"
  local R2="${INPUT_DIR%/}/${base}_R2_001.fastq.gz"
  if [ ! -f "$R2" ]; then
    echo "Warning: matching R2 file not found for ${base}; sample skipped." >&2
    return 0
  fi

  local out_txt="${OUTPUT_DIR%/}/${base}_kraken2_output.txt"
  local out_report="${OUTPUT_DIR%/}/${base}_kraken2_report.txt"
  local log_file="${LOG_DIR%/}/${base}.kraken2.log"

  {
    echo "== $(date) =="
    echo "Running Kraken2: ${base}"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "DB: $DB_DIR"
    echo "Threads: $THREADS"
  } >> "$log_file"

  set +e
  kraken2 \
    --db "$DB_DIR" \
    --paired "$R1" "$R2" \
    --threads "$THREADS" \
    --output "$out_txt" \
    --report "$out_report" \
    --use-names \
    >> "$log_file" 2>&1
  status=$?
  set -e

  if [ $status -ne 0 ]; then
    echo "Error: Kraken2 failed for ${base}. See ${log_file}" >&2
    return $status
  fi

  # Consolidate reports into a dedicated folder
  mv -f "$out_report" "${REPORT_DIR%/}/" 2>>"$log_file"

  echo "OK: ${base}"
}

export -f run_one
export INPUT_DIR

if [ "$PARALLEL_JOBS" -gt 1 ]; then
  command -v parallel >/dev/null 2>&1 || { echo "Error: PARALLEL_JOBS=${PARALLEL_JOBS} but GNU parallel is not in PATH" >&2; exit 1; }
  printf '%s\0' "${R1_FILES[@]}" | parallel -0 -j "${PARALLEL_JOBS}" --will-cite \
    run_one {} "${THREADS}" "${DB_DIR}" "${OUTPUT_DIR}" "${LOG_DIR}" "${REPORT_DIR}"
else
  for R1 in "${R1_FILES[@]}"; do
    run_one "$R1" "${THREADS}" "${DB_DIR}" "${OUTPUT_DIR}" "${LOG_DIR}" "${REPORT_DIR}"
  done
fi

echo "Done. Reports in: ${REPORT_DIR}"
echo "Logs in: ${LOG_DIR}"
