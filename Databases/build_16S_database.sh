#!/bin/bash
#
# build_16S_database.sh
#
# Downloads a pre-built Kraken2 database for 16S rRNA taxonomic
# classification, based on the RDP (Ribosomal Database Project) 16S
# reference collection distributed via the genome-idx project. This
# database is already built (no --build step is needed), so this script
# only downloads it, extracts it, and (optionally) removes intermediate
# library files.
#
# Output: kraken2_db_16S/ — pass this path to the -d flag of
# 04_taxonomic_classification.sh.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Requirements: kraken2-build in PATH, wget, tar
set -euo pipefail

DB_NAME="kraken2_db_16S"
RDP_URL="https://genome-idx.s3.amazonaws.com/kraken/16S_RDP11.5_20200326.tgz"
RDP_ARCHIVE="16S_RDP11.5_20200326.tgz"
# Name of the folder contained inside the downloaded archive.
RDP_EXTRACTED_DIR="16S_RDP_k2db"

command -v kraken2-build >/dev/null 2>&1 || { echo "Error: kraken2-build not found in PATH" >&2; exit 1; }
command -v wget >/dev/null 2>&1 || { echo "Error: wget not found in PATH" >&2; exit 1; }

echo "Downloading pre-built RDP 16S Kraken2 database..."
wget -O "$RDP_ARCHIVE" "$RDP_URL"

echo "Extracting database..."
tar -xzf "$RDP_ARCHIVE"

if [ ! -d "$RDP_EXTRACTED_DIR" ]; then
  echo "Error: expected extracted folder '$RDP_EXTRACTED_DIR' not found." >&2
  echo "Check the contents of $RDP_ARCHIVE and adjust RDP_EXTRACTED_DIR accordingly." >&2
  exit 1
fi

mv "$RDP_EXTRACTED_DIR" "$DB_NAME"

echo "Removing intermediate library files (database is already built)..."
kraken2-build --clean --db "$DB_NAME"

echo "Verifying database..."
kraken2-inspect --db "$DB_NAME"

echo "Done. Database ready at: $DB_NAME"
