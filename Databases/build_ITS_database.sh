#!/bin/bash
#
# build_ITS_database.sh
#
# Builds a Kraken2 database for ITS (fungal) taxonomic classification,
# combining the UNITE ITS reference database with the NCBI fungi genome
# library.
#
# Output: kraken2_db_ITS/ — pass this path to the -d flag of
# 04_taxonomic_classification.sh.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Requirements: kraken2-build in PATH, tar
#
# IMPORTANT: this script does NOT download the UNITE database
# automatically. Download it manually beforehand (UNITE general FASTA
# release, "dynamic" files) from:
#   https://doi.plutof.ut.ee/doi/10.15156/BIO/2959332
# and place the archive (e.g. sh_general_release_04.04.2024.tgz) in this
# directory before running this script.
set -euo pipefail

DB_NAME="kraken2_db_ITS"
UNITE_ARCHIVE="sh_general_release_04.04.2024.tgz"
UNITE_FASTA="sh_general_release_dynamic_04.04.2024.fasta"

command -v kraken2-build >/dev/null 2>&1 || { echo "Error: kraken2-build not found in PATH" >&2; exit 1; }

if [ ! -f "$UNITE_ARCHIVE" ]; then
  echo "Error: UNITE archive not found: $UNITE_ARCHIVE" >&2
  echo "Download it manually from https://doi.plutof.ut.ee/doi/10.15156/BIO/2959332" >&2
  echo "and place it in this directory before running this script." >&2
  exit 1
fi

mkdir -p "$DB_NAME"

echo "Extracting UNITE database..."
tar -xf "$UNITE_ARCHIVE" -C ./

if [ ! -f "$UNITE_FASTA" ]; then
  echo "Error: expected FASTA file not found after extraction: $UNITE_FASTA" >&2
  echo "Check the contents of $UNITE_ARCHIVE and adjust UNITE_FASTA accordingly." >&2
  exit 1
fi

echo "Downloading NCBI taxonomy (this takes a few minutes)..."
kraken2-build --download-taxonomy --db "$DB_NAME"

echo "Adding UNITE sequences to the library..."
kraken2-build --add-to-library "$UNITE_FASTA" --db "$DB_NAME"

echo "Downloading fungi genome library..."
kraken2-build --download-library fungi --db "$DB_NAME"

echo "Building the database (this takes a few minutes)..."
kraken2-build --build --db "$DB_NAME"

echo "Removing intermediate library files..."
kraken2-build --clean --db "$DB_NAME"

echo "Verifying database..."
kraken2-inspect --db "$DB_NAME"

echo "Done. Database ready at: $DB_NAME"
