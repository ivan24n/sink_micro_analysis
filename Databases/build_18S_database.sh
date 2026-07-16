#!/bin/bash
#
# build_18S_database.sh
#
# Builds a custom Kraken2 database for 18S rRNA (eukaryotic) taxonomic
# classification, combining three reference sequence sets (EukRibo,
# SILVA, PR2) plus the NCBI protozoa genome library.
#
# Output: kraken2_db_18S/ — pass this path to the -d flag of
# 04_taxonomic_classification.sh.
#
# Part of the analysis pipeline for:
#   "Spatiotemporal dynamics of multi-kingdom microbial communities in
#    hospital sinks"
#
# Requirements: kraken2-build in PATH, wget, gunzip
set -euo pipefail

DB_NAME="kraken2_db_18S"

command -v kraken2-build >/dev/null 2>&1 || { echo "Error: kraken2-build not found in PATH" >&2; exit 1; }
command -v wget >/dev/null 2>&1 || { echo "Error: wget not found in PATH" >&2; exit 1; }

mkdir -p "$DB_NAME"

# ==========================================
# 1. Download reference sequences
# ==========================================

echo "Downloading EukRibo..."
wget -O 46346_EukRibo-02_full_seqs_2022-07-22.fas.gz \
  "https://zenodo.org/records/6896896/files/46346_EukRibo-02_full_seqs_2022-07-22.fas.gz?download=1"
gunzip -f 46346_EukRibo-02_full_seqs_2022-07-22.fas.gz

echo "Downloading SILVA..."
wget -O silva_database.fasta.gz \
  "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
gunzip -f silva_database.fasta.gz

echo "Downloading PR2..."
wget -O pr2_version_5.0.0_SSU_taxo_long.fasta.gz \
  "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz"
gunzip -f pr2_version_5.0.0_SSU_taxo_long.fasta.gz

# ==========================================
# 2. NCBI taxonomy and Kraken2 library
# ==========================================

# --use-ftp avoids rsync connections being blocked on some networks.
echo "Downloading NCBI taxonomy (via FTP/HTTP)..."
kraken2-build --download-taxonomy --use-ftp --db "$DB_NAME"

echo "Adding downloaded sequences to the library..."
kraken2-build --add-to-library pr2_version_5.0.0_SSU_taxo_long.fasta --db "$DB_NAME"
kraken2-build --add-to-library silva_database.fasta --db "$DB_NAME"
kraken2-build --add-to-library 46346_EukRibo-02_full_seqs_2022-07-22.fas --db "$DB_NAME"

echo "Downloading protozoa library..."
kraken2-build --download-library protozoa --use-ftp --db "$DB_NAME"

# ==========================================
# 3. Build the database
# ==========================================
echo "Building the database (this will take a while)..."
kraken2-build --build --db "$DB_NAME"

# Uncomment to remove intermediate library files after building (keeps
# only the files needed for classification):
# kraken2-build --clean --db "$DB_NAME"

echo "Verifying database..."
kraken2-inspect --db "$DB_NAME"

echo "Done. Database ready at: $DB_NAME"
