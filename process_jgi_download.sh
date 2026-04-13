#!/bin/bash
# Process JGI downloaded zip file

BASE="/data/habib/EcoSurfBGC"
ZIPFILE="${BASE}/genomes/download.20260315.045748.zip"
JGI_DIR="${BASE}/genomes/jgi_download"
FASTA_DIR="${BASE}/genomes/fasta"
LOGFILE="${BASE}/logs/jgi_processing.log"

echo "=== Starting JGI download processing at $(date) ===" | tee "${LOGFILE}"

# Step 1: Check zip size and unzip
echo "Zip file size: $(ls -lh ${ZIPFILE} | awk '{print $5}')" | tee -a "${LOGFILE}"
echo "Unzipping..." | tee -a "${LOGFILE}"

unzip -q "${ZIPFILE}" -d "${JGI_DIR}"
echo "Unzip complete. Contents:" | tee -a "${LOGFILE}"
ls "${JGI_DIR}" | head -20 | tee -a "${LOGFILE}"
echo "Total items: $(ls ${JGI_DIR} | wc -l)" | tee -a "${LOGFILE}"

# Step 2: Find all FASTA files recursively
echo "" | tee -a "${LOGFILE}"
echo "=== Finding FASTA files ===" | tee -a "${LOGFILE}"
find "${JGI_DIR}" -name "*.fna" -o -name "*.fasta" -o -name "*.fa" | \
  head -20 | tee -a "${LOGFILE}"
echo "Total FASTA-like files: $(find ${JGI_DIR} \
  -name '*.fna' -o -name '*.fasta' -o -name '*.fa' | wc -l)" | \
  tee -a "${LOGFILE}"

# Step 3: Check what file types are present
echo "" | tee -a "${LOGFILE}"
echo "=== File types in download ===" | tee -a "${LOGFILE}"
find "${JGI_DIR}" -type f | sed 's/.*\.//' | sort | uniq -c | \
  sort -rn | head -20 | tee -a "${LOGFILE}"

echo "" | tee -a "${LOGFILE}"
echo "=== Directory structure (2 levels) ===" | tee -a "${LOGFILE}"
find "${JGI_DIR}" -maxdepth 2 -type d | head -30 | tee -a "${LOGFILE}"
