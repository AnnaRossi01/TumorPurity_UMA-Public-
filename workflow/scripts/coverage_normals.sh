#!/bin/bash

# --- CONFIG PATHS ---
export PURECN="/usr/local/lib/R/site-library/PureCN/extdata"
export INTERVALS="/mnt/drive02/PureCN/output/intervals_ref/baits_hg19_intervals.txt"

# --- NORMALS PATHS ---
export NORMALS_DIR="/mnt/drive02/PureCN/input_file/bam/normals"
export OUT_NORMALS="/mnt/drive02/PureCN/output/coverage_output/normals_output"
export NORMALS_LIST="/mnt/drive02/PureCN/normals_bam.list"

mkdir -p "$OUT_NORMALS"

echo "Looking for the BAM files of the normal ones..."
find "$NORMALS_DIR" -name "*.bam" > "$NORMALS_LIST"
NUM_NORMALS=$(wc -l < "$NORMALS_LIST")
echo "Found $NUM_NORMALS normals. Starting coverage extraction in batch..."

# --- PURECN COVERAGE.R EXECUTION ---
# Il flag --bam accetta il file .list e --cores velocizza il processo
Rscript "$PURECN/Coverage.R" \
  --out-dir "$OUT_NORMALS" \
  --bam "$NORMALS_LIST" \
  --intervals "$INTERVALS" \
  --cores 8

echo "Coverage extraction for normals completed."