#!/bin/bash

# --- CONFIG PATHS ---
export PURECN="/usr/local/lib/R/site-library/PureCN/extdata"
export INTERVALS="/mnt/drive02/PureCN/output/intervals_ref/baits_hg19_intervals.txt"

# --- TUMORS PATHS ---
# CORRETTO: Ora punta direttamente alla tua nuova cartella isolata
export TUMORS_DIR="/mnt/drive02/PureCN/input_file/bam/tumors"
export BASE_OUT="/mnt/drive02/PureCN/output/coverage_output/samples"
export TUMORS_LIST="/mnt/drive02/PureCN/tumors_bam.list"

mkdir -p "$BASE_OUT"

echo "Looking for the BAM files of the tumors..."
# Cerca tutti i .bam SOLO al primo livello della cartella tumors/
find "$TUMORS_DIR" -maxdepth 1 -name "*.bam" > "$TUMORS_LIST"
NUM_TUMORS=$(wc -l < "$TUMORS_LIST")
echo "Found $NUM_TUMORS tumor samples to process."

# --- PROCESSING LOOP ---
while read -r BAM; do
    # Estrae il nome pulito del campione (es. rimuove path e .bam)
    SAMPLE=$(basename "$BAM" .bam)
    
    echo "---------------------------------------------------"
    echo " Start of sample processing: $SAMPLE "
    echo "---------------------------------------------------"
    
    # Create the specific folder for this patient
    SAMPLE_OUT="$BASE_OUT/$SAMPLE"
    mkdir -p "$SAMPLE_OUT"
    
    # Lancia lo script salvando nella cartuccia dedicata
    Rscript "$PURECN/Coverage.R" \
      --out-dir "$SAMPLE_OUT" \
      --bam "$BAM" \
      --intervals "$INTERVALS"
      
done < "$TUMORS_LIST"

echo "Processing of all tumor samples completed."