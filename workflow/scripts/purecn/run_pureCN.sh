#!/bin/bash

# --- PERCORSI DI SISTEMA E RIFERIMENTI FISSI ---
export PURECN="/usr/local/lib/R/site-library/PureCN/extdata"
export INTERVALS="/mnt/drive02/PureCN/output/intervals_ref/baits_hg19_intervals.txt"
export NORMALDB="/mnt/drive02/PureCN/output/normalDB/normalDB_hg19.rds"

# --- LISTA INPUT E ROOT OUTPUT ---
export TUMORS_LIST="/mnt/drive02/PureCN/tumors_bam.list"
export COV_BASE="/mnt/drive02/PureCN/output/coverage_output/samples"
export RESULTS_BASE="/mnt/drive02/PureCN/output/purecn_results"

echo "Starting PureCN batch processing..."

# Reads the previously created BAM list and processes each sample iteratively
while read -r BAM_PATH; do
    # 1. Extracts the full name from the path (e.g. remDup_recal_6959_S2)
    FULL_NAME=$(basename "$BAM_PATH" .bam)
    
    # 2. Optimizes SAMPLEID: removes "remDup_recal_" to have only the number (e.g. 6959_S2) using sed
    SAMPLEID=$(echo "$FULL_NAME" | sed 's/remDup_recal_//')
    
    echo "------------------------------------------------------------"
    echo "Sample Analysis: $SAMPLEID (File: $FULL_NAME)"
    echo "------------------------------------------------------------"
    
    # 3. Defines the necessary file paths based on the full name
    # Coverage is located in the folder created in Step 2
    TUMOR_COVERAGE="$COV_BASE/$FULL_NAME/${FULL_NAME}_coverage.txt.gz"
    
    # Assumendo che i VCF siano in questa cartella e abbiano il nome completo
    VCF_FILE="/mnt/drive02/PureCN/input_file/vcf/mutect2_${SAMPLEID}.vcf.gz"
    
    # 4. Crea la cartella di output specifica per il SAMPLEID pulito
    OUT_DIR="$RESULTS_BASE/$SAMPLEID"
    mkdir -p "$OUT_DIR"
    
    # 5. Esecuzione PureCN
    # Nota: Usiamo $SAMPLEID pulito per il flag --sampleid e per il nome del file --out
    Rscript "$PURECN/PureCN.R" \
        --out "$OUT_DIR/${SAMPLEID}" \
        --tumor "$TUMOR_COVERAGE" \
        --sampleid "$SAMPLEID" \
        --vcf "$VCF_FILE" \
        --normaldb "$NORMALDB" \
        --intervals "$INTERVALS" \
        --genome hg19 \
        --model betabin \
        --fun-segmentation PSCBS \
        --post-optimize \
        --force \
        --seed 123 \
      

done < "$TUMORS_LIST"

echo "Batch processing completed successfully."