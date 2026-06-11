#!/usr/bin/env python3

# LAST UPDATE: 25/05/2026

# AIM: generate a manifest file (manifest.csv) for Snakemake, mapping the sample IDs from the Excel database to their corresponding BAM, VCF, and ichorCNA result files on the HPC server. 
#      The script checks for the existence of these files and creates a CSV file with absolute paths for use in downstream analysis.

# INPUTS: - flow_cytometry_db.xlsx (Excel file with sample information, specifically from the 'db_purity_0.90' sheet)
#         - BAM files located in /mnt/prunity_share/SEQUENCING/AIO_PANEL/BAM/run_{RUN}/remDup_recal_{ID}.bam
#         - VCF files located in /mnt/drive02/UMA/run_{RUN}/Mutect2/mutect2_{ID}.vcf.gz
#         - ichorCNA result files located in /mnt/drive02/ichorCNA_results/run_{RUN}/{ID}_ichor.txt 

# OUTPUT: - manifest.csv (CSV file with columns: sample_id, bam, vcf, ichor_res)


# NOTE: This script is designed to be run in an environment where the specified directories and files exist -- HPC server with the given structure


import pandas as pd
from pathlib import Path

def main():
    print(" Manifest Generator for Snakemake - Starting...")

    # 1. Configuring basic routes on HPC server
    base_dir = Path("/mnt")
    bam_base_folder = base_dir / 'prunity_share' / 'SEQUENCING' / 'AIO_PANEL' / 'BAM'
    vcf_base_folder = base_dir / 'drive02' / 'UMA'
    ichor_base_folder = base_dir / 'drive02' / 'UMA'       # DA MODIFICARE

    # Path to the input Excel file
    excel_file = Path("/mnt/drive02/Tumor_purity_pipeline/flow_cytometry_db.xlsx")
    if not excel_file.exists():
        print(f"ERROR: no excel file found at {excel_file.resolve()}")
        return

    # 2. Reading the Excel database (specific sheet)
    df = pd.read_excel(excel_file, sheet_name='db_purity_no1')
    print(f" Uploaded {len(df)} excel's samples. Processing...")

    # List that will contain the rows of our final manifest
    manifest_rows = []

    # 3. Iteration on samples to map real paths 
    for index, row in df.iterrows():
        sample_id = str(row['ID'])
        run = str(row['RUN'])

        # Building absolute paths using server structure
        bam_path = bam_base_folder / f"run_{run}" / f"remDup_recal_{sample_id}.bam"
        vcf_path = vcf_base_folder / f"run_{run}" / "Mutect2" / f"mutect2_{sample_id}.vcf.gz"
        ichor_path = ichor_base_folder / f"run_{run}" / "ichorCNA" / f"{sample_id}.params.txt"

        # Check the existence of the files 
        if not bam_path.exists():
            print(f"[WARNING] BAM not found for sample {sample_id} in {bam_path}")
            continue
        if not vcf_path.exists():
            print(f"[WARNING] VCF not found for sample {sample_id} in {vcf_path}")
            continue
        if not ichor_path.exists():
            print(f"[WARNING] ichorCNA result not found for sample {sample_id} in {ichor_path}")
            continue

        # If the core files exist, we save the paths in our dictionary
        manifest_rows.append({
            "sample_id": sample_id,
            "bam": str(bam_path.resolve()),
            "vcf": str(vcf_path.resolve()),
            "ichor_res": str(ichor_path.resolve())
        })

    # 4. Creation of the final DataFrame and export to CSV (Comma-separated)
    manifest_df = pd.DataFrame(manifest_rows)
    out_dir = Path("/mnt/drive02/Tumor_purity_pipeline/config/")
    output_csv = "manifest.csv"
    manifest_df.to_csv(out_dir / output_csv, sep=",", index=False)
    
    print(f"Successfully generated manifest! Created file '{output_csv}' with {len(manifest_df)} samples ready.")

if __name__ == "__main__":
    main()