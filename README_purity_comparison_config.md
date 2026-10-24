# Quick guide: updating the manifest and the flow cytometry database

If you want to run the workflow on a different cohort or use a different flow cytometry Excel file, the main changes are usually limited to the manifest path and the input table used for the final purity comparison.

## 1) Using a different manifest

### Files to update
- `Snakefile`
- `workflow/scripts/purity_comparisonOOP.py`

### Variables to change

In `Snakefile`:
- update the initial manifest loaded at the beginning of the workflow
  - example: `pd.read_csv("config/manifest_normali.csv")`
- update the manifest used in the `purity_comparison` rule
  - example: `manifest = "config/manifest_normali.csv"`

In `workflow/scripts/purity_comparisonOOP.py`:
- the script reads the manifest using `pd.read_csv(self.manifest_path)`
- it verifies that the `sample_id` column is present
- it ensures that `quality_flag` exists, or creates it automatically if missing

### Required columns in the manifest
The manifest should contain at least these columns:

- `sample_id`
- `bam`
- `vcf`
- `ichor_res`

If your file uses different column names, update the parsing logic in the Python script accordingly.

---

## 2) Using a different flow cytometry Excel file

### Files to update
- `Snakefile`
- `workflow/scripts/purity_comparisonOOP.py`

### Variables to change

In `Snakefile`:
- update the Excel path used in the `purity_comparison` rule
  - example: `flow_db = "flowdb_normali.xlsx"`

In `workflow/scripts/purity_comparisonOOP.py`:
- if the correct sheet is not the first one, specify it explicitly with `pd.read_excel(..., sheet_name=...)`
- verify that sample IDs in the Excel file match the manifest sample IDs
- ensure the purity column is selected and renamed to `purity_flow`

### Required columns in the flow cytometry file
The Excel file should contain at least:

- `ID`
- `purity`

If the column names differ, update the script to match the new names.

---

## 3) Minimal example structure

### Manifest CSV
```csv
sample_id,bam,vcf,ichor_res
sample001,/path/to/sample001.bam,/path/to/sample001.vcf.gz,/path/to/sample001.params.txt
sample002,/path/to/sample002.bam,/path/to/sample002.vcf.gz,/path/to/sample002.params.txt
