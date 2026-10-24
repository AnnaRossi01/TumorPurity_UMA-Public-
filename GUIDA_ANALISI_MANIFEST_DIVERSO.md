# Quick guide: run analysis with a new manifest

This is a short English guide to run the tumor purity workflow with a new cohort or a different sample manifest.

## 1. Prepare the inputs

Make sure the following exist:
- a manifest CSV in `config/`
- a flow cytometry Excel file
- BAM and VCF files for each sample
- ichorCNA result files
- the reference configuration in `config/config.yaml`

The manifest should contain at least:

```csv
sample_id,bam,vcf,ichor_res
sample001,/path/to/sample001.bam,/path/to/sample001.vcf.gz,/path/to/sample001.params.txt
```

## 2. Generate the manifest

If your cohort metadata is stored in an Excel file, generate the manifest with:

```bash
python workflow/scripts/generate_manifest.py
```

Check the output file in `config/` and verify that:
- `sample_id` values are unique,
- BAM, VCF and ichorCNA paths exist,
- IDs match the flow-cytometry table,
- the number of samples matches your expectation.

## 3. Update the workflow

Edit `Snakefile` and replace the default manifest path:

```python
manifest_df = pd.read_csv("config/manifest_95.csv")
```

with your new file, for example:

```python
manifest_df = pd.read_csv("config/manifest_25.csv")
```

Also update the manifest used in the `purity_comparison` rule and set a new output directory if needed.

## 4. Update the flow cytometry reference

In `Snakefile`, make sure the flow cytometry file path points to the correct Excel file:

```python
params:
    flow_db = "db_last.xlsx"
```

The Excel file should contain at least:

```text
ID,purity
sample001,0.65
```

## 5. Dry run and execution

Before running the full workflow, check the DAG:

```bash
snakemake --dry-run --cores 4
```

Then run the analysis:

```bash
snakemake --cores 4
```

The final summary files are typically written under:

```text
results/purity_comparison_<name>/
```

## 6. Main outputs

Look for:

```text
results/purity_comparison_<name>/purity_comparison_<name>.csv
results/purity_comparison_<name>/plots/
results/purity_comparison_<name>/quality_filter_report.txt
```

These files contain the merged tumor purity comparison and the associated QC information.
