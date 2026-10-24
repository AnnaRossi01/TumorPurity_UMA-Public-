# Quick guide: run PureCN on a single sample

This is a short English guide for running PureCN on a single sample without executing the full tumor-purity comparison workflow.

## 1. Create a minimal manifest

Create a manifest with one row only, for example:

```csv
sample_id,bam,vcf,ichor_res
sample001,/path/to/sample001.bam,/path/to/sample001.vcf.gz,/path/to/sample001.params.txt
```

Save it in `config/`, for example:

```text
config/manifest_sample001.csv
```

## 2. Check the configuration

Open `config/config.yaml` and verify that the following paths are valid:

```yaml
refGen: ./resources/reference/hg19/hs37d5.fa
bed_file: ./resources/reference/hg19/targets/targets_hg19.bed
normal: ./resources/reference/hg19/donor/remDup_recal_DONOR.bam
baits_clean: ./resources/purecn/all_probes_targets_hg19_clean.bed
mappability: ./resources/purecn/mappabilityAlign100mer.bigWig
```

These values define the reference genome, panel targets, mappability track, and the normal reference used by the workflow.

## 3. Point the workflow to the single-sample manifest

In `Snakefile`, replace the default manifest with your single-sample manifest:

```python
manifest_df = pd.read_csv("config/manifest_sample001.csv")
```

If needed, also update the `manifest` field used in the `purity_comparison` rule.

## 4. Run a dry run

From the project root, run:

```bash
snakemake --dry-run --cores 4 results/purecn/sample001/sample001.csv
```

If the dry run looks correct, launch the job:

```bash
snakemake --cores 4 results/purecn/sample001/sample001.csv
```

## 5. Main output

The main output for the sample is:

```text
results/purecn/sample001/sample001.csv
```

You should also check:

```text
results/purecn/coverage/sample001/sample001_coverage.txt.gz
logs/purecn/sample001_run.log
```

## 6. Interpret the result

The key column is `purity`.

If the file contains `NA`, check the log file because it usually indicates that PureCN did not converge or did not find sufficient reliable variant signal.

## 7. Minimal command summary

```bash
snakemake --dry-run --cores 4 results/purecn/sample001/sample001.csv
snakemake --cores 4 results/purecn/sample001/sample001.csv
cat results/purecn/sample001/sample001.csv
less logs/purecn/sample001_run.log
```

This gives a minimal and controlled way to run PureCN on a single sample while keeping the rest of the workflow unchanged.
