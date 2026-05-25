# !/bin/bash

# --- PURECN PACKAGE PATH ---
export PURECN="/usr/local/lib/R/site-library/PureCN/extdata"

# --- CONFIG PATHS ---
export GENOME="/mnt/drive02/PureCN/input_file/intervals/hs37d5.fa"
export BAITS="/mnt/drive02/PureCN/input_file/intervals/all_probes_targets_hg19_clean.bed"
export MAP="/mnt/drive02/PureCN/input_file/intervals/mappabilityAlign100mer.bigWig"
export OUT_REF="/mnt/drive02/PureCN/output/intervals_ref"


# --- CREATE OUTPUT DIR ---
mkdir -p $OUT_REF


# --- PURECN INTERVALS.R EXECUTION ---
Rscript $PURECN/IntervalFile.R \
	--in-file $BAITS \
	--fasta $GENOME \
	--out-file $OUT_REF/baits_hg19_intervals.txt \
	--off-target \
	--genome hg19 \
	--mappability $MAP \
	--export $OUT_REF/baits_optimized_h19.bed

