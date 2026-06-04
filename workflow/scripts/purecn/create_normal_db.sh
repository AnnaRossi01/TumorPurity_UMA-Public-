# 1. Create a list of the fully processed GC-normalized coverage files
ls /mnt/drive02/PureCN/output/coverage_output/normals_output/*_coverage.txt.gz > /mnt/drive02/PureCN/normal_coverages.list

# 2. Build the Database
export PURECN="/usr/local/lib/R/site-library/PureCN/extdata"
export NORMAL_LIST="/mnt/drive02/PureCN/normal_coverages.list"
export OUT_DIR="/mnt/drive02/PureCN/output/normalDB"

mkdir -p "$OUT_DIR"

Rscript "$PURECN/NormalDB.R" \
  --out-dir "$OUT_DIR" \
  --coverage-files "$NORMAL_LIST" \
  --genome hg19