# !/usr/bin/r

# ASCAT.sc v.0.1 - compatible with xgboost v.1.7.3.1: install.packages("xgboost", repos = "https://packagemanager.posit.co/cran/2023-03-01/")

# ASCAT.sc function
#   run_sc_sequencing: for ULP (from <0.1X to >10X); count reads are read at SNP positions, not bins for the logR
#   run_targeted_sequencing: leverages off-target reads from exome/targeted panel to derive genome-wide CN profiles, should complement focal CN profiles. need to
#   derive genomic bins and GC-content of each bin before running the algorithm

# INPUT
#   tumor BAM files: 
#   normal BAM files: not required to be coupled, but females or normals without sexual chrs are preferred (the larger the better)
#   refGen FASTA format
#   BED file with regions to be excluded from the analyses - derived from normal by excluding targeted regions


# data format
# default purity [vector, list]: lower bound of 10% - can be changed if required
# ploidy [vector, list]

# IMPORTANTE --> run_targeted_sequencing.R > printResults_all.R > plot_Results.R

library(ASCAT.sc)
library(data.table)
library(tidyverse)

# ---- TAKE THE ARGUMENTS FROM THE SNAKEFILE ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Error: 4 arguments needed: Patient_ID, Path_BAM, Path_Normal, Out_Dir")
}

# ---- SETUP VARIABLES ----
sample_id <- args[1]    # Patient ID, e.g., "6922_S6"
bam_path  <- args[2]    # Path to the tumor BAM file, e.g., "C:/Users/rossi/ASCAT_UMA/bams_purity0.9/remDup_recal_6922_S6.bam"
normal_bam <- args[3]   # Path to the normal BAM file, e.g., "C:/Users/rossi/ASCAT_UMA/normal_bam/remDup_recal_DONOR.bam"
outDir    <- args[4]    # Output directory, e.g., "C:/Users/rossi/ASCAT_UMA/output_6922_S6"
bed       <- "./targets_hg19.bed"

message("--- Execution for Sample: ", sample_id, " ---")
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

projectName <- paste0('project_', sample_id)

# Running ASCAT.sc (ON ONE BAM ONLY)
res <- run_targeted_sequencing(tumour_bams = bam_path, bed_file = bed, build = "hg19",
                               normal_bams = normal_bam, projectname = projectName, 
                               outdir = outDir, predict_refit=FALSE, print_results = FALSE)

printResults_all(res, outdir = outDir, projectname = projectName)

# ---- PLOT THE RESULTS ----
profile_file <- list.files(
  path = file.path(outDir, basename(bam_path)),
  pattern = paste0(".*", sample_id, ".*_project.*\\.ASCAT\\.scprofile\\.txt$"), 
  full.names = TRUE
)[1]

if (!is.na(profile_file) && file.exists(profile_file)) {
  plot_ascat_sc_segments(profile_file, sample_id, outDir)
} else {
  warning("WARNING: .scprofile file not found")
}