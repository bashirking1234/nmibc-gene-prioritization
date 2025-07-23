# ------------------------------------------------------------------------------
# Script Name: eqtl_lookup.R
#
# Description:
#   This script performs a lookup to identify overlapping genetic variants between
#   a TensorQTL eQTL file and a metaGWAS file. It harmonizes variant identifiers,
#   merges the datasets based on chromosome, position, and alleles, and outputs
#   the matched variants with relevant annotation.
#
# Usage:
#   Rscript eqtl_lookup.R <eqtl_input> <gwas_input> <output_file>
#     - <eqtl_input>:    Path to the TensorQTL eQTL input file (tab-delimited).
#     - <gwas_input>:    Path to the metaGWAS input file (tab-delimited).
#     - <output_file>:   Path to the output file for matched variants (tab-delimited).
#
# Arguments:
#   eqtl_input   : Input file containing eQTL summary statistics.
#   gwas_input   : Input file containing metaGWAS summary statistics.
#   output_file  : Output file to write the filtered, matched variants.
#
# Details:
#   - The script expects both input files to contain variant information in the
#     format: Chr:Position:Allele1:Allele2 (with or without 'chr' prefix).
#   - It harmonizes the variant identifiers, splits them into separate columns,
#     and merges the datasets on these columns.
#   - Only variants present in both datasets are retained in the output.
#   - The output file contains selected columns, reordered if present.
#
# Dependencies:
#   - data.table
#   - R.utils
#
# Example:
#   Rscript eqtl_lookup.R eqtl_data.txt gwas_data.txt matched_variants.txt
#
#' @author Bashir Hussein 
#' @date 2025-01-07
#'
# ------------------------------------------------------------------------------
#!/usr/bin/env Rscript

library(data.table)
library(R.utils)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript eqtl_lookup.R <eqtl_input> <gwas_input> <output_file>")
}

eqtl_input <- args[1]
gwas_input <- args[2]
output_file <- args[3]

# Load TensorQTL eQTL file
eqtl <- fread(eqtl_input, sep = "\t", header = TRUE)
eqtl[, c("Chr", "Position", "Allele1", "Allele2") := tstrsplit(variant_id, ":", keep = c(1, 2, 3, 4))]
eqtl[, Chr := sub("chr", "", Chr)]
eqtl[, Position := as.integer(Position)]

# Load metaGWAS file
gwas <- fread(gwas_input, sep = "\t", header = TRUE)
gwas[, c("Chr", "Position", "Allele1", "Allele2") := tstrsplit(MarkerName, ":", keep = c(1, 2, 3, 4))]
gwas[, Chr := sub("chr", "", Chr)]
gwas[, Position := as.integer(Position)]

# Harmonize types
gwas[, Chr := as.character(Chr)]
eqtl[, Chr := as.character(Chr)]

# Merge to find matching variants
common <- merge(gwas[, .(Chr, Position, Allele1, Allele2)],
                eqtl[, .(Chr, Position, Allele1, Allele2)],
                by = c("Chr", "Position", "Allele1", "Allele2"))
cat("Number of overlapping variants found:", nrow(common), "\n")

# Merge full data
merged <- merge(eqtl,
                gwas[, .(MarkerName, Chr, Position, Allele1, Allele2)],
                by = c("Chr", "Position", "Allele1", "Allele2"),
                all.x = TRUE)

# Rename to avoid confusion
setnames(merged, "MarkerName", "MarkerName_data")

# Filter out NAs
filtered <- merged[!is.na(MarkerName_data)]

# Write output
fwrite(filtered, output_file, sep = "\t")
