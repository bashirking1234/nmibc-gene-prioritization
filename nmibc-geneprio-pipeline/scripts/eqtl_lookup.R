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
eqtl[, variant_id := gsub("_", ":", variant_id)]
eqtl[, variant_id := sub(":b38", "", variant_id)]
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
cat("Aantal overlappende varianten gevonden:", nrow(common), "\n")

# Merge full data
merged <- merge(eqtl,
                gwas[, .(MarkerName, Chr, Position, Allele1, Allele2)],
                by = c("Chr", "Position", "Allele1", "Allele2"),
                all.x = TRUE)

# Rename to avoid confusion
setnames(merged, "MarkerName", "MarkerName_data")

# Reorder columns if ze bestaan
desired_cols <- c("Cancer_type", "Gene_name", "RS_id", "Variant", "Chr",
                  "MarkerName_data", "Position", "Allele1", "Allele2",
                  "Credible_Set", "Credible_Set_Size", "PIP")
existing_cols <- desired_cols[desired_cols %in% colnames(merged)]
setcolorder(merged, existing_cols)

# Filter out NAs
filtered <- merged[!is.na(MarkerName_data)]

# Write output
fwrite(filtered, output_file, sep = "\t")
