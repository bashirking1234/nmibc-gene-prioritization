.libPaths("Q:\\Stagiairs\\Bashir Hussein.W\\R")


# The GRCh37 genome assembly is used for the coordinates of genetic variants and genes in this database.
# https://hanlaboratory.com/PancanQTLv2//Documents.html

# The TCGA data were filtered for PIP & FRD
# Posterior inclusion probabilities (PIPs) >= 90% for all SNP's. it means that the sum of the PIPs for all SNPs in the set is at least 0.90
# This means that the researchers selected SNPs based on their PIP scores and included enough of them so that their combined PIPs add up to 0.90 or more.
# PIPs help researchers identify and prioritize genetic variants that are most likely to have a significant impact on traits or diseases
# false discovery rate (FDR) < 0.05
# FDR is used to control for false positives, This means they selected genes where the likelihood of false positives was controlled to be less than 5%
#
#

# ------------------------------------------------
# Load necessary library
# ------------------------------------------------
library(data.table)

# ------------------------------------------------
# Load TCGA data
# ------------------------------------------------
tcga_file_path <- "Q:/Stagiairs/Bashir Hussein.W/Data/TCGA/BLCA.cis.susie_with_converted_variant.txt"
tcga_eqtl <- fread(tcga_file_path, sep = "\t", header = TRUE)

# ------------------------------------------------
# Load metaGWAS data (reused from previous context)
# ------------------------------------------------
metaGWAS_file_path <- "Z:/genepi/Bashir_Hussein/metaGWAS_outputbestand/MA_STDERR_gwas_prognose_update_maf005_r208_recur_1.tbl"
metaGWAS_data <- fread(metaGWAS_file_path, sep = "\t", header = TRUE)

# ------------------------------------------------
# Preprocessing TCGA eQTL data
# ------------------------------------------------

# Standardize variant format and extract components
tcga_eqtl[, converted_variant_b38 := gsub("_", ":", converted_variant_b38)]
tcga_eqtl[, converted_variant_b38 := paste0("chr", converted_variant_b38)]
tcga_eqtl[, c("Chr", "Position", "Allele1", "Allele2") := tstrsplit(converted_variant_b38, ":", keep = c(1, 2, 3, 4))]
tcga_eqtl[, Chr := sub("chr", "", Chr)]
tcga_eqtl[, Position := as.integer(Position)]

# ------------------------------------------------
# Preprocessing metaGWAS data
# ------------------------------------------------

# Split MarkerName into components
metaGWAS_data[, c("Chr", "Position", "Allele1", "Allele2") := tstrsplit(MarkerName, ":", keep = c(1, 2, 3, 4))]
metaGWAS_data[, Chr := sub("chr", "", Chr)]
metaGWAS_data[, Position := as.integer(Position)]


# ------------------------------------------------
# Ensure matching data types
# ------------------------------------------------
metaGWAS_data[, Chr := as.character(Chr)]
tcga_eqtl[, Chr := as.character(Chr)]
metaGWAS_data[, Position := as.integer(Position)]
tcga_eqtl[, Position := as.integer(Position)]

# ------------------------------------------------
# Merge Step 1: Merge on Chr and Position only
# ------------------------------------------------
merged_data_chr_pos <- merge(
  tcga_eqtl, 
  metaGWAS_data[, .(Chr, Position, MarkerName, Allele1, Allele2)], 
  by = c("Chr", "Position"), 
  all = FALSE
)

# ------------------------------------------------
# Step 2: Filter rows where alleles match between both datasets
# ------------------------------------------------
merged_data_filtered_alleles <- merged_data_chr_pos[
  (Allele1.x == Allele1.y & Allele2.x == Allele2.y) |
    (Allele1.x == Allele2.y & Allele2.x == Allele1.y)
]

# ------------------------------------------------
# Rename MarkerName to variant_id_metaGWAS
# ------------------------------------------------
setnames(merged_data_filtered_alleles, "MarkerName", "variant_id_metaGWAS")

# ------------------------------------------------
# Optional: Reorder columns for clarity
# ------------------------------------------------
setcolorder(
  merged_data_filtered_alleles, 
  c("Cancer_type", "Gene_name", "RS_id", "Variant", "variant_id_metaGWAS", "Position", "Allele1.x", "Allele2.x", "Credible_Set", "Credible_Set_Size", "PIP")
)


# ------------------------------------------------
# Summary
# ------------------------------------------------
cat("Total shared SNPs on Chr:Pos:", nrow(merged_data_chr_pos), "\n")
cat("After allele matching:", nrow(merged_data_filtered_alleles), "\n")

# ------------------------------------------------
# Optional: Filter out rows with NA in variant_id_metaGWAS
# ------------------------------------------------
filtered_final_tcga <- merged_data_filtered_alleles[!is.na(variant_id_metaGWAS)]
