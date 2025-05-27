.libPaths("Q:\\Stagiairs\\Bashir Hussein.W\\R")

library(data.table)
library(R.utils)
library(dplyr)



#-------------------------------------------------
# Load the data for GTEx & metaGWAS
# https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/README_eQTL_v10.txt
#-------------------------------------------------

bladder_eqtl_file <- "Q:/Stagiairs/Bashir Hussein.W/Data/GTEx/GTEx_Analysis_v10_eQTL_updated/GTEx_Analysis_v10_eQTL_updated/Bladder.v10.eGenes.txt.gz"
gtex_eqtl <- fread(file = bladder_eqtl_file)

metaGWAS_file_path <- "Z:/genepi/Bashir_Hussein/metaGWAS_outputbestand/MA_STDERR_gwas_prognose_update_maf005_r208_recur_1.tbl"
metaGWAS_data <- fread(metaGWAS_file_path, sep = "\t", header = TRUE)

# -----------------------------------------------
# Preprocessing GTEx eQTL data
# -----------------------------------------------

# Format variant_id: change underscores to colons and remove ":b38"
gtex_eqtl[, variant_id := gsub("_", ":", variant_id)]
gtex_eqtl[, variant_id := sub(":b38$", "", variant_id)]

# Convert 'af' column to numeric and filter for MAF ≥ 5%
gtex_eqtl$af <- as.numeric(gtex_eqtl$af)
gtex_eqtl <- gtex_eqtl[af >= 0.05]

# Split variant_id into Chr, Position, Allele1, Allele2
gtex_eqtl[, c("Chr", "Position", "Allele1", "Allele2") := tstrsplit(variant_id, ":", keep = c(1, 2, 3, 4))]
gtex_eqtl[, Chr := sub("chr", "", Chr)]
gtex_eqtl[, Position := as.integer(Position)]

# -----------------------------------------------
# Preprocessing metaGWAS data
# -----------------------------------------------

# Split MarkerName into Chr, Position, Allele1, Allele2
metaGWAS_data[, c("Chr", "Position", "Allele1", "Allele2") := tstrsplit(MarkerName, ":", keep = c(1, 2, 3, 4))]
metaGWAS_data[, Chr := sub("chr", "", Chr)]
metaGWAS_data[, Position := as.integer(Position)]

# Ensure matching data types
metaGWAS_data[, Chr := as.character(Chr)]
gtex_eqtl[, Chr := as.character(Chr)]

gtex_eqtl[, Position := as.integer(Position)]
metaGWAS_data[, Position := as.integer(Position)]


# -----------------------------------------------
# Merge Step 1: Merge on Chr and Position only
# -----------------------------------------------

gtex_eqtl[, Chr := as.character(Chr)]
metaGWAS_data[, Chr := as.character(Chr)]
gtex_eqtl[, Position := as.integer(Position)]
metaGWAS_data[, Position := as.integer(Position)]

# Clean merge with only shared positions
merged_data_chr_pos <- merge(
  gtex_eqtl, 
  metaGWAS_data[, .(Chr, Position, MarkerName, Allele1, Allele2)], 
  by = c("Chr", "Position"), 
  all = FALSE  # INNER merge only
)



# -----------------------------------------------
# Step 2: Filter rows where alleles match between both datasets
# -----------------------------------------------
merged_data_filtered_alleles <- merged_data_chr_pos[
  (Allele1.x == Allele1.y & Allele2.x == Allele2.y) |
    (Allele1.x == Allele2.y & Allele2.x == Allele1.y)
]



# -----------------------------------------------
# Rename MarkerName to variant_id_metaGWAS
# -----------------------------------------------
setnames(merged_data_filtered_alleles, "MarkerName", "variant_id_metaGWAS")
setnames(merged_data_chr_pos, "MarkerName", "variant_id_metaGWAS")


# -----------------------------------------------
# Final output
# -----------------------------------------------
cat("Total shared SNPs on Chr:Pos:", nrow(merged_data_chr_pos), "\n")
cat("After allele matching:", nrow(merged_data_filtered_alleles), "\n")



# Optional: Save to file
# fwrite(merged_data_filtered, "Q:/Stagiairs/Bashir Hussein.W/Data/GTEx/merged_filtered_data.txt", sep = "\t")






#-----------------------------------------------------------------------------------------------------------
# laad de GTEx 
tar_file <- "Q:/Stagiairs/Bashir Hussein.W/Data/GTEx/GTEx_Analysis_v10_eQTL (1).tar"
untar(tar_file, exdir = "Q:/Stagiairs/Bashir Hussein.W/Data/GTEx/GTEx_Analysis_v10_eQTL_updated")

