## 

.libPaths("Q:\\Stagiairs\\Bashir Hussein.W\\R")

## loading the eQTL output from Tensorqtl, formatting it to the desired format & extracting gene positions

# load the eqtl output file from tensorqtl
output_tensorqtl <- "Q:/Stagiairs/Bashir Hussein.W/Data/OUTPUT_sampledata_eqtl/output_prefix.cis_qtl.txt/output_prefix.cis_qtl.txt"
output_tensorqtl_table <- fread(output_tensorqtl, sep = "\t",  header = TRUE)

# Transform the variant_id column from _ to :
output_tensorqtl_table[, variant_id := gsub("_", ":", variant_id)]
# Remove the ":b38" suffix from all entries in the variant_id column
output_tensorqtl_table[, variant_id := sub(":b38", "", variant_id)]

# Extract chromosome, position, and alleles from variant_id column
output_tensorqtl_table[, c("Chr", "Position", "Allele1", "Allele2") := tstrsplit(variant_id, ":", keep = c(1, 2, 3, 4))]
output_tensorqtl_table[, Chr := sub("chr", "", Chr)]
output_tensorqtl_table[, Position := as.integer(Position)]

## loading the metaGWAS output file & extracting gene positions

# Load the metaGWAS data
metaGWAS_file_path <- "Z:/genepi/Bashir_Hussein/metaGWAS_outputbestand/MA_STDERR_gwas_prognose_update_maf005_r208_recur_1.tbl"
metaGWAS_data <- fread(metaGWAS_file_path, sep = "\t", header = TRUE)

metaGWAS_data[, c("Chr", "Position", "Allele1", "Allele2") := tstrsplit(MarkerName, ":", keep = c(1, 2, 3, 4))]
metaGWAS_data[, Chr := sub("chr", "", Chr)]
metaGWAS_data[, Position := as.integer(Position)]


# Make sure both columns are of the same type
metaGWAS_data[, Chr := as.character(Chr)]
output_tensorqtl_table[, Chr := as.character(Chr)]

# Check for common values between the datasets based on Chr, Position, Allele1, and Allele2
common_values_tensorqtl_data <- merge(metaGWAS_data[, .(Chr, Position, Allele1, Allele2)], output_tensorqtl_table[, .(Chr, Position, Allele1, Allele2)], by = c("Chr", "Position", "Allele1", "Allele2"))
print("Common values between data and tcga_eqtl:")
print(common_values_tensorqtl_data)

# Merge the datasets based on Chr, Position, Allele1, and Allele2 columns
merged_data_tensorqtl <- merge(output_tensorqtl_table, metaGWAS_data[, .(MarkerName, Chr, Position, Allele1, Allele2)], by = c("Chr", "Position", "Allele1", "Allele2"), all.x = TRUE)

# Rename the MarkerName column to MarkerName_data to make sure there is no confusion.
setnames(merged_data_tensorqtl, "MarkerName", "MarkerName_data")

# Reorder columns to place MarkerName_data next to Variant
setcolorder(merged_data_tensorqtl, c("Cancer_type", "Gene_name", "RS_id", "Variant", "Chr", "MarkerName_data", "Position", "Allele1", "Allele2", "Credible_Set", "Credible_Set_Size", "PIP"))

# Filter out rows with NA in MarkerName_data
filtered_merged_data_tcga <- merged_data_tensorqtl[!is.na(MarkerName_data)]
