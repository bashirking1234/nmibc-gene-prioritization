.libPaths("Q:\\Stagiairs\\Bashir Hussein.W\\R")

install.packages("coloc")
library(coloc)

# https://cran.r-project.org/web/packages/coloc/vignettes/a03_enumeration.html

# Read the TSV files
eqtl_results <- read.table("Q:/Stagiairs/Bashir Hussein.W/Data/GTEx/filtered_merged_data_eqtl_gtex_results", header = TRUE, sep = "\t")
metaGWAS_data <- read.table("Z:/genepi/Bashir_Hussein/metaGWAS_outputbestand/MA_STDERR_gwas_prognose_update_maf005_r208_recur_1.tbl", header = TRUE, sep = "\t")


# Rename SNP columns
colnames(eqtl_results)[which(names(eqtl_results) == "variant_id")] <- "snp"
colnames(metaGWAS_data)[which(names(metaGWAS_data) == "MarkerName")] <- "snp"

# Merge datasets on SNP column
merged_data <- merge(eqtl_results, metaGWAS_data, by = "snp")
merged_data <- merged_data[!duplicated(merged_data$snp), ]

# Convert columns to numeric where needed
merged_data$beta_eqtl <- as.numeric(merged_data$slope)  # Effect size (regression slope)
merged_data$varbeta_eqtl <- as.numeric(merged_data$slope_se)^2  # Variance = SE^2 (slope_se)
merged_data$sample_size_eqtl <- as.numeric(merged_data$num_var)  # Assuming num_var corresponds to sample size, double-check

# Assign GWAS data columns (already correctly identified in your previous code)
merged_data$beta_gwas <- as.numeric(merged_data$Effect)  # Effect size from GWAS
merged_data$varbeta_gwas <- as.numeric(merged_data$StdErr)^2  # Squared SE for variance
merged_data$sample_size_gwas <- as.numeric(merged_data$TotalSampleSize)  # Sample size for GWAS

# Calculate number of cases and controls
merged_data$cases_gwas <- as.numeric(merged_data$TotalEvents)
merged_data$controls_gwas <- merged_data$sample_size_gwas - merged_data$cases_gwas


# Calculate global s from the full study
total_cases <- sum(merged_data$cases_gwas, na.rm = TRUE)
total_controls <- sum(merged_data$controls_gwas, na.rm = TRUE)
overall_s <- total_cases / (total_cases + total_controls)


# Calculate sdY (standard deviation of the trait)
sdY <- sd(merged_data$slope, na.rm = TRUE)

# Create lists for coloc analysis
eqtl_results <- list(
  snp = merged_data$snp,  # SNP column
  beta = merged_data$beta_eqtl,  # Effect size from eQTL
  varbeta = merged_data$varbeta_eqtl,  # Variance from eQTL
  N = merged_data$sample_size_eqtl,  # Sample size for eQTL
  type = "quant",
  sdY = sdY  # Standard deviation of the trait
)



# Then use this single value:
metaGWAS_data <- list(
  snp = merged_data$snp,
  beta = merged_data$beta_gwas,
  varbeta = merged_data$varbeta_gwas,
  N = merged_data$sample_size_gwas,
  s = overall_s,
  type = "cc"
)

# Run coloc analysis
result <- coloc.abf(eqtl_results, metaGWAS_data)

print(result)

# extract the more likely causal variants by
subset(result$results,SNP.PP.H4>0.01)

# Extract the 95% credible set
o <- order(result$results$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(result$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
result$results[o,][1:w,]$snp

