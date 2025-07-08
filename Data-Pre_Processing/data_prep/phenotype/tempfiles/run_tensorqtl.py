import pandas as pd
import torch
import tensorqtl
from tensorqtl import cis, pgen
import time
import psutil
import numpy as np

print(f" Memory available: {psutil.virtual_memory().available / 1024 ** 3:.2f} GB")

# Define file paths
plink_prefix_path = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/genotype/output_dir/chr22_data"
expression_bed = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/phenotype_tensorqtl_ready.bed.gz"
covariates_file = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/covariates_tensorqtl_aligned.txt"
prefix = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/tensorqtl_output"

# Load phenotype expression and positions
print(" Loading phenotypes...")
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)

# Ensure chromosome labels start with 'chr'
phenotype_pos_df['chr'] = phenotype_pos_df['chr'].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)
print("Chromosomes found in phenotype BED:", phenotype_pos_df['chr'].unique())

# Subset to chromosome 22 only
phenotype_df = phenotype_df.loc[phenotype_pos_df['chr'] == 'chr22']
phenotype_pos_df = phenotype_pos_df.loc[phenotype_pos_df['chr'] == 'chr22']
print(f" Filtered to {len(phenotype_df)} phenotypes on chr22")

# Clean sample names (strip whitespace)
original_sample_columns = phenotype_df.columns[4:]
phenotype_df.columns = list(phenotype_df.columns[:4]) + [s.strip() for s in original_sample_columns]

# Convert sample columns to numeric
sample_columns = phenotype_df.columns.difference(['chr', 'start', 'end', 'gene_id'])
phenotype_df[sample_columns] = phenotype_df[sample_columns].apply(pd.to_numeric, errors='coerce')

# Load covariates
print(" Loading covariates...")
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
print(" Covariates loaded.")

# Load genotype data
print(" Loading genotypes...")
start = time.time()
print("start of reading genotypes")
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
print("finished reading genotypes")
variant_df = pgr.variant_df
end = time.time()
print(f" Genotypes loaded. Time taken: {end - start:.2f} seconds")

# Ensure chromosome labels are strings
phenotype_pos_df['chr'] = phenotype_pos_df['chr'].astype(str)
variant_df['chrom'] = variant_df['chrom'].astype(str)

# Debug info
print(" Phenotype sample value dtypes:")
print(phenotype_df.iloc[:, 4:].dtypes)

# Collect sample IDs
geno_samples = list(genotype_df.columns)
pheno_samples = list(phenotype_df.columns)
cov_samples = list(covariates_df.index)

print("Sample counts:")
print(f"- Genotype samples: {len(geno_samples)}")
print(f"- Phenotype samples: {len(pheno_samples)}")
print(f"- Covariate samples: {len(cov_samples)}")

# Determine shared samples
shared_samples = list(set(geno_samples) & set(pheno_samples) & set(cov_samples))
print(f"Number of shared samples: {len(shared_samples)}")

if geno_samples != pheno_samples:
    print("Warning: genotype and phenotype sample order do not match")
if geno_samples != cov_samples:
    print("Warning: genotype and covariate sample order do not match")

# Reorder to common sample list
common = sorted(shared_samples)
genotype_df = genotype_df[common]
phenotype_df = phenotype_df[common]
covariates_df = covariates_df.loc[common]

# Check for missing values
print("Missing values in phenotype data:", phenotype_df.isnull().sum().sum())
print("Missing values in covariate data:", covariates_df.isnull().sum().sum())

# ============================
# 🔍 VALIDATION STEP
# ============================
print(" Validating phenotype matrix before regression...")

min_samples = 10   # Minimum valid (non-NaN) samples per gene
min_std = 1e-5     # Minimum variation required

non_na_counts = phenotype_df.notna().sum(axis=1)
stds = phenotype_df.std(axis=1)
valid_genes = (non_na_counts >= min_samples) & (stds > min_std)
n_valid = valid_genes.sum()

print("\n📊 QC Summary for phenotypes:")
print(f"- Total genes: {len(phenotype_df)}")
print(f"- Genes with ≥10 valid samples: {(non_na_counts >= min_samples).sum()}")
print(f"- Genes with std > {min_std}: {(stds > min_std).sum()}")


if n_valid == 0:
    raise ValueError("❌ ERROR: No usable phenotypes (genes) for regression. All rows have too few samples or no variance. Aborting.")
else:
    print(f"✅ {n_valid} valid phenotypes detected. Proceeding with TensorQTL...")

# Filter phenotypes and positions
phenotype_df = phenotype_df.loc[valid_genes]
phenotype_pos_df = phenotype_pos_df.loc[valid_genes]

# ============================
# Run TensorQTL
# ============================
print(" Running TensorQTL (cis)...")
start = time.time()
cis.map_nominal(
    genotype_df, variant_df,
    phenotype_df, phenotype_pos_df,
    prefix, covariates_df=covariates_df
)
end = time.time()
print(f" TensorQTL run completed. Time taken: {end - start:.2f} seconds")
