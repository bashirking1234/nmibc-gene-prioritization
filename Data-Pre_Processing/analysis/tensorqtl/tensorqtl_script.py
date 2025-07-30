"""
This script performs cis-eQTL mapping using TensorQTL on genotype, phenotype, and covariate data.

Workflow:
1. Loads phenotype (expression) data from a BED file and filters for a chromosome.
2. Loads covariate data from a tab-delimited text file.
3. Loads genotype data from PLINK2 .pgen files for chromosome of your choice.
4. Harmonizes samples across genotype, phenotype, and covariate datasets to ensure only shared samples are analyzed.
5. Validates the presence of missing values in phenotype and covariate data.
6. Runs TensorQTL's cis-eQTL mapping using the harmonized data and saves results to the specified prefix.

Inputs:
- plink_prefix_path: Path prefix for PLINK2 .pgen, .pvar, and .psam files.
- expression_bed: Path to gzipped BED file containing phenotype (expression) data.
- covariates_file: Path to tab-delimited file containing covariate data.
- prefix: Output directory and file prefix for TensorQTL results.

Dependencies:
- pandas
- numpy
- tensorqtl
- psutil
- time

Note:
- The script assumes all input files are properly formatted and available at the specified paths. 
"""

# === Import libraries ===
import pandas as pd
import numpy as np
import tensorqtl
from tensorqtl import cis, pgen
import time
import psutil

print(f" Memory available: {psutil.virtual_memory().available / 1024 ** 3:.2f} GB")

# === Input paths ===
plink_prefix_path = ""
expression_bed = ""
covariates_file = ""


# === Output path ===
# Note: The prefix should not include the file extension. TensorQTL will add it automatically.
# Example of output will be: tensorqtl_result.cis_eqtl_pairs.chr22.parquet
# the path is only: 
# Note: Adjust the prefix to your desired output directory and file name.
prefix = "" 


# === Load expression data ===
print(" Loading phenotypes...")
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
phenotype_pos_df['chr'] = phenotype_pos_df['chr'].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)
phenotype_df.columns = list(phenotype_df.columns[:4]) + [s.strip() for s in phenotype_df.columns[4:]]
phenotype_df = phenotype_df.loc[phenotype_pos_df['chr'] == 'chr22']
phenotype_pos_df = phenotype_pos_df.loc[phenotype_pos_df['chr'] == 'chr22']

# Convert expression values to numeric
sample_columns = phenotype_df.columns.difference(['chr', 'start', 'end', 'gene_id'])
phenotype_df[sample_columns] = phenotype_df[sample_columns].apply(pd.to_numeric, errors='coerce')

# === Load covariates ===
print(" Loading covariates...")
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# === Load genotypes ===
print(" Loading genotypes...")
start = time.time()
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df
print(f" Genotypes loaded. Time taken: {time.time() - start:.2f} seconds")

# === Harmonize samples ===
geno_samples = list(genotype_df.columns)
pheno_samples = list(phenotype_df.columns)
cov_samples = list(covariates_df.index)
shared_samples = list(set(geno_samples) & set(pheno_samples) & set(cov_samples))

common = sorted(shared_samples)
genotype_df = genotype_df[common]
phenotype_df = phenotype_df[common]
covariates_df = covariates_df.loc[common]

# === Validate data ===
print("Missing values in phenotype data:", phenotype_df.isnull().sum().sum())
print("Missing values in covariate data:", covariates_df.isnull().sum().sum())


# === Run TensorQTL ===
print(" Running TensorQTL (cis)...")
start = time.time()
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df)
print(f" TensorQTL run completed. Time taken: {time.time() - start:.2f} seconds")
