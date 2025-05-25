import pandas as pd

# Load clean phenotype
pheno = pd.read_csv("phenotypes_80_final_clean.bed.gz", sep='\t', compression='gzip')
pheno_samples = pheno.columns[4:].astype(str).str.strip()

# Load covariates (transposed version)
cov = pd.read_csv("covariates_tensorqtl_fixed.txt", sep='\t', index_col=0)
cov.columns = cov.columns.astype(str).str.strip()

# Reorder covariates to match phenotype sample order
cov = cov[pheno_samples]

# Save aligned covariates
cov.to_csv("covariates_tensorqtl_aligned.txt", sep='\t')

