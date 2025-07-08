import pandas as pd

# Load phenotype file
pheno = pd.read_csv('phenotypes_80_final.bed.gz', sep='\t', compression='gzip')

# Get sample IDs from phenotype (columns 5+)
pheno_samples = list(pheno.columns[4:])

# Load covariates file (transposed already)
covar = pd.read_csv('covariates_tensorqtl_ready.txt', sep='\t', index_col=0)

# Get sample IDs from covariates (row index)
covar_samples = list(covar.index)

# Compare
missing_in_covar = [s for s in pheno_samples if s not in covar_samples]
missing_in_pheno = [s for s in covar_samples if s not in pheno_samples]

print("Missing in covariates:", missing_in_covar)
print("Missing in phenotype:", missing_in_pheno)

