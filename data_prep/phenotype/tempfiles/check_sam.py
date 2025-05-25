import pandas as pd

# Load phenotype file
pheno = pd.read_csv("phenotypes_80_final.bed.gz", sep='\t', compression='gzip')
pheno_samples = pheno.columns[4:].astype(str).str.strip()  # exclude chr/start/end/gene_id

# Load covariates file (already transposed)
cov = pd.read_csv("covariates_tensorqtl_fixed.txt", sep='\t', index_col=0)
cov_samples = cov.columns.astype(str).str.strip()

# Compare
print("❌ Missing in covariates:", sorted(set(pheno_samples) - set(cov_samples)))
print("❌ Missing in phenotype:", sorted(set(cov_samples) - set(pheno_samples)))

