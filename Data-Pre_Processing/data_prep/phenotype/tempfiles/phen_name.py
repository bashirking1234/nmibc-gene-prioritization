import pandas as pd

pheno = pd.read_csv('phenotypes_80_final.bed.gz', sep='\t', compression='gzip')
print(pheno.columns.tolist())

covar = pd.read_csv('covariates_tensorqtl_ready.txt', sep='\t', index_col=0)
print("Covariate columns:")
print(covar.columns.tolist())

