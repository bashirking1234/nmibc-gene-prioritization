import pandas as pd

# Load your current covariates
cov = pd.read_csv("covariates_tensorqtl_ready.txt", sep="\t", index_col=0)

# Transpose to make sample IDs the column headers
cov_t = cov.transpose()

# Save in TensorQTL-compatible format
cov_t.to_csv("covariates_tensorqtl_fixed.txt", sep="\t")

