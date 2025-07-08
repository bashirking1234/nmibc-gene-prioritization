import pandas as pd

# Load phenotype file
pheno = pd.read_csv("phenotypes_80_final.bed.gz", sep="\t", compression="gzip")

# Drop all unnamed columns (like 'Unnamed: 81')
pheno = pheno.loc[:, ~pheno.columns.str.contains("^Unnamed")]

# Save clean version
pheno.to_csv("phenotypes_80_final_clean.bed.gz", sep="\t", compression="gzip", index=False)

