
"""
This script splits a BED-format gene expression file into two separate files based on sample sex.

Inputs:
    - bed_file: Path to the input BED-format expression file (tab-separated, gzipped).
    - man_ids_file: Path to a text file containing male sample IDs (one per line).
    - vrouw_ids_file: Path to a text file containing female sample IDs (one per line).

Outputs:
    - phenotype_chrX_male.bed.gz: BED-format file containing only male samples.
    - phenotype_chrX_female.bed.gz: BED-format file containing only female samples.
"""
import pandas as pd


# 1. Loads lists of male and female sample IDs from text files.
bed_file = ""
male_ids_file = "male_samples.txt"
female_ids_file ="female_samples.txt"

# 2. Reads a BED-format gene expression file (with gene info in the first 4 columns and sample IDs as column headers).
with open(male_ids_file) as f:
    male_ids = set(line.strip() for line in f if line.strip())
with open(female_ids_file) as f:
    female_ids = set(line.strip() for line in f if line.strip())

# read expression data
df = pd.read_csv(bed_file, sep="\t")

# Geninfo = First 4 columns
geninfo = df.iloc[:, :4]
expr = df.iloc[:, 4:]

# 3. Filters the expression data to retain only columns (samples) corresponding to males and females, respectively.
expr_male = expr.loc[:, expr.columns.isin(male_ids)]
expr_female = expr.loc[:, expr.columns.isin(female_ids)]

# 4. Combines the filtered expression data with the gene information columns.
df_male = pd.concat([geninfo, expr_male], axis=1)
df_female = pd.concat([geninfo, expr_female], axis=1)

# 5. Writes the resulting male and female-specific BED files to disk, compressed with gzip.
df_male.to_csv("phenotype_chrX_male.bed", sep="\t", index=False, compression="gzip")
df_female.to_csv("phenotype_chrX_female.bed", sep="\t", index=False, compression="gzip")
