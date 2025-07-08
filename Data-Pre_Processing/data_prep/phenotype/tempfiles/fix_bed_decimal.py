import pandas as pd

# Load the corrupted BED
df = pd.read_csv("phenotype_synced.bed.gz", sep="\t", dtype=str)

# Clean only the sample columns
meta = df.iloc[:, :4]
samples = df.iloc[:, 4:]

# Step 1: remove thousand separators ('.')
samples = samples.replace(r"\.", "", regex=True)

# Step 2: convert to float
samples = samples.applymap(lambda x: float(x.replace(",", ".")) if pd.notna(x) else x)

# Recombine and write cleaned BED
df_fixed = pd.concat([meta, samples], axis=1)
df_fixed.to_csv("phenotype_fixed.bed", sep="\t", index=False)

print(" Cleaned BED file written to: phenotype_fixed.bed")
print(" Compress it with: bgzip phenotype_fixed.bed")
