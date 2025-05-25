"""
Transforms a BED file containing phenotype expression data into a format compatible with tensorQTL.

This script performs the following steps:
1. Loads a BED file with all columns as strings to preserve formatting.
2. Splits the metadata (first 4 columns) from the expression data.
3. Cleans the expression data by removing thousand separators, replacing commas with dots, and converting to float.
4. Applies log2(TPM + 1) transformation to the expression values.
5. Combines the metadata and transformed expression data.
6. Writes the result to a new BED file ready for tensorQTL analysis.

Input:
    - phenotype_synced.bed.gz: Input BED file with phenotype data.

Output:
    - phenotype_tensorqtl_ready.bed: Output BED file with transformed expression values.

Note:
    - The output file should be compressed (e.g., using bgzip) before use with tensorQTL.
"""

import numpy as np
import pandas as pd


# Path to your existing BED file
input_path = "phenotype_synced.bed.gz"
output_path = "phenotype_tensorqtl_ready.bed"

# Step 1: Load the BED file with all columns as strings
df = pd.read_csv(input_path, sep="\t", dtype=str)

# Step 2: Split metadata and expression data
# Assuming the first 4 columns are chr,start,end & gene_ids and the rest are samples_ids with (normalised expression values)
meta = df.iloc[:, :4]
expr = df.iloc[:, 4:]

# Stap 3: remove duizendtallen, vervang komma's door punten en converteer naar float
expr = expr.replace(r"\.", "", regex=True)
expr = expr.replace(",", ".", regex=True)
expr = expr.applymap(lambda x: float(x) if pd.notna(x) and x != "" else np.nan)

# Step 4: log2(TPM + 1)
# Note: This assumes that the expression data is in TPM format.
# Log2 transformation makes the data more symmetric, so regression works better in TensorQTL.
expr_log = np.log2(expr + 1)

# Step 5: Combine metadata and transformed expression data
df_out = pd.concat([meta, expr_log], axis=1)

# Step 6: Write the result to a new BED file
# Note: The output file should be compressed (e.g., using bgzip) before use with tensorQTL.
df_out.to_csv(output_path, sep="\t", index=False)

