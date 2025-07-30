import pandas as pd
import os
import pyarrow.parquet as pq
import glob
""" 
This script combines cis-QTL result files for chromosomes 1–22 into a single Parquet file.

Workflow:
- Searches a specified directory for Parquet files matching the pattern for chromosomes 1–22.
- Excludes any files related to chromosome X.
- Reads each file into a pandas DataFrame.
- Concatenates all DataFrames into one.
- Saves the combined DataFrame as a new Parquet file.

Dependencies:
    - pandas
    - os
    - pyarrow
    - glob
    - pyarrow.parquet
Args:
    None

Returns:
    None

Output:
    A single Parquet file containing the combined results for chromosomes 1–22.

"""


# === Path to tensrorwtl outputs ===
base_dir = ""

# === collect all chr1–22 files ===
base_dir = ""
parquet_files = sorted(
    glob.glob(os.path.join(base_dir, "tensorqtl_result.cis_qtl_pairs.chr*.parquet"))
)
parquet_files = [f for f in parquet_files if "chrX" not in f]

print(f"📂 found {len(parquet_files)} files of (chr1–22)...")

# === combine the files ===
dfs = []
for f in parquet_files:
    df = pd.read_parquet(f)
    dfs.append(df)

combined_df = pd.concat(dfs, ignore_index=True)

# === save ===
out_path = os.path.join(base_dir, "tensorqtl_results_chr1_22_combined.parquet")
combined_df.to_parquet(out_path)
print(f"results saved: {out_path}")
