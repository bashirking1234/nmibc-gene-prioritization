import pandas as pd
import pyarrow.parquet as pq

"""
Reads a Parquet file containing cis-QTL pair results into a Pandas DataFrame and displays the first few rows and column names.

Steps:
1. Imports necessary libraries: pandas and pyarrow.parquet.
2. Specifies the path to the Parquet results file.
3. Reads the Parquet file into a Pandas DataFrame.
4. Prints the first few rows of the DataFrame.
5. Prints the column names of the DataFrame.

dependencies:
    - pandas
    - pyarrow
    - pyarrow.parquet
Args:
    None
Returns:
    None
Output:
    Displays the first few rows and column names of the DataFrame read from the Parquet file.       
"""

# Pad naar je resultaatbestand
parquet_file = "tensorqtl_results.cis_qtl_pairs.chr22.parquet"

# Lees het bestand in als Pandas DataFrame
df = pq.read_table(parquet_file).to_pandas()

# Toon de eerste paar rijen
print(df.head())

# Bekijk kolomnamen
print(df.columns)
