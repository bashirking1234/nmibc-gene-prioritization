import pandas as pd

"""
This script processes a covariates file for tensorQTL analysis.
It performs the following steps:
1. Loads the covariates file.
2. Checks the orientation of the matrix (samples in rows or columns).
3. Filters samples based on a provided list.
4. Drops samples with missing values.
5. Replaces the index with PC labels (PC1, PC2, ...).
6. Ensures all values are floats.
7. Inserts covariate names as the first column and adds a blank top-left cell.
8. Saves the final cleaned covariates file as a tab-separated file without an index.
"""

def main():
    # File paths
    covariates_path = "data_prep/covariates/covariates_4PCAs.txt"
    keep_samples_path = "data_prep/covariates/samples_to_keep.txt" # Optional change to None if not needed
    output_path = "data_prep/covariates/covariates_tensorqtl_ready2.txt"

    print("Step 1: Loading covariates file")
    cov = pd.read_csv(covariates_path, sep="\t", index_col=0, decimal=",")


    if cov.shape[0] > cov.shape[1]:
        print("Transposing matrix: samples appear to be in rows")
        cov = cov.transpose()
    else:
        print("Matrix orientation looks correct (samples as columns)")

    print(f"Matrix dimensions: {cov.shape[0]} covariates x {cov.shape[1]} samples")

    # Filter sample IDs
    if keep_samples_path is not None:
        print(f"Filtering samples using list: {keep_samples_path}")
        with open(keep_samples_path, "r") as f:
            keep_ids = [line.strip() for line in f if line.strip()]
        cov.columns = cov.columns.astype(str).str.strip()
        keep_ids = [s.strip() for s in keep_ids]
        cov = cov[[s for s in cov.columns if s in keep_ids]]
        print(f"Samples retained after filtering: {len(cov.columns)}")

    # Drop columns (samples) with missing values
    cov_clean = cov.dropna(axis=1, how="any")
    removed = cov.shape[1] - cov_clean.shape[1]
    if removed > 0:
        print(f"Removed {removed} samples with missing values (NA)")
    else:
        print("No missing values detected")

    # Replace index with PC labels (PC1, PC2, ...)
    cov_clean.index = [f"PC{i+1}" for i in range(cov_clean.shape[0])]

    # Make sure all values are floats
    cov_clean = cov_clean.astype(float)

    # Insert covariate names as first column and add blank top-left cell
    cov_final = cov_clean.copy()
    cov_final.insert(0, "", cov_final.index)

    # Save as tab-separated file without index
    cov_final.to_csv(output_path, sep="\t", index=False)

    print(f"Final cleaned covariates file saved to: {output_path}")

if __name__ == "__main__":
    main()
    print("Script completed successfully.")
    print("Note: Check the output file for any discrepancies.")
    print("Note: Ensure that the covariates are in the correct format for tensorQTL.")
    print("Note: The covariates file is now ready for tensorQTL analysis.")
    print("Note: If you have any questions, please refer to the documentation.")
    print("Note: Thank you for using this script.")

   