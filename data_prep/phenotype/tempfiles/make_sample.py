import pandas as pd

# Laad bestand zonder index (sample-ID’s staan in kolommen)
cov = pd.read_csv("covariates_tensorqtl_aligned.txt", sep="\t", index_col=0)

# Transponeer zodat sample-ID’s in index komen
cov = cov.transpose()

# Opslaan als correcte covariatenmatrix
cov.index.name = "IID"
cov.to_csv("covariates_transposed.txt", sep="\t")

# Schrijf PLINK2 sample keep-lijst (1 kolom: IID)
with open("samples_to_keep.txt", "w") as f:
    for s in cov.index:
        f.write(f"{s}\n")


