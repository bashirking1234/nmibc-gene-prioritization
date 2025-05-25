import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")

# Paths
plink_prefix_path = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/genotype/output_dir/chr22_data"
expression_bed = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/phenotypes_80_final_clean.bed.gz"
covariates_file = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/covariates_tensorqtl_aligned.txt"
prefix = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/tensorqtl_output"


# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
print(phenotype_pos_df['chr'].unique())
# Ensure all chromosome labels are strings and start with 'chr'
phenotype_pos_df['chr'] = phenotype_pos_df['chr'].astype(str)
phenotype_pos_df['chr'] = phenotype_pos_df['chr'].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)


covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df

# 1. Find common samples
common_samples = sorted(list(set(genotype_df.columns) & set(phenotype_df.columns) & set(covariates_df.index)))
print(f"Number of shared samples: {len(common_samples)}")

# 2. Subset and reorder all matrices to common samples
genotype_df = genotype_df[common_samples]
phenotype_df = phenotype_df[common_samples]
covariates_df = covariates_df.loc[common_samples]

# 3. Ensure numeric types
phenotype_df = phenotype_df.apply(pd.to_numeric, errors='coerce')
covariates_df = covariates_df.apply(pd.to_numeric, errors='coerce')

# 4. Check for missing values
print("Missing values in phenotype data:", phenotype_df.isnull().sum().sum())
print("Missing values in covariate data:", covariates_df.isnull().sum().sum())

print("Phenotype DataFrame shape:", phenotype_df.shape)
print(phenotype_df.head())
print(phenotype_df.isnull().sum())


# map all cis-associations (results for each chromosome are written to file)

# all genes
# cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df)


# genes on chr22
cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr'] == 'chr22'],
                phenotype_pos_df.loc[phenotype_pos_df['chr'] == 'chr22'],
                prefix, covariates_df=covariates_df)

# load results
pairs_df = pd.read_parquet(f'{prefix}.cis_qtl_pairs.chr22.parquet')
pairs_df.head()
