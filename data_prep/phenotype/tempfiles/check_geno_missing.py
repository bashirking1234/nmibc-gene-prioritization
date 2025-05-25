import pandas as pd
from tensorqtl import pgen

# === Invoerpad naar PLINK2 .pgen-bestand (zonder extensie) ===
plink_prefix = plink_prefix = plink_prefix = "/home/hbashir1/metaGWASPipeline/rnaseq_eQTL/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18"


# === Genotype data inladen ===
print("📥 Loading genotype data...")
pgr = pgen.PgenReader(plink_prefix)
genotype_df = pgr.load_genotypes()
print("✅ Genotypes loaded.")

# === Totale missingness berekenen ===
total_missing = genotype_df.isna().mean().mean()
print(f"\n➡️ Mean genotype missingness (alle samples × SNPs): {total_missing:.6f}")

# === Verdeling per variant (SNP) ===
print("\n📊 Missingness per SNP (variant):")
print(genotype_df.isna().mean(axis=1).describe())

# === Verdeling per sample ===
print("\n👤 Missingness per sample:")
print(genotype_df.isna().mean(axis=0).describe())
