#!/bin/bash
#SBATCH --job-name=plink2_chr22
#SBATCH --output=/home/hbashir1/nmibc_eqtl_analysis/data_prep/genotype/logs/plink2.out
#SBATCH --error=/home/hbashir1/nmibc_eqtl_analysis/data_prep/genotype/logs/plink2.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

# Input files
VCF=/home/hbashir1/Data/SNPdata_UroLife_NBCS2/urolife_nbcs2_samen_chr22.dose.vcf.gz
OUTDIR=/home/hbashir1/nmibc_eqtl_analysis/data_prep/genotype/output_dir
KEEP=/home/hbashir1/nmibc_eqtl_analysis/data_prep/genotype/samples_to_keep.txt

# Run PLINK2
~/miniconda3/envs/myenv/bin/plink2 --vcf $VCF --output-chr chrM --rm-dup exclude-mismatch --maf 0.01 --keep $KEEP --make-pgen --out $OUTDIR/chr22c_data

