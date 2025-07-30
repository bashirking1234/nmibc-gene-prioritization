#!/bin/bash
#SBATCH --job-name=plink2_chr22
#SBATCH --output=/logs/plink2.out
#SBATCH --error=/logs/plink2.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

# Input files. full path
VCF=null
OUTDIR=null
KEEP=null

# Run PLINK2
~/miniconda3/envs/myenv/bin/plink2 --vcf $VCF --output-chr chrM --rm-dup exclude-mismatch --maf 0.01 --keep $KEEP --make-pgen --out $OUTDIR/chr22c_data

