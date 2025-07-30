import pandas as pd
import os

"""
This script prepares eQTL (expression Quantitative Trait Loci) data for EasyStrata analysis by performing the following steps for each specified genomic region ('par' and 'nonpar'):

1. Reads male and female eQTL results from Parquet files.
2. Ensures 'phenotype_id' is present as a column in both datasets.
3. Merges male and female datasets on 'phenotype_id' and 'variant_id'.
4. Extracts and renames relevant columns for EasyStrata input.
5. Writes the combined data to a space-separated text file for EasyStrata.
6. Generates an EasyStrata script to perform P-difference and joint tests between sexes.
7. Creates a SLURM batch script to run the EasyStrata analysis on a computing cluster.

Directories for output, scripts, and SLURM jobs are created if they do not exist.

Dependencies:
    - pandas
    - os

Inputs:
    - Parquet files containing eQTL results for males and females for each region.

Outputs:
    - EasyStrata input files (space-separated .txt)
    - EasyStrata script files (.txt)
    - SLURM batch scripts (.sh)

Usage:
    Run this script directly to generate all necessary files for EasyStrata analysis and SLURM job submission.
import os

"""  

# === Config ===
regions = ['par', 'nonpar']
base_path = ""
input_suffix = ".cis_qtl_pairs.chrX.parquet"
output_dir = ""
script_dir = ""
slurm_dir = ""

# make directories if needed
os.makedirs(output_dir, exist_ok=True)
os.makedirs(script_dir, exist_ok=True)
os.makedirs(slurm_dir, exist_ok=True)

for region in regions:
    print(f"🔄 Processing region: {region.upper()}")

    # === read files
    male_file = os.path.join(base_path, f"tensorqtl_result_{region}_male{input_suffix}")
    female_file = os.path.join(base_path, f"tensorqtl_result_{region}_female{input_suffix}")

    male = pd.read_parquet(male_file)
    female = pd.read_parquet(female_file)

    # make sure that phenotype_id is a column
    if 'phenotype_id' not in male.columns:
        male = male.reset_index()
    if 'phenotype_id' not in female.columns:
        female = female.reset_index()

    # === Merge op phenotype_id + variant_id
    merged = male.merge(female, on=['phenotype_id', 'variant_id'], suffixes=('.men', '.women'))

    # === EasyStrata inputfiles generation

    easy_input_path = os.path.join(output_dir, f"chrX_{region}_eQTL_combined.txt")

    easy_df = merged[['variant_id', 'slope.men', 'slope_se.men', 'slope.women', 'slope_se.women']]
    easy_df = easy_df.rename(columns={
    'variant_id': 'RSID',
    'slope.men': 'beta.men',
    'slope_se.men': 'se.men',
    'slope.women': 'beta.women',
    'slope_se.women': 'se.women'})

    easy_df.to_csv(easy_input_path, sep=' ', index=False)

    print(f"EasyStrata input written: {easy_input_path}")

    # === EasyStrata script genereren
    easy_script_content = f"""
####################################################

DEFINE  --acolIn RSID;beta.men;se.men;beta.women;se.women
        --acolInClasses character;numeric;numeric;numeric;numeric
        --pathOut {output_dir}
        --strSeparator SPACE

EASYIN  --fileIn {easy_input_path}

START EASYSTRATA

CALCPDIFF   --acolBETAs beta.men;beta.women
            --acolSEs se.men;se.women
            --colOutPdiff p_diff_sex

JOINTTEST   --acolBETAs beta.men;beta.women
            --acolSEs se.men;se.women
            --colOutPjoint p_joint_sex

WRITE   --strMode txt
        --strPrefix chrX_{region}_eQTL_EasyStrata.
        --strSep TAB

STOP EASYSTRATA

####################################################
"""
    easy_script_path = os.path.join(script_dir, f"run_chrX_{region}_eQTL_easystrata.txt")
    with open(easy_script_path, "w") as f:
        f.write(easy_script_content.strip())
    print(f"EasyStrata script written: {easy_script_path}")

    # === generate SLURM-script
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=EasyStrata_{region}
#SBATCH --output={slurm_dir}/chrX_{region}_easystrata.out
#SBATCH --error={slurm_dir}/chrX_{region}_easystrata.err
#SBATCH --time=00:10:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

module load EasyStrata

EasyStrata {easy_script_path}
"""

    slurm_path = os.path.join(slurm_dir, f"run_chrX_{region}_easystrata.sh")
    with open(slurm_path, "w") as f:
        f.write(slurm_script.strip())
    print(f"SLURM script written: {slurm_path}\n")
