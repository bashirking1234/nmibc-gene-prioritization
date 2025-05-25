#!/bin/bash
#SBATCH --job-name=read_parquet
#SBATCH --output=read_parquet.log
#SBATCH --error=read_parquet.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# Load necessary modules (if any)
module load python/3.8  # Adjust the Python version as needed

# Run the Python script
python3 /home/hbashir1/metaGWASPipeline/rnaseq_eQTL/read_parquet.py 
