#!/bin/bash
#SBATCH --job-name=tensorqtl_run
#SBATCH --output=/tensorqtl_run.out
#SBATCH --error=/tensorqtl_run.err
#SBATCH --time=00:20:00          
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=rome


# Run your TensorQTL Python script
python3 tensorqtl_script.py




