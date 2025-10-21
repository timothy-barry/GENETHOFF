#!/bin/bash
#SBATCH --job-name=process_spry_1620_data
#SBATCH --output=logs/process_spry_1620_data_%j.out
#SBATCH --error=logs/process_spry_1620_data_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH -c 6
#SBATCH --partition short

# load modules
mamba init
mamba activate crisprde-venv

# run command
snakemake -s ../00-pipeline/genethOFF.snakemake \
  -j 6 \
  -k \
  --use-conda \ 
  --report-after-run --report runtime_report.html  \
  -n