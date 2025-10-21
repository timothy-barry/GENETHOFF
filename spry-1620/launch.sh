#!/bin/bash
#SBATCH --job-name=process_spry_1620_data
#SBATCH --output=logs/process_spry_1620_data_%j.out
#SBATCH --error=logs/process_spry_1620_data_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH -c 12
#SBATCH --partition short

# source research config
pipeline_dir=$HOME"research_code/GENETHOFF/"
work_dir=$HOME"research_offsite/work"
cp $pipeline_dir"spry-1620/configuration.yml" $work_dir

# load modules
module load conda/miniforge3/24.11.3-0
conda init
conda activate crisprde-venv

# run command
snakemake -n \
-s $pipeline_dir"00-pipeline/genethOFF.snakemake" \
-j 12 \
--use-conda \
--directory $workdir
-n
