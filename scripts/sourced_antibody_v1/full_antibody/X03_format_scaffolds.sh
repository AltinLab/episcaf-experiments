#!/bin/bash
#SBATCH --job-name=sourced_antibody_v1_format_scaffolds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --output=logs/sourced_antibody_v1_format_scaffolds/slurm.%j.log

conda run -n episcaf-experiments --live-stream python \
    ./workflows/bin/format_scaffolds.py \
    --input_pq data/Xsourced_antibody_v1/full_antibody/scaffold_summary.parquet \
    --output data/Xsourced_antibody_v1/full_antibody/assay_scaffold.parquet
