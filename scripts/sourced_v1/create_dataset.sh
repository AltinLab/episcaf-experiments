#!/bin/bash
#SBATCH --job-name=sourced_antibody
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=8:00:00
#SBATCH --output=logs/sourced_antibody/slurm.%j.log

export NXF_LOG_FILE=logs/sourced_antibody/.nextflow.log
export NXF_CACHE_DIR=logs/sourced_antibody/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/sourced_v1/create_dataset.nf \
    --input data/abdb/raw/abdb_20240706 \
    -output-dir data/abdb \
    -profile gemini