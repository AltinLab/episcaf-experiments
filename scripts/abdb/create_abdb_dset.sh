#!/bin/bash
#SBATCH --job-name=create_abdb_dset
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=8:00:00
#SBATCH --output=logs/abdb/slurm.%j.log

export NXF_LOG_FILE=logs/abdb/.nextflow.log
export NXF_CACHE_DIR=logs/abdb/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/pipelines/abdb/create_dataset.nf \
    --input data/abdb/raw/abdb_20240706 \
    -output-dir data/abdb \
    -profile gemini \
    -resume