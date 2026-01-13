#!/bin/bash
#SBATCH --job-name=sourced_antibody_v1_predict
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/sourced_antibody_v1_predict/slurm.%j.log

export NXF_LOG_FILE=logs/sourced_antibody_v1_predict/.nextflow.log
export NXF_CACHE_DIR=logs/sourced_antibody_v1_predict/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/pipelines/sourced_antibody_v1/predict.nf \
    --input data/sourced_antibody_v1/assay_scaffold_summary.parquet \
    -output-dir data/sourced_antibody_v1 \
    -profile gemini \
    -with-report logs/sourced_antibody_v1_predict/reports/report.html \
    -with-timeline logs/sourced_antibody_v1_predict/reports/timeline.html