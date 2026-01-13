#!/bin/bash
#SBATCH --job-name=sourced_antibody_v1_predict
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/sourced_antibody_v1/no_antibody/predict/slurm.%j.log

export NXF_LOG_FILE=logs/sourced_antibody_v1/no_antibody/predict/.nextflow.log
export NXF_CACHE_DIR=logs/sourced_antibody_v1/no_antibody/predict/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/pipelines/sourced_antibody_v1/predict.nf \
    --input data/sourced_antibody_v1/no_antibody/assay_scaffold.parquet \
    --inf_dir data/sourced_antibody_v1/no_antibody/af3_predictions \
    -output-dir data/sourced_antibody_v1/no_antibody \
    -profile gemini \
    -with-report logs/sourced_antibody_v1/no_antibody/predict/reports/report.html \
    -with-timeline logs/sourced_antibody_v1/no_antibody/predict/reports/timeline. \
    -resume