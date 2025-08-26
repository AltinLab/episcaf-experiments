#!/bin/bash
#SBATCH --job-name=sourced_antibody_scaffold
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=8:00:00
#SBATCH --output=logs/sourced_antibody_scaffold/slurm.%j.log

export NXF_LOG_FILE=logs/sourced_antibody_scaffold/.nextflow.log
export NXF_CACHE_DIR=logs/sourced_antibody_scaffold/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/sourced_v1/scaffold.nf \
    --input data/abdb/complexes/chosen_complex.parquet \
    --pdb_dir data/abdb/complex_pdbfiles/cleaned \
    -output-dir data/abdb/scaffold \
    -profile gemini \
    -with-report logs/sourced_antibody_scaffold/reports/report.html \
    -with-timeline logs/sourced_antibody_scaffold/reports/timeline.html