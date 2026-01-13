#!/bin/bash
#SBATCH --job-name=sourced_antibody_v1_scaffold
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/sourced_antibody_v1_scaffold/slurm.%j.log

export NXF_LOG_FILE=logs/sourced_antibody_v1_scaffold/.nextflow.log
export NXF_CACHE_DIR=logs/sourced_antibody_v1_scaffold/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/pipelines/sourced_antibody_v1/scaffold.nf \
    --input data/sourced_antibody_v1/abdb_contigs.parquet \
    --pdb_dir data/abdb/complex_pdbfiles/cleaned \
    --rfd_num_designs 8 \
    --mpnn_seqs_per_struct 8 \
    -output-dir data/sourced_antibody_v1 \
    -profile gemini \
    -with-report logs/sourced_antibody_v1_scaffold/reports/report.html \
    -with-timeline logs/sourced_antibody_v1_scaffold/reports/timeline.html