#!/bin/bash
#SBATCH --job-name=sourced_antibody_v1_scaffold
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH -c 1
#SBATCH --time=3:00:00
#SBATCH --output=logs/sourced_antibody_v1_scaffold/slurm.%j.log

export NXF_LOG_FILE=logs/sourced_antibody_v1_scaffold/.nextflow.log
export NXF_CACHE_DIR=logs/sourced_antibody_v1_scaffold/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/pipelines/sourced_antibody_v1/scaffold.nf \
    --input data/Xsourced_antibody_v1/full_antibody/Xabdb_contigs.parquet \
    --pdb_dir data/abdb/complex_pdbfiles/cleaned \
    --rfd_num_designs 10 \
    --mpnn_seqs_per_struct 8 \
    -output-dir data/Xsourced_antibody_v1/full_antibody \
    -profile gemini \
    -with-report logs/sourced_antibody_v1_scaffold/reports/report.html \
    -with-timeline logs/sourced_antibody_v1_scaffold/reports/timeline.html \
    -resume