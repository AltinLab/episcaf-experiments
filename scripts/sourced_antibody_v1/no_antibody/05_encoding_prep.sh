#!/bin/bash
#SBATCH --job-name=sourced_antibody_v1_no_antibody_simple
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH --time=8:00:00
#SBATCH --output=logs/sourced_antibody_v1/no_antibody/extract_metrics/slurm.%j.log

conda run -n episcaf-experiments --live-stream python \
    workflows/bin/extract_metrics.py \
    --input_pq data/sourced_antibody_v1/no_antibody/assay_scaffold_simple.parquet \
    --true_dir data/abdb/complex_pdbfiles/cleaned \
    --mpnn_dir data/sourced_antibody_v1/no_antibody/proteinmpnn \
    --af3_dir data/sourced_antibody_v1/no_antibody/af3_predictions \
    --output_pq data/sourced_antibody_v1/no_antibody/assay_scaffold_simple_metrics.parquet \
    --output_fasta data/sourced_antibody_v1/no_antibody/dp2/sourced_antibody_v1_no_antibody_simple.fasta 


# additional step performed in notebooks/sourced_antibody_v1/sel_simple.ipynb