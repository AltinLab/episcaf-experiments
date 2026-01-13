#!/bin/bash
#SBATCH --job-name=sourced_antibody_v1_gen_contig
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=8:00:00
#SBATCH --output=logs/sourced_antibody_v1/no_antibody/gen_contig/slurm.%j.log

conda run -n episcaf-experiments --live-stream python \
    workflows/bin/generate_contigs.py \
    --topology_path data/abdb/complex_pdbfiles/cleaned \
    --input_pq data/abdb/complexes/chosen_complex.parquet \
    --contig_per_pdb 100 \
    --target_length 104 \
    --no_antibody \
    --output_pq data/sourced_antibody_v1/no_antibody/abdb_contigs.parquet
