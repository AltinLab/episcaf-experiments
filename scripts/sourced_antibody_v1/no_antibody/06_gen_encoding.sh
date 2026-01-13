#!/bin/bash
#SBATCH --job-name=sourced_antibody_v1_no_antibody_simple_gen_encoding
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH --time=8:00:00
#SBATCH --output=logs/sourced_antibody_v1/no_antibody/gen_encoding/slurm.%j.log

conda run -n pepseq_encoding --live-stream python \
    /tgen_labs/altin/alphafold3/workspace/Library-Design/oligo_encoding/fasta_to_encodable.py \
        -i data/sourced_antibody_v1/no_antibody/sourced_antibody_v1_no_antibody_simple_403.fasta \
        -o data/sourced_antibody_v1/no_antibody/dp2/sourced_antibody_v1_no_antibody_simple_encodable.csv \
        -p "dp2_RFD" \
        -m data/sourced_antibody_v1/no_antibody/dp2/sourced_antibody_v1_no_antibody_simple_encodable_map.csv

cd /tgen_labs/altin/alphafold3/workspace/Library-Design/oligo_encoding

conda run -n pepseq_encoding --live-stream python \
    oligo_encoding.py \
        -m DeepLearning_model_R_1539970074840_1_20181019 \
        -r /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody/dp2/sourced_antibody_v1_no_antibody_simple_ratio_file.csv \
        -n 1 \
        -o /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody/dp2/sourced_antibody_v1_no_antibody_encoded.csv \
        -s /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody/dp2/sourced_antibody_v1_no_antibody_oligo_encoding_sequences.csv \
        -i /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody/dp2/sourced_antibody_v1_no_antibody_simple_encodable.csv \
        -p codon_weights.csv