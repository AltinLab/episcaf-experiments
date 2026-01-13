#!/bin/bash

conda run -n epibind --live-stream python \
    ./workflows/bin/dp2_fmt.py \
        --encoded_csv data/sourced_antibody_v1/no_antibody/dp2/sourced_antibody_v1_no_antibody_encoded.csv \
        --metadata_parquet data/sourced_antibody_v1/no_antibody/assay_scaffold_simple_metrics_403.parquet \
        --offset 597 \
        --output_encoded data/sourced_antibody_v1/no_antibody/dp2/dp2_rfd_encoded.csv \
        --output_annot data/sourced_antibody_v1/no_antibody/dp2/dp2_rfd_annot.csv