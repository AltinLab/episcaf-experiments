#!/usr/bin/env python

import MDAnalysis as mda
import argparse
import polars as pl
from pathlib import Path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--encoded_csv",
        type=str,
    )
    parser.add_argument(
        "--metadata_parquet",
        type=str,
    )
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument(
        "--output_annot",
        type=str,
    )
    parser.add_argument(
        "--output_encoded",
        type=str,
    )

    args = parser.parse_args()

    offset = args.offset

    encoded = pl.read_csv(args.encoded_csv)

    meta = pl.read_parquet(args.metadata_parquet)

    encoded = encoded.join(
        meta, left_on="AA Peptide", right_on="scaffolded_epitope_seq"
    )

    all = (
        encoded.with_row_index()
        .with_columns(
            (
                pl.concat_str(
                    pl.lit("DP2_"),
                    (pl.col("index") + 1 + offset).cast(pl.String).str.zfill(3),
                )
            ).alias("Name"),
            pl.concat_str(
                "id",
                pl.lit("_"),
                "contig_id",
                pl.lit("_"),
                "rfd_id",
                pl.lit("_"),
                "mpnn_id",
            ).alias("Design_ID"),
            pl.lit("RFD").alias("Model"),
            pl.col("id").alias("Target"),
            pl.lit(None).alias("Kd"),
            pl.lit(104).alias("BinderSequenceLength"),
            pl.lit(None).alias("toOrder"),
            pl.concat_str("heavy_seq", pl.lit(":"), "light_seq").alias("Sequence"),
        )
        .with_columns(pl.col("Name").alias("Seq ID"))
        .rename(
            {
                "AA Peptide": "BinderSequence",
            }
        )
    ).select(
        "Name",
        "Sequence",
        "Model",
        "Design_ID",
        "Target",
        "Kd",
        "BinderSequence",
        "BinderSequenceLength",
        "toOrder",
        "Seq ID",
        "Nucleotide Encoding w/ Adapters",
    )

    all.select(
        "Name",
        "Sequence",
        "Model",
        "Design_ID",
        "Target",
        "Kd",
        "BinderSequence",
        "BinderSequenceLength",
        "toOrder",
    ).write_csv(args.output_annot)

    all.select(
        "Seq ID",
        "Nucleotide Encoding w/ Adapters",
    ).write_csv(args.output_encoded)
