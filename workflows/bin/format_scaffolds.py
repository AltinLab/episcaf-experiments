#!/usr/bin/env python
"""
Given a dataframe with at least the following columns:
- scaffolded_epitope_seq: sequence of residues of the rfd+mpnn predicted scaffold (including the epitope chunks)
- scaffolded_epitope_chunk_resindices: 0-indexed indices for residues in the scaffolded_epitope_seq which
  are part of the epitope chunks
- scaffolded_epitope_resindices: 0-indexed indices for residues in the scaffolded_epitope_seq which
  are part of the true epitope

Output a dataframe with all input columns as well as:
- assay_scaffolded_epitope_seq: scaffolded_epitope_seq with fixed N- and C-terminus attached
- assay_scaffolded_epitope_id: AF3-compatible identifier for the row computed from assay_scaffolded_epitope_seq
- assay_scaffolded_epitope_chunk_resindices: scaffolded_epitope_chunk_resindices adjusted for terminii offsets
- scaffolded_epitope_resindices scaffolded_epitope_resindices adjusted for terminii offsets
"""
import polars as pl
import argparse
import hashlib

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_pq",
        type=str,
    )
    parser.add_argument("--output_pq", type=str)
    parser.add_argument("--n_terminal", default="GA", type=str)
    parser.add_argument("--c_terminal", default="GDSLSWLLRLLNGNA", type=str)
    args = parser.parse_args()

    scaffolds = pl.read_parquet(args.input_pq)

    n_pad = len(args.n_terminal)

    assay_scaffolds = scaffolds.with_columns(
        (
            pl.lit(args.n_terminal)
            + pl.col("scaffolded_epitope_seq")
            + pl.lit(args.c_terminal)
        ).alias("assay_scaffolded_epitope_seq")
    ).with_columns(
        pl.col("assay_scaffolded_epitope_seq")
        .map_elements(
            lambda s: hashlib.md5(s.encode("utf-8")).hexdigest(), return_dtype=pl.String
        )
        .alias("assay_scaffolded_epitope_id"),
        pl.col("scaffolded_epitope_chunk_resindices")
        .list.eval(pl.element() + pl.lit(n_pad))
        .alias("assay_scaffolded_epitope_chunk_resindices"),
        pl.col("scaffolded_epitope_resindices")
        .list.eval(pl.element() + pl.lit(n_pad))
        .alias("assay_scaffolded_epitope_resindices"),
    )

    assay_scaffolds.write_parquet(args.output_pq)
