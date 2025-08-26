#!/usr/bin/env python
""" """
from mdaf3.FeatureExtraction import serial_apply

import MDAnalysis as mda
import numpy as np
import polars as pl
import argparse
from pathlib import Path


def annot_epitopes(row, top_path):

    u = mda.Universe(top_path / (row["complexed_antibody_id"] + ".pdb"))

    epitope_res = u.select_atoms(
        "segid A and (not name H*) and around 4 ((segid B or segid C) and (not name H*))"
    ).residues

    ep_res_idx = epitope_res.resindices
    row["epitope_resindices"] = list(ep_res_idx)
    row["epitope_seq"] = epitope_res.sequence(format="string")

    bm = np.zeros(len(u.residues), dtype=bool)
    bm[epitope_res.resindices] = True
    row["epitope_boolmask"] = list(bm)

    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--complex_pq",
        type=str,
    )
    parser.add_argument("--topology_path", type=str)
    parser.add_argument("--output_pq", type=str)
    args = parser.parse_args()

    complexes = pl.read_parquet(args.complex_pq)
    top_path = Path(args.topology_path)

    complexes = serial_apply(complexes, annot_epitopes, top_path)

    # now that we know the epitope residues, ensure uniqueness in the dataset:
    # if two rows have the same antigen, heavy seq, light seq, AND epitope residues,
    # keep the one with the better resolution
    exclude = ["antigen_seq", "heavy_seq", "light_seq", "epitope_seq"]
    columns = complexes.columns
    complexes = (
        complexes.group_by(exclude)
        .agg([pl.col(colname) for colname in columns if colname not in exclude])
        .with_columns(pl.col("resolution").list.arg_min().alias("keep_idx"))
        .with_columns(
            [
                pl.col(colname).list.get(pl.col("keep_idx")).alias(colname)
                for colname in columns
                if colname not in exclude
            ]
        )
    ).select(pl.exclude("keep_idx"))

    complexes.write_parquet(args.output_pq)
