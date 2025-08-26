#!/usr/bin/env python
import argparse
import polars as pl


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--complex_pq", type=str)
    parser.add_argument("--antigen_rep_tsv", type=str)
    parser.add_argument("--heavy_rep_tsv", type=str)
    parser.add_argument("--light_rep_tsv", type=str)
    parser.add_argument("--output_chosen", type=str)
    parser.add_argument("--output_discard", type=str)
    args = parser.parse_args()

    orig_pq = pl.read_parquet(args.complex_pq)

    antigen_df = (
        pl.read_csv(
            args.antigen_rep_tsv,
            separator="\t",
            has_header=False,
            new_columns=["antigen_rep", "antigen_member"],
        )
        .group_by("antigen_rep")
        .agg(pl.col("antigen_member"))
        .sort(by=pl.col("antigen_member").list.len(), descending=True)
    )
    heavy_df = (
        pl.read_csv(
            args.heavy_rep_tsv,
            separator="\t",
            has_header=False,
            new_columns=["heavy_rep", "heavy_member"],
        )
        .group_by("heavy_rep")
        .agg(pl.col("heavy_member"))
        .sort(by=pl.col("heavy_member").list.len(), descending=True)
    )
    light_df = (
        pl.read_csv(
            args.light_rep_tsv,
            separator="\t",
            has_header=False,
            new_columns=["light_rep", "light_member"],
        )
        .group_by("light_rep")
        .agg(pl.col("light_member"))
        .sort(by=pl.col("light_member").list.len(), descending=True)
    )

    keep_ids = []

    for row in antigen_df.iter_rows(named=True):

        id = row["antigen_rep"]

        focal_heavy = heavy_df.filter(pl.col("heavy_member").list.contains(id))

        if focal_heavy.height == 0:
            continue

        heavy_df = heavy_df.filter(
            pl.col("heavy_rep") != focal_heavy.select("heavy_rep").item()
        )

        keep_ids.append(id)

    tmp_id = pl.DataFrame({"complexed_antibody_id": keep_ids})

    chosen = orig_pq.join(tmp_id, on="complexed_antibody_id")

    discard = orig_pq.join(tmp_id, on="complexed_antibody_id", how="anti")

    chosen.write_parquet("chosen_complex.parquet")
    discard.write_parquet("discard_complex.parquet")
