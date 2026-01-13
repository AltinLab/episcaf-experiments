#!/usr/bin/env python
"""
Perform 3 steps:
- remove HETATM and collapse occupancies
- Reorder the topologies so that Antibody = chain A, Heavy chain = Chain B, light chain = chain C
- Collapse down the heavy and light chains to just their CDR loops
"""
import anarci
import numpy as np
import MDAnalysis as mda
import polars as pl
import argparse
from pathlib import Path
from collections import OrderedDict


def annotate_human_antibody(seq, resindices):
    results = anarci.run_anarci(
        [("AAA", seq)],
        scheme="martin",
        allowed_species=["human"],
    )
    if results[1][0] is None:
        raise ValueError(f"No domain found for sequence '{seq}'")
    elif len(results[1][0]) > 1:
        raise ValueError("Multiple domains found for sequence")

    numbering = results[1][0][0]

    sub_start = numbering[1]
    sub_stop = numbering[2] + 1

    # this entire substring will be numbered
    # however, there may be gaps in the sequence
    # which are given an IMGT number
    numbered_substring = seq[sub_start:sub_stop]
    resindices_slice = resindices[sub_start:sub_stop]
    imgt_num = np.zeros((len(range(sub_start, sub_stop)),), dtype=np.int32)
    imgt_tuples = numbering[0]
    j = 0
    for i in range(len(numbered_substring)):
        aa = numbered_substring[i]
        while j < len(imgt_tuples) and imgt_tuples[j][1] != aa:
            j += 1
        if j < len(imgt_tuples):
            imgt_num[i] = imgt_tuples[j][0][0]
            j += 1

    # zero is not an IMGT number, so we use this as a quick error cehck
    n_zeroes = np.count_nonzero(imgt_num == 0)
    if n_zeroes != 0:
        raise ValueError(
            f"0 is not an IMGT number. Numbering failed for antibody chain '{seq}'"
        )

    return resindices_slice, imgt_num, (sub_start, sub_stop)


def reorder_universe(universe, reorder_map):

    chain_us = []

    reorder_map_srt = OrderedDict(sorted(reorder_map.items(), key=lambda x: x[1]))

    for from_segid, to_segid in reorder_map_srt.items():

        if from_segid not in universe.segments.segids:
            raise ValueError

        chain_u = mda.Merge(universe.select_atoms(f"segid {from_segid}"))
        chain_u.segments.segids = to_segid
        chain_u.atoms.chainIDs = [to_segid] * len(chain_u.atoms)
        chain_us.append(chain_u.atoms)

    return mda.Merge(*chain_us)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--complex_pq",
        type=str,
    )
    parser.add_argument("--topology_path", type=str)
    parser.add_argument("--output_dir", type=str)
    args = parser.parse_args()

    complexes = pl.read_parquet(args.complex_pq)
    top_path = Path(args.topology_path)
    outdir = Path(args.output_dir)

    for row in complexes.iter_rows(named=True):

        orig_u = mda.Universe(top_path / (row["complexed_antibody_id"] + ".pdb"))

        reorder_map = {
            row["antigen_segid"]: "A",
            row["heavy_segid"]: "B",
            row["light_segid"]: "C",
        }
        fmt_u = reorder_universe(orig_u, reorder_map)

        # only atoms (not HETATM) and choose single occupancy
        fmt_u = fmt_u.select_atoms(
            "record_type ATOM and (altloc A or not altloc [!?])"
        ).universe

        with mda.Writer(outdir / (row["complexed_antibody_id"] + ".cleaned.pdb")) as W:
            W.write(fmt_u.atoms)
