#!/usr/bin/env python
"""
Given a dataframe with at least the following columns:
- id: unique identifier of topology/row
- epitope_resindices: 0-indexed indices for residues in the PDB file which are epitope residues

And a PDB directory containing a file for each id

Generate args.contig_per_pdb randomly generated contig strings for each topology file/row

The strategy for random generation is as follows:

- Aggregate epitope resindices into a series of epitope "chunks"
  - if two subsequent epitope resindices are 1. within 10 residues or 2. part of the same beta sheet or alpha helix
    as one another, group them into the same chunk
  - else, consider there to be a gap between them which needs to have a scaffold diffused
- Each gap is allocated some number of residues to diffuse using this process
    - given the total number of target residues, subtract the number of residues already allocated to "epitope chunks"
    - assign a minimum number of residues for each gap: 2 * (GAP_DIST / ANGSTROM_PER_LINEAR_AA) (min_res)
      where GAP_DIST is the distance between the last calpha of the opening chunk and the first calpha of the closing chunk (first/last by sequence)
      this quantity is the minimum found from looking at RFdiffusion bencmarks (see notebooks/paper_analysis/gap_len.ipynb)
    - for all gaps (excluding the opening/closing "gaps" of the scaffold, which are handled later):
      - randomly select a gap from a list of gaps
      - assign it a random number of residues somewhere between min_res and min((15 * (min_res // 2)), remaining_res), where
        remaining_res is the number of residues not already "spent" on epitope chunks or other gaps
    - after all non-opening and closing gaps have been filled, evenly distribute the remaining_res to the opening and closing gaps

This results in a dataframe with the columns of args.input_pq +
- epitope_chunk_resindices, List(int): 0-indexed indices for residues in the PDB file which are considered part of epitope chunks
- contig_string: the RFdiffusion-formatted contig string
- epitope_chunks: the number of distinct chunks of epitopes identified
- contig_gaps, List(int): the number of residues allocated to each gap, including the opening/closing gaps
- contig_min, List(int): the minimum number of residues that could have been allocated to each gap by definition above
- contig_id, List(int): list of unique integers identifying each contig (unique per id)
- contig_length, str: RFdiffusion-formatted length string (i.e. 104-104)
- scaffolded_epitope_chunk_resindices: 0-indexed indices for residues in the scaffolded_epitope_seq which are part of the epitope chunks
- scaffolded_epitope_resindices: 0-indexed indices for residues in the scaffolded_epitope_seq which are part of the epitope

Each id in args.input_pq will have args.contig_per_pdb rows in the output parquet file. If after 100000 attempts, no new contig strings
are found but the number found is less than args.contig_per_pdb, the number already found is accepted.
"""
from mdaf3.FeatureExtraction import serial_apply
import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP
import polars as pl
import argparse

from scipy.spatial.distance import squareform
from pathlib import Path

from MDAnalysis.analysis.distances import self_distance_array
import math
import random
import numpy as np

import warnings

random.seed(42)

# https://www.researchgate.net/publication/6767198_Contour_Length_and_Refolding_Rate_of_a_Small_Protein_Controlled_by_Engineered_Disulfide_Bonds?_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6InF1ZXN0aW9uIiwicGFnZSI6InF1ZXN0aW9uIn19
ANGSTROM_PER_LINEAR_AA = 3.6
ANGSTROM_PER_HELICAL_AA = 1.5


def fill_gaps(contig_str, n_gaps, max_res, consumed_res, min_bridgeres):
    remaining_res = max_res - consumed_res
    if remaining_res - sum(min_bridgeres) <= 0:
        return None

    # replace all except opening and closing gaps
    # with min_bridgres
    for i in range(2, n_gaps):
        # must match the slash
        contig_str = contig_str.replace(
            f"GAP{i}/", f"{min_bridgeres[i-2]}-{min_bridgeres[i-2]}/", 1
        )

    # evenly split remaining residues between opening and closing gap
    # if odd, give the extra residue to opener
    remaining_res -= sum(min_bridgeres)
    if remaining_res % 2 == 0:
        opening_res = remaining_res // 2
        closing_res = remaining_res // 2
    else:
        opening_res = (remaining_res // 2) + 1
        closing_res = remaining_res // 2

    if closing_res > 0:
        contig_str = contig_str.replace("GAP1/", f"{opening_res}-{opening_res}/", 1)
        contig_str = contig_str.replace(
            f"GAP{n_gaps}/", f"{closing_res}-{closing_res}/", 1
        )
    else:
        contig_str.replace("GAP1/", "", 1)
        contig_str = contig_str.replace(f"GAP{n_gaps}/0", "0", 1)
    return contig_str


def fill_gaps_random(
    contig_str, n_gaps, max_res, consumed_res, min_bridgeres, no_antibody
):
    """
    Takes a string containing tokens of the form:

    GAP<integer>

    the string contains n_gaps such tokens, where the integer of each gap token monotonically increases from 1 to n_gaps (inclusive)

    Replaces gaps with some random integer number of residues s.t. the total number of residues added to the string = max_res - consumed_res
    """
    remaining_res = max_res - consumed_res
    if remaining_res - sum(min_bridgeres) <= 0:
        return None

    gap_ints = list(range(2, n_gaps))
    remaining_min_bridgeres = sum(min_bridgeres)
    rand_bridgeres_ls = [0 for _ in range(len(min_bridgeres))]

    while gap_ints:
        i = random.choice(gap_ints)
        gap_ints.remove(i)

        min_res = min_bridgeres[i - 2]
        remaining_min_bridgeres -= min_res
        max_res = min(remaining_res - remaining_min_bridgeres, 15 * (min_res // 2))

        if max_res - min_res < 0:
            raise ValueError("max_res - min_res <= 0")

        rand_bridgeres = random.randint(min_res, max_res)
        rand_bridgeres_ls[i - 2] = rand_bridgeres

        remaining_res -= rand_bridgeres
        contig_str = contig_str.replace(
            f"GAP{i}/", f"{rand_bridgeres}-{rand_bridgeres}/", 1
        )

    # evenly split remaining residues between opening and closing gap
    # if odd, give the extra residue to opener
    if remaining_res < 0:
        raise ValueError("Can't have less than 0 remaining res")
    if remaining_res % 2 == 0:
        opening_res = remaining_res // 2
        closing_res = remaining_res // 2
    else:
        opening_res = (remaining_res // 2) + 1
        closing_res = remaining_res // 2

    if closing_res > 0:
        if not no_antibody:
            contig_str = contig_str.replace(
                f"GAP{n_gaps}/", f"{closing_res}-{closing_res}/", 1
            )
        else:
            contig_str = contig_str.replace(
                f"GAP{n_gaps}/0 ", f"{closing_res}-{closing_res}", 1
            )
    else:
        if not no_antibody:
            contig_str = contig_str.replace(f"GAP{n_gaps}/0", "0", 1)
        else:
            contig_str = contig_str.replace(f"/GAP{n_gaps}/0 ", "", 1)
    if opening_res > 0:
        contig_str = contig_str.replace("GAP1/", f"{opening_res}-{opening_res}/", 1)

    else:
        contig_str = contig_str.replace("GAP1/", "", 1)

    rand_bridgeres_ls = [opening_res] + rand_bridgeres_ls + [closing_res]

    return contig_str, rand_bridgeres_ls


def apply_scaffold_index_offsets(ep_chunk_ls, bridge_ls, is_original_ls):

    cumulative_offset = bridge_ls[0]
    # offsets.append(cumulative_offset)
    ep_chunk = []
    ep = []

    i = 1
    for is_original, chunk in zip(is_original_ls, ep_chunk_ls):
        chunk_adjusted = [res - chunk[0] + cumulative_offset for res in chunk]
        ep_chunk.extend(chunk_adjusted)
        ep.extend(np.array(chunk_adjusted)[is_original].tolist())
        cumulative_offset += bridge_ls[i]
        cumulative_offset += len(chunk)
        i += 1

    return ep_chunk, ep


def gen_contig_string(row, top_path, target_length, contig_per_pdb, no_antibody):
    u = mda.Universe(top_path / (row["id"] + ".pdb"))

    antigen_ag = u.select_atoms("segid A")

    sec_struct_results = DSSP(antigen_ag).run()

    sec_struct = sec_struct_results.results.dssp_ndarray[0]

    ep_res_idx = row["epitope_resindices"]

    dist = squareform(self_distance_array(antigen_ag.select_atoms("name CA")))
    # ANTIGEN RESINDEX OFFSET
    # SHOULD BE 0 BUT JUST TO BES SAFE WHEN ACCESSING DISTANCE ARRAY
    # DO dist[resindex - aro]
    aro = u.select_atoms("segid A").residues[0].resindex

    # keep_resindices
    # # use a few heuristics to get good residue painting
    # # 1. if two residues are a part of the same secondary structure element,
    # # no matter their distance in AAs, keep the residues between them
    # # 2. if two residues are within 10 residues, keep the residues between them
    # # 3. otherwise, don't keep the residues between them
    start_i = None
    contig_str = ""
    consumed_res = 0
    min_bridgeres = []
    n_chunks = 0
    gap_token = f"GAP{n_chunks + 1}"
    epitope_chunk_resindices = []
    residx_cache = []
    is_original_ls = []

    for i in range(len(ep_res_idx) - 1):
        ri1 = ep_res_idx[i]
        ri2 = ep_res_idx[i + 1]

        # both are part of same sheet or helix
        # or are close in sequence
        if (
            sec_struct[ri1 : ri2 + 1][:, 1].all()
            or sec_struct[ri1 : ri2 + 1][:, 2].all()
            or ri2 - ri1 <= 10
        ):
            residx_cache.append(ri1)
            if start_i is None:
                start_i = ri1

        else:
            if start_i is not None:
                contig_str += (
                    f"{gap_token}/A{u.residues[start_i].resid}-{u.residues[ri1].resid}/"
                )
                ep_range = list(range(start_i, ri1 + 1))
                epitope_chunk_resindices.append(ep_range)
                is_original_ls.append(np.isin(np.array(ep_range), residx_cache + [ri1]))

                consumed_res += ri1 - start_i + 1
            else:
                contig_str += (
                    f"{gap_token}/A{u.residues[ri1].resid}-{u.residues[ri1].resid}/"
                )
                epitope_chunk_resindices.append([ri1])
                is_original_ls.append(np.array([True]))
                consumed_res += 1
            min_bridgeres.append(
                math.ceil(dist[ri1 - aro, ri2 - aro] / ANGSTROM_PER_LINEAR_AA) * 2
            )
            n_chunks += 1
            gap_token = f"GAP{n_chunks + 1}"
            start_i = None
            residx_cache = []

    # handle last residue
    last_ri = ep_res_idx[-1]
    gap_tok2 = f"GAP{n_chunks + 2}"
    if start_i is not None:
        contig_str += f"{gap_token}/A{u.residues[start_i].resid}-{u.residues[last_ri].resid}/{gap_tok2}/0 "
        ep_range = list(range(start_i, last_ri + 1))
        epitope_chunk_resindices.append(ep_range)
        is_original_ls.append(np.isin(np.array(ep_range), residx_cache + [last_ri]))
        consumed_res += last_ri - start_i + 1
    else:
        contig_str += f"{gap_token}/A{u.residues[last_ri].resid}-{u.residues[last_ri].resid}/{gap_tok2}/0 "
        epitope_chunk_resindices.append([last_ri])
        is_original_ls.append(np.array([True]))
        consumed_res += 1
    n_chunks += 1

    if not no_antibody:
        for i, segid in enumerate(u.segments.segids[1:]):
            contig_str += f"{segid}1-{len(u.select_atoms(f'segid {segid}').residues)}"
            if i != len((u.segments)) - 2:
                contig_str += "/0 "

    contig_str_varaints = []
    contig_gaps = []
    scaffolded_epitope_chunk_resindices_ls = []
    scaffolded_epitope_resindices_ls = []

    attempts = 0
    while len(contig_str_varaints) < contig_per_pdb and attempts < 100000:

        result = fill_gaps_random(
            contig_str,
            n_chunks + 1,
            target_length,
            consumed_res,
            min_bridgeres,
            no_antibody,
        )
        if result is None:
            contig_str_varaints = None
            break

        candidate_contig_str, candidate_contig_gaps = result

        if candidate_contig_str in contig_str_varaints:
            attempts += 1
            continue
        else:
            contig_str_varaints.append(candidate_contig_str)
            contig_gaps.append(candidate_contig_gaps)

            scaffolded_epitope_chunk_resindices, scaffolded_epitope_resindices = (
                apply_scaffold_index_offsets(
                    epitope_chunk_resindices, candidate_contig_gaps, is_original_ls
                )
            )
            scaffolded_epitope_chunk_resindices_ls.append(
                scaffolded_epitope_chunk_resindices
            )
            scaffolded_epitope_resindices_ls.append(scaffolded_epitope_resindices)

    if contig_str_varaints is None:
        warnings.warn(
            f"Target length {target_length} too short for {row['id']}, skipping"
        )
        row["contig_string"] = None
        row["epitope_chunks"] = None
        row["contig_gaps"] = None
        row["contig_min"] = None
        row["contig_id"] = None
        row["contig_length"] = None
        row["epitope_chunk_resindices"] = None
        row["scaffolded_epitope_chunk_resindices"] = None
        row["scaffolded_epitope_resindices"] = None
    else:
        row["contig_string"] = contig_str_varaints
        row["epitope_chunks"] = n_chunks
        row["contig_gaps"] = contig_gaps
        row["contig_min"] = [0] + min_bridgeres + [0]
        row["contig_id"] = list(range(len(contig_str_varaints)))
        row["contig_length"] = f"{target_length}-{target_length}"
        row["epitope_chunk_resindices"] = [
            x for sublist in epitope_chunk_resindices for x in sublist
        ]
        row["scaffolded_epitope_chunk_resindices"] = (
            scaffolded_epitope_chunk_resindices_ls
        )
        row["scaffolded_epitope_resindices"] = scaffolded_epitope_resindices_ls

    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_pq",
        type=str,
    )
    parser.add_argument("--topology_path", type=str)
    parser.add_argument("--output_pq", type=str)
    parser.add_argument("--target_length", type=int)
    parser.add_argument("--contig_per_pdb", type=int)
    parser.add_argument("--no_antibody", action="store_true")
    args = parser.parse_args()

    complexes = pl.read_parquet(args.input_pq).sort(by="id")
    top_path = Path(args.topology_path)
    no_antibody = args.no_antibody

    complexes = serial_apply(
        complexes,
        gen_contig_string,
        top_path,
        args.target_length,
        args.contig_per_pdb,
        no_antibody,
    )

    # some drop out due to implausibility
    valid_complexes = complexes.filter(pl.col("contig_string").is_not_null())

    valid_complexes = valid_complexes.explode(
        [
            "contig_string",
            "contig_gaps",
            "contig_id",
            "scaffolded_epitope_chunk_resindices",
            "scaffolded_epitope_resindices",
        ]
    )

    valid_complexes.write_parquet(args.output_pq)
