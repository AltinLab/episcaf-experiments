#!/usr/bin/env python
from mdaf3.FeatureExtraction import serial_apply
import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP
import polars as pl
import argparse

from scipy.spatial.distance import squareform
from pathlib import Path

from MDAnalysis.analysis.distances import self_distance_array
import math

"""
Intuition:

min num residues that connect two epitope residues is (distance between res Calpha) / (avg len AA) 

NOTE: modify me to produce 100 contigs per seq
"""
# https://www.researchgate.net/publication/6767198_Contour_Length_and_Refolding_Rate_of_a_Small_Protein_Controlled_by_Engineered_Disulfide_Bonds?_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6InF1ZXN0aW9uIiwicGFnZSI6InF1ZXN0aW9uIn19
ANGSTROM_PER_LINEAR_AA = 3.6
ANGSTROM_PER_HELICAL_AA = 1.5


def fill_gaps(contig_str, n_gaps, max_res, consumed_res, min_bridgeres):
    """
    Takes a string containing tokens of the form:

    GAP<integer>

    the string contains n_gaps such tokens, where the integer of each gap token monotonically increases from 1 to n_gaps (inclusive)

    Replaces gaps with some random integer number of residues s.t. the total number of residues added to the string = max_res - consumed_res
    """
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


def gen_contig_string(row, top_path, max_res):
    u = mda.Universe(top_path / (row["complexed_antibody_id"] + ".pdb"))

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
            if start_i is None:
                start_i = ri1
        else:
            if start_i is not None:
                contig_str += (
                    f"{gap_token}/A{u.residues[start_i].resid}-{u.residues[ri1].resid}/"
                )
                consumed_res += ri1 - start_i + 1
            else:
                contig_str += (
                    f"{gap_token}/A{u.residues[ri1].resid}-{u.residues[ri1].resid}/"
                )
                consumed_res += 1
            min_bridgeres.append(
                math.ceil(dist[ri1 - aro, ri2 - aro] / ANGSTROM_PER_HELICAL_AA)
            )
            n_chunks += 1
            gap_token = f"GAP{n_chunks + 1}"
            start_i = None

    # handle last residue
    last_ri = ep_res_idx[-1]
    gap_tok2 = f"GAP{n_chunks + 2}"
    if start_i is not None:
        contig_str += f"{gap_token}/A{u.residues[start_i].resid}-{u.residues[last_ri].resid}/{gap_tok2}/0 "
        consumed_res += last_ri - start_i + 1
    else:
        contig_str += f"{gap_token}/A{u.residues[last_ri].resid}-{u.residues[last_ri].resid}/{gap_tok2}/0 "
        consumed_res += 1
    n_chunks += 1

    contig_str += f"B1-{len(u.select_atoms('segid B').residues)}/0 "
    contig_str += f"C1-{len(u.select_atoms('segid C').residues)}"

    contig_str = fill_gaps(
        contig_str, n_chunks + 1, max_res, consumed_res, min_bridgeres
    )

    if contig_str is None:
        row["contig_string"] = None
        row["epitope_chunks"] = None
    else:
        row["contig_string"] = contig_str
        row["epitope_chunks"] = n_chunks

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

    complexes = serial_apply(complexes, gen_contig_string, top_path, 104)

    # some drop out due to implausibility
    complexes = complexes.filter(pl.col("contig_string").is_not_null())

    complexes.write_parquet(args.output_pq)
