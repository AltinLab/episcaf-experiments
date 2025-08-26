#!/usr/bin/env python
"""
For a site to be an attractive candidate for binding it should have >~3 hydrophobic residues for the binder to interact with.
Binding to charged polar sites is still quite hard.
Binding to sites with glycans close to them is also hard since they often become ordered upon binding and you will take an energetic hit for that.
Historically, binder design has also avoided unstructured loops, it is not clear if this is still a
requirement as RFdiffusion has been used to bind unstructured peptides which share a lot in common with unstructured loops.

RFdiffusion scales in runtime as O(N^2) where N is the number of residues in your system.
As such, it is a very good idea to truncate large targets so that your computations are not unnecessarily expensive.

Hotspots are a feature that we integrated into the model to allow for the control of the site on the target which the binder will interact with.
In the paper we define a hotspot as a residue on the target protein which is within 10A Cbeta distance of the binder.

Of all of the hotspots which are identified on the target 0-20% of these hotspots are actually provided to the model and the rest are masked.
This is important for understanding how you should pick hotspots at inference time.; the model is expecting to have to make more contacts than you specify.
We normally recommend between 3-6 hotspots, you should run a few pilot runs before generating thousands of designs to make sure the number of hotspots you are providing will give results you like.

For all in silico benchmarks in this paper, we use the AF2 structure prediction network17
for validation and define an in silico 'success' as an RFdiffusion output for which the AF2 structure
predicted from a single sequence is (1) of high confidence (mean predicted aligned error (pAE), less than five),
(2) globally within a 2 Å backbone root mean-squared deviation (r.m.s.d.) of the designed structure and (3)
within 1 Å backbone r.m.s.d. on any scaffolded functional site (Supplementary Methods).
This measure of in silico success has been found to correlate with experimental success4,7,26 and is
significantly more stringent than template modelling (TM)-score-based metrics used elsewhere5,16,27,28,29
(Supplementary Fig. 2c,d).
"""

from mdaf3.FeatureExtraction import serial_apply, split_apply_combine
from MDAnalysis.lib.util import convert_aa_code
import MDAnalysis as mda
import requests
from io import StringIO
import polars as pl
import re
from pathlib import Path
import argparse
import shutil
from collections import OrderedDict
import numpy as np


def parse_abdb_cluster_txt(path):
    clust_df = pl.read_csv(
        path, comment_prefix="#", separator="\t", new_columns=["blob"]
    )

    clust_df = clust_df.with_columns(
        pl.col("blob")
        .str.split_exact(":", 1)
        .struct.rename_fields(["free_antibody_id_blob", "complexed_antibody_id_blob"])
        .alias("fields")
    ).unnest("fields")

    clust_df = clust_df.with_columns(
        pl.when(pl.col("complexed_antibody_id_blob") == "")
        .then(pl.lit([]))
        .otherwise(pl.col("complexed_antibody_id_blob").str.split(","))
        .alias("complexed_antibody_id_list"),
        pl.when(pl.col("free_antibody_id_blob") == "")
        .then(pl.lit([]))
        .otherwise(pl.col("free_antibody_id_blob").str.split(","))
        .alias("free_antibody_id_list"),
    )

    clust_df = clust_df.select(["free_antibody_id_list", "complexed_antibody_id_list"])
    return clust_df


def pdb_has_LH(row, abdb_dir):
    with open(abdb_dir / ("pdb" + row["complexed_antibody_id"] + ".faa")) as f:
        for line in f:
            if line.startswith(">"):
                if "L_H" in line:
                    row["has_LH"] = True
                    return row

    row["has_LH"] = False
    return row


def filter_for_LH(df, abdb_dir):

    df = serial_apply(df, pdb_has_LH, abdb_dir)

    return df.filter(pl.col("has_LH")).select(pl.exclude("has_LH"))


def parse_pdb_metadata(row, abdb_dir):

    # SEQRES used to get full sequence http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_35.html
    file_rows = []
    pat = re.compile(r"^REMARK\s+950\s+CHAIN\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")

    errdict = {
        "resolution": None,
        "light_segid": None,
        "light_original_segid": None,
        "heavy_segid": None,
        "heavy_original_segid": None,
        "antigen_segid": None,
        "antigen_original_segid": None,
        "light_seq": None,
        "heavy_seq": None,
        "antigen_seq": None,
    }

    try:
        with open(abdb_dir / ("pdb" + row["complexed_antibody_id"] + ".mar"), "r") as f:
            for line in f:
                if not line.startswith("REMARK"):
                    if line.startswith("MODRES"):
                        # just removed antibodies with modified residues
                        row.update(errdict)
                        return row
                    else:
                        break

                file_rows.append(line)

        res = re.compile(r".*(\d+(?:\.\d+))A").match(file_rows[0]).group(1)
        row["resolution"] = float(res)

        pat = re.compile(r"^REMARK\s+950\s+CHAIN\s+(\S+)\s+(\S+)\s+(\S+)\s*$")

        light_segid = None
        heavy_segid = None
        antigen_segid = None

        for dat in file_rows[2:]:
            m = pat.match(dat)
            chain_type, segid, original_segid = m.groups()

            if chain_type == "L":
                if light_segid is not None:
                    row.update(errdict)
                    return row  # multiple light chains, not supported
                light_segid = segid
                light_original_segid = original_segid
            elif chain_type == "H":
                if heavy_segid is not None:
                    row.update(errdict)
                    return row
                heavy_segid = segid
                heavy_original_segid = original_segid
            elif chain_type == "A":
                if antigen_segid is not None:
                    row.update(errdict)
                    return row
                antigen_segid = segid
                antigen_original_segid = original_segid

        row["light_segid"] = light_segid
        row["light_original_segid"] = light_original_segid
        row["heavy_segid"] = heavy_segid
        row["heavy_original_segid"] = heavy_original_segid
        row["antigen_segid"] = antigen_segid
        row["antigen_original_segid"] = antigen_original_segid

        res_rows = []
        with open(abdb_dir / ("pdb" + row["complexed_antibody_id"] + ".mar"), "r") as f:
            for line in f:
                if line.startswith("SEQRES"):
                    res_rows.append(line)

        pat = re.compile(r"^SEQRES\s+\d+\s+(\S+)\s+(\d+)\s+(.+)")

        seq_dict = {
            row["light_segid"]: "",
            row["heavy_segid"]: "",
            row["antigen_segid"]: "",
        }

        for seq_str in res_rows:
            m = pat.match(seq_str)
            chain_type, res_count, seq = m.groups()

            seq = seq.strip().split(" ")

            try:
                seq = "".join([convert_aa_code(aa) for aa in seq])
            except ValueError:
                row.update(errdict)
                return row

            try:
                seq_dict[chain_type] += seq
            except KeyError:
                print("should not happen")

        row["light_seq"] = seq_dict[row["light_segid"]]
        row["heavy_seq"] = seq_dict[row["heavy_segid"]]
        row["antigen_seq"] = seq_dict[row["antigen_segid"]]

    # some PDB IDs only have a ".bad" extension
    # I'm assuming these were incorrectly parsed by PDBD intake for some reason
    # and annotated as such- filtering this out
    except FileNotFoundError:
        row.update(errdict)
        return row

    return row


POLYMER_QUERY = """
    query($id: String!) {
      entry(entry_id: $id) {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            entity_id
            asym_ids
            auth_asym_ids
          }
          entity_poly {
            pdbx_seq_one_letter_code_can
          }
          rcsb_entity_source_organism {
            ncbi_taxonomy_id
          }
        }
      }
    }
    """

R_FACTOR_QUERY = """
    query($id: String!) {
      entry(entry_id: $id) {
        refine {
            ls_R_factor_R_work
            ls_R_factor_R_free
        }
      }
    }
"""


def query_pdb_metadata(row):
    pdb_id = row["complexed_antibody_id"][:4].upper()
    r = requests.post(
        "https://data.rcsb.org/graphql",
        json={"query": R_FACTOR_QUERY, "variables": {"id": pdb_id}},
        timeout=120,
    )
    r.raise_for_status()

    refine = r.json()["data"]["entry"]["refine"]

    if refine is None:
        row["r_work"] = None
        row["r_free"] = None

    else:
        row["r_work"] = refine[0]["ls_R_factor_R_work"]
        row["r_free"] = refine[0]["ls_R_factor_R_free"]

    r = requests.post(
        "https://data.rcsb.org/graphql",
        json={"query": POLYMER_QUERY, "variables": {"id": pdb_id}},
        timeout=120,
    )
    r.raise_for_status()
    data = r.json()["data"]["entry"]["polymer_entities"]

    dat_dict = {}

    for polymer in data:
        for auth_asym_id in polymer["rcsb_polymer_entity_container_identifiers"][
            "auth_asym_ids"
        ]:
            dat_dict[auth_asym_id] = polymer

    for chain in ["antigen", "light", "heavy"]:
        asym_id = row[f"{chain}_original_segid"]

        if asym_id not in dat_dict:
            # some cases ABDB didn't parse correctly
            # i.e. 7T3M
            # row[f"{chain}_seq"] = None
            row[f"{chain}_ncbi_taxonomy_id"] = None

        else:

            polymer_data = dat_dict[asym_id]
            # seq data better parsed from SEQRES
            # row[f"{chain}_seq"] = polymer_data["entity_poly"][
            #     "pdbx_seq_one_letter_code_can"
            # ]

            spec = (
                polymer_data["rcsb_entity_source_organism"][0]["ncbi_taxonomy_id"]
                if polymer_data["rcsb_entity_source_organism"] is not None
                else None
            )

            if spec is None:
                row[f"{chain}_ncbi_taxonomy_id"] = None
            else:
                row[f"{chain}_ncbi_taxonomy_id"] = int(spec)

    return row


# def write_pdbfile(row, abdb_dir, output_dir):
#     original = Path(abdb_dir / ("pdb" + row["complexed_antibody_id"] + ".mar"))
#     new = Path(output_dir / (row["complexed_antibody_id"] + ".pdb"))
#     shutil.copy2(original, new)
#     return row


def reorder_universe(universe, reorder_map):
    """
    Change the segment IDs in universe using reorder map
    which gives from_segid: to_segid pairs

    Chains will be added to the new universe in alphabetical
    order of to_segids

    Custom resids are removed and replaced with monotonically
    increasing integers- this strips i.e. antibody numbering
    (which MDAnalysis can't properly parse anyways since it requires
    letters in resids)
    """

    chain_us = []

    reorder_map_srt = OrderedDict(sorted(reorder_map.items(), key=lambda x: x[1]))

    offset = 0
    for from_segid, to_segid in reorder_map_srt.items():

        if from_segid not in universe.segments.segids:
            raise ValueError

        chain_u = mda.Merge(universe.select_atoms(f"segid {from_segid}"))
        # segindex, resindex automatically update
        chain_u.segments.segids = to_segid
        chain_u.residues.resids = np.arange(1, len(chain_u.residues) + 1)
        chain_u.residues.resnums = np.arange(1, len(chain_u.residues) + 1)

        # atom_resnums = []
        # for i, res in enumerate(chain_u.residues):
        #     atom_resnums.extend([offset + i + 1] * len(res.atoms))

        chain_u.atoms.chainIDs = [to_segid] * len(chain_u.atoms)
        # chain_u.atoms.resnums = np.array(atom_resnums)
        # chain_u.atoms.resids = np.array(atom_resnums)
        chain_us.append(chain_u.atoms)

        offset += len(chain_u.residues)

    return mda.Merge(*chain_us)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--abdb_snapshot_dir",
        type=str,
    )
    parser.add_argument("--output_parquet", type=str)
    parser.add_argument("--output_pdb_dir", type=str)

    args = parser.parse_args()

    abdb_dir = Path(args.abdb_snapshot_dir)

    clust_df = parse_abdb_cluster_txt(abdb_dir / "AbClusters.txt")

    complex_df = (
        clust_df.filter(pl.col("complexed_antibody_id_list").list.len() > 0)
        .select("complexed_antibody_id_list")
        .explode("complexed_antibody_id_list")
        .rename({"complexed_antibody_id_list": "complexed_antibody_id"})
    )

    # keep only complexes that are bound to ONLY a protein
    complex_df = complex_df.filter(
        pl.col("complexed_antibody_id").str.contains(r"_\d+P$")
    )

    # keep only L+H complexes
    complex_df = filter_for_LH(complex_df, abdb_dir)

    # extract metadata from ABDB annotations
    # label the L,H, and (A)ntigen chain
    complex_df = split_apply_combine(
        complex_df, parse_pdb_metadata, abdb_dir, chunksize=15
    )

    # filters out BAD files (see parse_pdb_metadata), modified res antibodies, duplicate chain entries
    complex_df = complex_df.filter(pl.col("resolution").is_not_null())

    # filter for resolution lower than 3 angstroms
    complex_df = complex_df.filter(pl.col("resolution") <= 3.0)

    # query PDB for additional metadata
    # (serial for network IO to avoid failure)
    complex_df = serial_apply(complex_df, query_pdb_metadata)

    # filter for only standard residues
    complex_df = complex_df.filter(
        ~pl.col("antigen_seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]"),
        ~pl.col("heavy_seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]"),
        ~pl.col("light_seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]"),
    )

    # filter for antigen seq len >= 39 (from bepipred3 paper)
    complex_df = complex_df.filter(pl.col("antigen_seq").str.len_chars() > 104)

    # keep only human complexes
    # filter for r-factor <= 0.3
    complex_df = complex_df.filter(
        pl.col("light_ncbi_taxonomy_id") == 9606,
        pl.col("heavy_ncbi_taxonomy_id") == 9606,
        pl.col("r_work") <= 0.3,
    )

    missing_backbone = []

    # minimally clean the PDBfiles and write them (no cropping yet)
    for row in complex_df.iter_rows(named=True):

        orig_pdb_path = Path(abdb_dir / ("pdb" + row["complexed_antibody_id"] + ".mar"))
        orig_u = mda.Universe(orig_pdb_path, topology_format="pdb")

        reorder_map = {
            row["antigen_segid"]: "A",
            row["heavy_segid"]: "B",
            row["light_segid"]: "C",
        }
        fmt_u = reorder_universe(orig_u, reorder_map)

        # only atoms (not HETATM) and choose single occupancy
        fmt_ag = fmt_u.select_atoms(
            "record_type ATOM and (altloc A or not altloc [!?])"
        )

        for res in fmt_ag.residues:
            an = list(res.atoms.names)
            if (
                an.count("N") != 1
                or an.count("C") != 1
                or an.count("O") != 1
                or an.count("CA") != 1
            ):
                missing_backbone.append(row["complexed_antibody_id"])
                continue

        with mda.Writer(
            Path(args.output_pdb_dir) / (row["complexed_antibody_id"] + ".pdb")
        ) as W:
            W.write(fmt_ag)

    tmp_df = pl.DataFrame({"complexed_antibody_id": missing_backbone})

    complex_df = complex_df.join(tmp_df, on="complexed_antibody_id", how="anti")

    complex_df.select(
        pl.exclude(["light_segid", "heavy_segid", "antigen_segid"])
    ).write_parquet(args.output_parquet)
