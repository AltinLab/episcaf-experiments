from mdaf3.AF3OutputParser import AF3Output
from mdaf3.FeatureExtraction import serial_apply, split_apply_combine
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.align import alignto
import numpy as np
import argparse
import polars as pl
from pathlib import Path

# compute RMSD of all MPNN res vs all AF3 res


def compute_epitope_chunk_rmsd_mpnn(row, mpnn_dir, af3_dir):
    try:
        af3_output = AF3Output(af3_dir / row["assay_scaffolded_epitope_id"])
        af3_u = af3_output.get_mda_universe()
    except Exception:
        row["overall_rmsd"] = None
        row["epitope_chunk_rmsd_vs_mpnn"] = None
        return row

    mpnn_u = mda.Universe(
        ((mpnn_dir / row["id"]) / str(row["contig_id"]))
        / (
            row["id"]
            + "_"
            + str(row["rfd_id"])
            + "_fixed_dldesign_"
            + str(row["mpnn_id"])
            + ".pdb"
        )
    )
    mpnn_ris = row["scaffolded_epitope_chunk_resindices"]
    af3_ris = row["assay_scaffolded_epitope_chunk_resindices"]

    assert mpnn_u.residues[mpnn_ris].sequence(format="string") == af3_u.residues[
        af3_ris
    ].sequence(format="string")

    start_idx = 2
    end_idx = start_idx + len(mpnn_u.select_atoms("segid A").residues)

    assert mpnn_u.residues.sequence(format="string") == af3_u.residues[
        start_idx:end_idx
    ].sequence(format="string")

    row["overall_rmsd"] = rms.rmsd(
        mpnn_u.select_atoms("segid A").select_atoms("backbone").positions,
        af3_u.residues[start_idx:end_idx].atoms.select_atoms("backbone").positions,
        superposition=True,
    )
    row["epitope_chunk_rmsd_vs_mpnn"] = rms.rmsd(
        mpnn_u.residues[mpnn_ris].atoms.select_atoms("backbone").positions,
        af3_u.residues[af3_ris].atoms.select_atoms("backbone").positions,
        superposition=True,
    )

    row["mean_pae"] = float(af3_output.get_pae_ndarr().mean())
    return row


# def superimpose_antibody(
#     None
# )


def num_clashing_interactions_mpnn(row, af3_dir, mpnn_dir, true_dir):
    try:
        af3_output = AF3Output(af3_dir / row["assay_scaffolded_epitope_id"])
        af3_u = af3_output.get_mda_universe()
    except Exception:
        row["mpnn_clash_resindices"] = None
        row["af3_clash_resindices"] = None
        return row

    af3_ris = row["assay_scaffolded_epitope_chunk_resindices"]

    mpnn_u = mda.Universe(
        ((mpnn_dir / row["id"]) / str(row["contig_id"]))
        / (
            row["id"]
            + "_"
            + str(row["rfd_id"])
            + "_fixed_dldesign_"
            + str(row["mpnn_id"])
            + ".pdb"
        )
    )
    mpnn_ris = row["scaffolded_epitope_chunk_resindices"]

    true_u = mda.Universe(true_dir / (row["id"] + ".pdb"))
    true_ris = row["epitope_chunk_resindices"]

    alignto(
        mpnn_u.residues[mpnn_ris].atoms.select_atoms("backbone"),
        true_u.residues[true_ris].atoms.select_atoms("backbone"),
    )

    # print(mpnn_u.residues[mpnn_ris].atoms.select_atoms("name CA").positions[:10])
    # print(true_u.residues[true_ris].atoms.select_atoms("name CA").positions[:10])

    mpnn_combined_u = mda.Merge(
        mpnn_u.atoms, true_u.atoms.select_atoms("segid B or segid C")
    )
    # mpnn_combined_u.atoms.write("mpnn_combined.pdb")

    tmp_boolmask = np.full(len(mpnn_u.residues), False)
    tmp_boolmask[mpnn_ris] = True
    unintended_res = mpnn_combined_u.select_atoms("segid A").residues[~tmp_boolmask]

    row["mpnn_clash_resindices"] = list(
        mpnn_combined_u.select_atoms(
            "((segid B or segid C) and (not name H*)) and around 4 group unintended",
            unintended=unintended_res.atoms,
        ).residues.resindices
    )

    alignto(
        af3_u.residues[af3_ris].atoms.select_atoms("backbone"),
        true_u.residues[true_ris].atoms.select_atoms("backbone"),
    )
    # af3_u.atoms.chainIDs = ["D"] * len(af3_u.atoms)
    # af3_combined_u = mda.Merge(
    #     af3_u.atoms, true_u.atoms
    # )
    # af3_combined_u.atoms.write("af3_combined.pdb")
    af3_combined_u = mda.Merge(
        af3_u.atoms, true_u.atoms.select_atoms("segid B or segid C")
    )
    tmp_boolmask = np.full(len(af3_u.residues), False)
    tmp_boolmask[af3_ris] = True
    unintended_res = af3_combined_u.select_atoms("segid A").residues[~tmp_boolmask]
    # print(unintended_res)
    # print(unintended_res.atoms.select_atoms("segid A and (not name H*)"))

    l = list(
        af3_combined_u.atoms.select_atoms(
            "((segid B or segid C) and (not name H*)) and around 4 group unintended",
            unintended=unintended_res.atoms,
        ).residues.resindices
    )

    # if len(l) > 0:
    #     print("Clash dectected")

    row["af3_clash_resindices"] = l

    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_pq",
        type=str,
    )
    parser.add_argument("--true_dir", type=str)
    parser.add_argument("--mpnn_dir", type=str)
    parser.add_argument("--af3_dir", type=str)
    parser.add_argument(
        "--output_pq",
        type=str,
    )
    parser.add_argument(
        "--output_fasta",
        type=str,
    )

    args = parser.parse_args()

    mpnn_dir = Path(args.mpnn_dir)
    true_dir = Path(args.true_dir)
    af3_dir = Path(args.af3_dir)

    rmsd = pl.read_parquet(args.input_pq)

    rmsd = split_apply_combine(
        rmsd, compute_epitope_chunk_rmsd_mpnn, mpnn_dir, af3_dir, chunksize=15
    )
    rmsd = split_apply_combine(
        rmsd,
        num_clashing_interactions_mpnn,
        af3_dir,
        mpnn_dir,
        true_dir,
        chunksize=15,
    )

    rmsd.write_parquet(args.output_pq)

    valid = rmsd.filter(pl.col("mean_pae").is_not_null()).filter(
        pl.col("mean_pae") < 5,
        pl.col("overall_rmsd") <= 2,
        pl.col("epitope_chunk_rmsd_vs_mpnn") <= 1,
        pl.col("af3_clash_resindices").list.len() == 0,
    )

    if args.output_fasta:
        with open(args.output_fasta, "w") as f:
            for row in valid.iter_rows(named=True):
                f.write(f">{row['scaffolded_epitope_seq']}\n")
                f.write(row["scaffolded_epitope_seq"] + "\n")
