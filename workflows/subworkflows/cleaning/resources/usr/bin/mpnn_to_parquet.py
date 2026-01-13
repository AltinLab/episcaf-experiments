#!/usr/bin/env python
"""
1iqd_0P_0_fixed_dldesign_0.pdb
"""
import MDAnalysis as mda
import polars as pl
import numpy as np
from pathlib import Path
import re
import pickle
import argparse

from io import StringIO

example_json_string = """{"antigen_seq":"CSMPLGMESKAISDAQITASSYFTNMFATWSPSKARLHLQGRSNAWRPQVNNPKEWLQVDFQKTMKVTGVTTQGVKSLLTSMYVKEFLISSSQDGHQWTLFFQNGKVKVFQGNQDSFTPVVNCLDPPLLTRYLRIHPQSWVHQIALRMEVLGCEAQ","heavy_seq":"QVQLVQSGAEVKKPGASVKVSCKVSGYTLTELPVHWVRQAPGKGLEWVGSFDPESGESIYAREFQGSVTMTADTSTNIAYMELSSLRSDDTAVYYCAVPDPDAFDIWGQGTMVTVSS","light_seq":"IALTQSPGTLSLSPGERATLSCRASQSFSSSYLAWYQQKPGQAPRLLIYGASTRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQKYGTSAITFGQGTRLEIKGTV","epitope_seq":"FTNMFRRQVSLLTMHQ","id":"1iqd_0P","resolution":2.0,"light_original_segid":"A","heavy_original_segid":"B","antigen_original_segid":"C","r_work":0.203,"r_free":0.253,"antigen_ncbi_taxonomy_id":9606,"light_ncbi_taxonomy_id":9606,"heavy_ncbi_taxonomy_id":9606,"epitope_resindices":[22,23,24,25,26,41,46,48,49,76,77,78,79,81,141,142],"epitope_boolmask":[false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false,true,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,true,true,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false],"contig_string":"14-14/A23-27/23-23/A42-50/24-24/A77-82/7-7/A142-143/14-14/0 B1-117/0 C1-110/0 ","epitope_chunk_resindices":[22, 23, 24, 25, 26, 41, 42, 43, 44, 45, 46, 47, 48, 49, 76, 77, 78, 79, 80, 81, 141, 142],"epitope_chunks":4,"contig_gaps":[14,7,23,24,14],"contig_min":[0,6,10,4,0],"contig_id":0}"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdb_dir",
        type=str,
    )
    parser.add_argument("--trb_dir", type=str)
    parser.add_argument("--meta_json", type=str)
    parser.add_argument("--output_pq", type=str)
    args = parser.parse_args()

    # fake_json = StringIO(example_json_string)

    meta_df = pl.read_json(args.meta_json)

    pdb_pat = re.compile(r"_(\d+)_fixed.*?dldesign_(\d+)\.pdb$", re.IGNORECASE)
    rows = []

    for f in Path(args.pdb_dir).iterdir():
        if f.is_file() and f.suffix.lower() == ".pdb":
            m = pdb_pat.search(f.name)
            if not m:
                raise ValueError

            rfdiffusion_id = int(m.group(1))
            mpnn_id = int(m.group(2))

            sequence = mda.Universe(f).segments[0].residues.sequence(format="string")

            rows.append(
                {
                    "rfd_id": rfdiffusion_id,
                    "mpnn_id": mpnn_id,
                    "scaffolded_epitope_seq": sequence,
                }
            )

    mpnn_df = pl.DataFrame(rows).with_columns(
        pl.lit(meta_df.select("id").item()).alias("id")
    )

    # trb_pat = re.compile(r"_(\d+)\.trb$", re.IGNORECASE)
    # trb_df = None

    # for f in Path(args.trb_dir).iterdir():
    #     if f.is_file() and f.suffix.lower() == ".trb":
    #         m = trb_pat.search(f.name)

    #         if not m:
    #             raise ValueError

    #         rfdiffusion_id = int(m.group(1))

    #         # snippet modified from: https://github.com/nrbennet/dl_binder_design/blob/main/helper_scripts/addFIXEDlabels.py
    #         data = np.load(f, allow_pickle=True)

    #         if "receptor_con_hal_pdb_idx" in data:
    #             # Identify the last residue number in A chain
    #             last_res_id = int(data["receptor_con_hal_pdb_idx"][0][1]) - 1

    #             # Identify where inpaint seq is True, ie, kept fixed
    #             indices = np.where(data["inpaint_seq"][:last_res_id])[0]
    #         else:
    #             # Identify where inpaint seq is True in the only chain
    #             indices = np.where(data["inpaint_seq"])[0]

    #         offset = (
    #             indices[0] - meta_df.select("epitope_chunk_resindices")[0].item()[0]
    #         )

    #         trb_df = pl.DataFrame(
    #             [
    #                 {
    #                     "scaffolded_epitope_chunk_resindices": list(indices),
    #                     "scaffolded_epitope_resindices": meta_df.select(
    #                         "epitope_resindices"
    #                     )[0].item()
    #                     + offset,
    #                 }
    #             ]
    #         ).with_columns(pl.lit(meta_df.select("id").item()).alias("id"))
    #         # all TRB files are the same for the same contig string
    #         break

    meta_df = meta_df.join(mpnn_df, on="id")

    # .join(trb_df, on="id")

    meta_df.write_parquet(args.output_pq)
