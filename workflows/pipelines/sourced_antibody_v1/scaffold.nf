nextflow.preview.output = true

/*
Parameters:
- input: parquet file with at least the columns:
  - id: identifier that corresponds to a <id>.pdb in params.pdb_dir
  - contig_id: integer contig number for row. in combination with id, creates unique ID for row
  - contig_string: RFdiffusion-formatted contigmap.contigs string
  - contig_length: RFdiffusion-formatted contigmap.length string
- rfd_num_designs: number of rfdiffusion designs per id:contig_id
- mpnn_seqs_per_struct: number of mpnn predicted sequences per id:contig_id:rfd_id
- pdb_dir: directory containing a pdb file for each id in input parquet

Outputs:
- summary_scaffold.parquet: a dataframe with all columns in input as well as:
  - rfd_id: integer rfd design number for id:contig_id combinantion
  - mpnn_id: integer protein mpnn design number for id:contig_id:rfd_id combination
  - scaffolded_epitope_seq: sequence of residues of the rfd+mpnn predicted scaffold (including the epitope chunks)

Each row in input will receive rfd_num_designs * mpnn_seqs_per_struct rows/designs in output
*/

include{ META_AS_JSON; MPNN_TO_PARQUET; CONCAT_PARQUETS } from '../../subworkflows/cleaning'
include { RFDIFFUSION; ADD_FIXED_LABELS; DL_INTERFACE_DESIGN } from '../../subworkflows/rfdiffusion'
include { splitParquet } from 'plugin/nf-parquet'




workflow {

    main:

    contig_pq = Channel.fromPath(params.input)
    input_pdb = Channel.fromPath(params.pdb_dir + "/*.pdb").toList()

    chosen_complex_tojoin = contig_pq.splitParquet().map {
        row -> 

            def id = row["id"]

            def meta = new LinkedHashMap(row)

            // println meta?.class
            // println groovy.json.JsonOutput.toJson(meta)

            tuple(id, meta)

    } 

    // chosen_complex_tojoin.view()

    cleaned_pdbfiles_tojoin = input_pdb.flatten().map { f ->
        def id = f.getSimpleName()
        tuple(id, f)
    }

    meta_pdb = chosen_complex_tojoin.combine(cleaned_pdbfiles_tojoin, by: 0).map { key, meta, pdb_file ->
        tuple(meta, pdb_file)
    }

    // meta_pdb.view()

    meta_diff_pdbs_trbs = RFDIFFUSION(meta_pdb, params.rfd_num_designs)

    meta_diff_pdbs_trbs_labelled = ADD_FIXED_LABELS(meta_diff_pdbs_trbs)

    meta_dl_pdbs_trbs = DL_INTERFACE_DESIGN(meta_diff_pdbs_trbs_labelled, params.mpnn_seqs_per_struct)

    meta_dl_pdbs_trbs_metajson = META_AS_JSON(meta_dl_pdbs_trbs)

    mpnn_pq = MPNN_TO_PARQUET(meta_dl_pdbs_trbs_metajson)

    output_pq = CONCAT_PARQUETS(mpnn_pq.toList())

    publish:
    contig_pq = contig_pq
    // meta_diff_pdbs_trbs = meta_diff_pdbs_trbs
    meta_dl_pdbs_trbs = meta_dl_pdbs_trbs
    output_pq = output_pq
}

output {
    contig_pq {}
    // meta_diff_pdbs_trbs { 
    //     path {meta, pdb, trb -> 
    //         "rfdiffusion/${meta.id}/${meta.contig_id}" 
    //     }
    // }
    meta_dl_pdbs_trbs { 
        path {meta, pdb, trb -> 
            "proteinmpnn/${meta.id}/${meta.contig_id}" 
        }
    }
    output_pq {}


}

