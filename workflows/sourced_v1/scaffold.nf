nextflow.preview.output = true

include { GENERATE_CONTIGS; RFDIFFUSION; ADD_FIXED_LABELS; DL_INTERFACE_DESIGN } from '../subworkflows/local/rfdiffusion'
include { splitParquet } from 'plugin/nf-parquet'

def hash_from(String seq, String algorithm) {
    def md = java.security.MessageDigest.getInstance(algorithm)
    md.update(seq.getBytes('UTF-8'))
    def bytes = md.digest()
    bytes.collect { String.format('%02x', it) }.join()
}


workflow {

    main:

    input_pq = Channel.fromPath(params.input)
    input_pdb = Channel.fromPath(params.pdb_dir + "/*.pdb").toList()

    contig_pq = GENERATE_CONTIGS(input_pq, input_pdb)

    chosen_complex_tojoin = contig_pq.splitParquet().map {
        row -> 

            def job_name = hash_from(row["antigen_seq"], "MD5")
            def msa_name = hash_from(row["antigen_seq"], "SHA-256")
            def id = row["complexed_antibody_id"]
            def meta = [
                id: id,
                contig: row["contig_string"],
                contig_length: "104-104",
                job_name: job_name,
                msa_name: msa_name,
                // epitope_positions: row["epitope_positions"],
                // light_seq: row["light_seq"],
                // heavy_seq: row["heavy_seq"],
                // antigen_seq: row["antigen_seq"],
            ]

            tuple(id, meta)

    } 

    cleaned_pdbfiles_tojoin = input_pdb.flatten().map { f ->
        def id = f.getSimpleName()
        tuple(id, f)
    }

    meta_pdb = chosen_complex_tojoin.join(cleaned_pdbfiles_tojoin).map { key, meta, pdb_file ->
        tuple(meta, pdb_file)
    }

    // meta_pdb.view()

    meta_diff_pdbs_trbs = RFDIFFUSION(meta_pdb)

    meta_diff_pdbs_trbs = ADD_FIXED_LABELS(meta_diff_pdbs_trbs)

    meta_dl_pdbs_trbs = DL_INTERFACE_DESIGN(meta_diff_pdbs_trbs)

    // PDB_TO_SCAFFOLD_FASTA(meta_diff_pdbs_trbs)

    // MSA_WORKFLOW(meta_diff_pdbs)

    // INFERENCE_WORKFLOW(meta_diff_pdbs)

    publish:
    contig_pq = contig_pq
    meta_diff_pdbs_trbs = meta_diff_pdbs_trbs
    meta_dl_pdbs_trbs = meta_dl_pdbs_trbs
}

output {
    contig_pq {}
    meta_diff_pdbs_trbs { path "rf" }
    meta_dl_pdbs_trbs { path "dl" }
}

