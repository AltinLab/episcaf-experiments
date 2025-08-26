nextflow.preview.output = true

include { FILTER_ANTIBODIES; EPITOPE_ANNOTATION } from '../subworkflows/local/cleaning'
include { CLUSTER_COMPLEXES } from '../subworkflows/local/clustering'


workflow {

    main:
    input_dir = Channel.fromPath(params.input)

    FILTER_ANTIBODIES(input_dir)
    candidate_complex = FILTER_ANTIBODIES.out.candidate_complex
    pdbfiles = FILTER_ANTIBODIES.out.pdbfiles

    candidate_complex = EPITOPE_ANNOTATION(candidate_complex, pdbfiles)

    CLUSTER_COMPLEXES(candidate_complex)

    heavy_tsv = CLUSTER_COMPLEXES.out.heavy_tsv
    light_tsv = CLUSTER_COMPLEXES.out.light_tsv
    antigen_tsv = CLUSTER_COMPLEXES.out.antigen_tsv
    
    chosen_complex = CLUSTER_COMPLEXES.out.chosen_complex
    discard_complex = CLUSTER_COMPLEXES.out.discard_complex


    // cleaned_pdbfiles = STANDARDIZE_PDBFILES(chosen_complex, pdbfiles)

    // chosen_complex_tojoin = ANNOTATE_EPITOPES.out.chosen_complex_annot.splitParquet().map {
    //     row -> 

    //         def id = row["complexed_antibody_id"]
    //         def meta = [
    //             id: id,
    //             epitope_positions: row["epitope_positions"],
    //             light_seq: row["light_seq"],
    //             heavy_seq: row["heavy_seq"],
    //             antigen_seq: row["antigen_seq"],
    //         ]

    //         tuple(id, meta)

    // } 

    // cleaned_pdbfiles_tojoin = cleaned_pdbfiles.map { f ->

    //     def id = f.getSimpleName()
    //     tuple(id, f)
    // }

    // meta_pdb = chosen_complex_tojoin.join(cleaned_pdbfiles_tojoin).map { key, meta, pdb_file ->
    //     tuple(meta, pdb_file)
    // }

    



    // map each complex row into meta, PDB format

    publish:
    candidate_complex = candidate_complex
    pdbfiles = pdbfiles

    heavy_tsv = heavy_tsv
    light_tsv = light_tsv
    antigen_tsv = antigen_tsv
    // cleaned_pdbfiles = cleaned_pdbfiles
    chosen_complex = chosen_complex
    discard_complex = discard_complex
}

output {
    candidate_complex {
        path "complexes"
    }
    pdbfiles {
        path "complex_pdbfiles/cleaned"
    }

    heavy_tsv {
        path "complexes"
    }
    light_tsv {
        path "complexes"
    }
    antigen_tsv {
        path "complexes"
    }
    // cleaned_pdbfiles {
    //     path "complex_pdbfiles/cleaned"
    // }
    chosen_complex {
        path "complexes"
    }
    discard_complex {
        path "complexes"
    }
}