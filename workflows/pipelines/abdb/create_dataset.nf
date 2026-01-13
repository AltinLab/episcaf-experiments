nextflow.preview.output = true

include { FILTER_ANTIBODIES; EPITOPE_ANNOTATION } from '../../subworkflows/cleaning'
include { CLUSTER_COMPLEXES } from '../../subworkflows/clustering'


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

    // add step to shorten PDB files to CDR region only for RFdiffusion speed

    publish:
    candidate_complex = candidate_complex
    pdbfiles = pdbfiles

    heavy_tsv = heavy_tsv
    light_tsv = light_tsv
    antigen_tsv = antigen_tsv
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
    chosen_complex {
        path "complexes"
    }
    discard_complex {
        path "complexes"
    }
}