nextflow.preview.output = true

/*
Parameters:
- input: parquet file with at least the columns:
  - assay_scaffolded_epitope_id
  - assay_scaffolded_epitope_seq

Outputs:
- $outputDir/af3_predictions: Directory containing a subdirectory <assay_scaffolded_epitope_id> for each
  row in input
*/

include { SEQ_LIST_TO_FASTA_FILT; COMPOSE_INFERENCE_JSON; INFERENCE; CLEAN_INFERENCE_DIR } from '../../modules/af3'
include { splitParquet } from 'plugin/nf-parquet'

workflow {

    main:

    scaffold_pq = Channel.fromPath(params.input)
    inf_dir = file("$params.inf_dir").toUriString()

    scaffold_ch = scaffold_pq.splitParquet().map {
        row -> 

            tuple([
                id: row["assay_scaffolded_epitope_id"]
            ], 
            [row["assay_scaffolded_epitope_seq"]])
    } 

    meta_fasta = SEQ_LIST_TO_FASTA_FILT(scaffold_ch, inf_dir)
    meta_json = COMPOSE_INFERENCE_JSON(meta_fasta)
    meta_inf = INFERENCE(meta_json)
    meta_inf = CLEAN_INFERENCE_DIR(meta_inf)


    publish:
    meta_inf = meta_inf
}

output {

    meta_inf { path "af3_predictions" }

}

