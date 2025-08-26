

include { PARQUET_TO_FASTA; 
            ANNOTATE_REPRESENTATIVES;
            CLUSTER_FASTA as CLUSTER_ANTIGEN_FASTA;
            CLUSTER_FASTA as CLUSTER_HEAVY_FASTA;
            CLUSTER_FASTA as CLUSTER_LIGHT_FASTA } from '../../../modules/local/clustering'

workflow CLUSTER_COMPLEXES {
    take:
    candidate_complex

    main:
    PARQUET_TO_FASTA(candidate_complex)
    antigen_tsv = CLUSTER_ANTIGEN_FASTA(PARQUET_TO_FASTA.out.antigen_fasta, 0.8)
    heavy_tsv = CLUSTER_HEAVY_FASTA(PARQUET_TO_FASTA.out.heavy_fasta, 0.95)
    light_tsv = CLUSTER_LIGHT_FASTA(PARQUET_TO_FASTA.out.light_fasta, 0.95)
    ANNOTATE_REPRESENTATIVES(candidate_complex, antigen_tsv, heavy_tsv, light_tsv)
    chosen_complex = ANNOTATE_REPRESENTATIVES.out.chosen_complex
    discard_complex = ANNOTATE_REPRESENTATIVES.out.discard_complex

    emit:
    chosen_complex = chosen_complex
    discard_complex = discard_complex
    heavy_tsv = heavy_tsv
    light_tsv = light_tsv
    antigen_tsv = antigen_tsv
}