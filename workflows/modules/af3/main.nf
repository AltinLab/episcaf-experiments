process SEQ_LIST_TO_FASTA {
    label "af3_process_local"

    input:
      tuple val(meta), val(seq_list)

    output:
      tuple val(meta), path("${meta.id}.fasta")

    script:
    """
    filename="${meta.id}.fasta"
    : > "\$filename"

    i=1
    for seq in ${seq_list.join(' ')}; do
        echo ">\$i"   >> "\$filename"
        echo "\$seq"  >> "\$filename"
        i=\$(( i + 1 ))
    done
    """
}

process SEQ_LIST_TO_FASTA_FILT {
    label "af3_process_local"

    input:
      tuple val(meta), val(seq_list)
      val(inf_dir)

    output:
      tuple val(meta), path("${meta.id}.fasta"), optional: true

    script:
    """
    if [ ! -e "${inf_dir}/${meta.id}" ]; then
        filename="${meta.id}.fasta"
        : > "\$filename"

        i=1
        for seq in ${seq_list.join(' ')}; do
            echo ">\$i"   >> "\$filename"
            echo "\$seq"  >> "\$filename"
            i=\$(( i + 1 ))
        done
    fi
    """
}

process COMPOSE_INFERENCE_JSON {
    label "af3_process_local"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json")


    script:
    """
    compose_inference_JSON.py \\
        --job_name "${meta.id}" \\
        --fasta "$fasta" \\
        --seeds 1
    """
 }

process INFERENCE {
    tag "inference"
    label "alphafold3_inference"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*")

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --json_path=$json \\
        --model_dir=$params.af3_model_dir \\
        --db_dir=$params.af3_db_dir \\
        --output_dir=. \\
        --norun_data_pipeline \\
        --num_diffusion_samples=1
    """
}

process CLEAN_INFERENCE_DIR {
    label "af3_process_local"
    
    input:
    tuple val(meta), path(inference_dir)

    output:
    tuple val(meta), path("*", includeInputs: true)

    script:
    """
    clean_inference_dir.py \\
        -i $inference_dir
    """
}
