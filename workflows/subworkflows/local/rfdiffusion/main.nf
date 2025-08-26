
process GENERATE_CONTIGS {
    label "episcaf_local"

    input:
    path complex_pq
    path pdb_dir

    output:
    path "*.parquet"

    script:
    """
    generate_contigs.py \\
        --complex_pq ${complex_pq} \\
        --topology_path . \\
        --output_pq "${complex_pq.getSimpleName()}.contigs.parquet" \\
    """
}



process RFDIFFUSION {
    label "rfdiffusion"

    input:
    tuple val(meta), path(pdb_file)

    output:
    tuple val(meta), path("*.pdb"), path("*.trb")


    script:
    """
    python3.9 /app/RFdiffusion/scripts/run_inference.py \\
        inference.model_directory_path=${params.rfdiffusion_models_path} \\
        inference.schedule_directory_path=${params.rfdiffusion_schedule_path} \\
        inference.input_pdb="${pdb_file}" \\
        'contigmap.length=${meta.contig_length}' \\
        'contigmap.contigs=[$meta.contig]' \\
        inference.output_prefix=./${meta.id}
    """
}

process ADD_FIXED_LABELS {
    executor 'local'
    conda "${moduleDir}/environment.yaml"

    input:
    tuple val(meta), path(pdb), path(trb)

    output:
    tuple val(meta), path("*.pdb"), path(trb)

    script:
    """
    addFIXEDlabels.py \\
        --pdbdir ${pdb} \\
        --trbdir ${trb} \\
        --verbose 
    """
}

process DL_INTERFACE_DESIGN {
    executor 'local'
    conda "${moduleDir}/environment.yaml"

    input:
    tuple val(meta), path(pdb), path(trb)

    output:
    tuple val(meta), path("*.pdb"), path(trb)

    script:
    """
    ${params.df_binder_design_path}/mpnn_fr/dl_interface_design.py \\
        -pdbdir ${pdb} \\
        -outpdbdir . \\
        -seqs_per_struct 8 \\
        -relax_cycles 0
    """
}