
process RFDIFFUSION {
    label "rfdiffusion"

    input:
    tuple val(meta), path(pdb_file)
    val rfd_num_designs

    output:
    tuple val(meta), path("*.pdb"), path("*.trb")


    script:
    """
    python3.9 /app/RFdiffusion/scripts/run_inference.py \\
        inference.model_directory_path=${params.rfdiffusion_models_path} \\
        inference.schedule_directory_path=${params.rfdiffusion_schedule_path} \\
        inference.input_pdb="${pdb_file}" \\
        'contigmap.length=${meta.contig_length}' \\
        'contigmap.contigs=[$meta.contig_string]' \\
        inference.num_designs=${rfd_num_designs} \\
        inference.output_prefix=./${meta.id}
    """
}

process ADD_FIXED_LABELS {
    label "proteinmpnn_local"

    input:
    tuple val(meta), path(pdb), path(trb)

    output:
    tuple val(meta), path("*_fixed.pdb"), path(trb)

    script:
    """
    addFIXEDlabels.py \\
        --pdbdir . \\
        --trbdir . \\
        --outdir . \\
        --verbose 
    """
}

process DL_INTERFACE_DESIGN {
    label "proteinmpnn_gpu"

    input:
    tuple val(meta), path(pdb), path(trb)
    val mpnn_seqs_per_struct

    output:
    tuple val(meta), path("*.pdb"), path(trb)

    script:
    """
    ${params.dl_binder_design_path}/mpnn_fr/dl_interface_design.py \\
        -pdbdir . \\
        -outpdbdir . \\
        -seqs_per_struct ${mpnn_seqs_per_struct} \\
        -temperature 0.1 \\
        -relax_cycles 0
    """
}