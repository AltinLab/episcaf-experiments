

process FILTER_ANTIBODIES {
    label "episcaf_local"

    input:
    path input_dir

    output:
    path "*.parquet", emit: candidate_complex
    path "*.pdb", emit: pdbfiles

    script:
    """
    filt_abdb.py \\
        --abdb_snapshot_dir ${input_dir} \\
        --output_parquet candidate_complex.parquet \\
        --output_pdb_dir .
    """
}

// process STANDARDIZE_PDBFILES {
//     label "episcaf_local"

//     input:
//     path chosen_complex
//     path pdbfiles

//     output:
//     path "*.pdb"

//     script:
//     """
    
//     standardize_pdbfiles.py \\
//         --complex_pq ${chosen_complex} \\
//         --topology_path . \\
//         --output_dir .
//     """
// }

process EPITOPE_ANNOTATION {
    label "episcaf_local"

    input:
    path chosen_complex
    path cleaned_pdbfiles

    output:
    path "*annot*.parquet"

    script:
    """
    annot_epitopes.py \\
        --complex_pq ${chosen_complex} \\
        --topology_path . \\
        --output_pq "${chosen_complex.getSimpleName()}.annot.parquet"
    """
}