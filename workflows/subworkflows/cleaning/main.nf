

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

process META_AS_JSON {
    label "episcaf_local"

    input:
    tuple val(meta), path(pdb), path(trb)

    output:
    tuple val(meta), path(pdb), path(trb), path("*.json")

    script:
    def meta_as_json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(meta))
    """
    cat <<'JSON' > "${meta.id}_${meta.contig_id}.json" 
${meta_as_json_string}
JSON
    """
}

process MPNN_TO_PARQUET {

    label "episcaf_local"

    input:
    tuple val(meta), path(pdb), path(trb), path(json)

    output:
    path "*.parquet"

    script: 
    """
    mpnn_to_parquet.py \\
        --pdb_dir . \\
        --trb_dir . \\
        --meta_json $json \\
        --output_pq "${meta.id}_${meta.contig_id}.parquet"
    """
}

process CONCAT_PARQUETS {

    label "episcaf_local"

    input:
    path pq

    output:
    path "*.parquet"

    script:
    """
    #!/usr/bin/env python

    from pathlib import Path
    import polars as pl

    dfs = []

    for f in Path(".").iterdir():
        if f.is_file() and f.suffix.lower() == ".parquet":
            dfs.append(pl.read_parquet(f))

    out_df = pl.concat(dfs)

    out_df.write_parquet("scaffold_summary.parquet")
    """

}