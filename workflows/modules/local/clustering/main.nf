process PARQUET_TO_FASTA {
  label "episcaf_local"
  
  input:
      path(parquet)
  
  output:
      path("*antigen.fasta"), emit: antigen_fasta
      path("*heavy.fasta"), emit: heavy_fasta
      path("*light.fasta"), emit: light_fasta
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  
  df = pl.read_parquet("${parquet}")
  
  with open("${parquet.getSimpleName()}_antigen.fasta", "w") as f:
      for row in df.iter_rows(named=True):
          f.write(f">{row['complexed_antibody_id']}\\n{row['antigen_seq']}\\n")

  with open("${parquet.getSimpleName()}_heavy.fasta", "w") as f:
      for row in df.iter_rows(named=True):
          f.write(f">{row['complexed_antibody_id']}\\n{row['heavy_seq']}\\n")

  with open("${parquet.getSimpleName()}_light.fasta", "w") as f:
    for row in df.iter_rows(named=True):
        f.write(f">{row['complexed_antibody_id']}\\n{row['light_seq']}\\n")
  """
}

process FASTA_TO_PARQUET {
  label "episcaf_local"

  
  input:
      path(fasta)
  
  output:
      path("*.parquet")
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  from episcaf.utils import fasta_to_polars
  
  df = fasta_to_polars("${fasta}").rename({"name" : "complexed_antibody_id", "seq" : "antigen_seq"})
  
  df.write_parquet("${fasta.getSimpleName()}.parquet")
  """
}


process CLUSTER_FASTA {
    conda "${moduleDir}/mmseqs2.yaml"
    executor 'local'
    publishDir "/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/tmp/DEBUG_CLUST", mode: 'copy'

    input:
        path(fasta)
        val(min_seq_id)

    output:
        path("*.tsv")

    script:
    """
    mmseqs createdb ${fasta} db
    mmseqs cluster db db_clu tmp \\
        --min-seq-id ${min_seq_id}
    mmseqs createtsv db db db_clu ${fasta.getSimpleName()}.tsv
    """
}


process ANNOTATE_REPRESENTATIVES {
    label "episcaf_local"

    input:
        path(orig_pq)
        path(antigen_rep_tsv)
        path(heavy_rep_tsv)
        path(light_rep_tsv)

    output:
        path("*chosen*.parquet"), emit: chosen_complex
        path("*discard*.parquet"), emit: discard_complex

    script:
    """
    pick_reps.py \\
        --complex_pq ${orig_pq} \\
        --antigen_rep_tsv ${antigen_rep_tsv} \\
        --heavy_rep_tsv ${heavy_rep_tsv} \\
        --light_rep_tsv ${light_rep_tsv} \\
        --output_chosen "${orig_pq.getSimpleName()}.chosen.parquet" \\
        --output_discard "${orig_pq.getSimpleName()}.discard.parquet"
    """
}