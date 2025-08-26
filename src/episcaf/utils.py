import polars as pl
from Bio import SeqIO
import polars as pl


def fasta_to_polars(fasta_path: str, desc_as_name: bool = False) -> pl.DataFrame:
    """
    Read a FASTA file and convert it into a Polars DataFrame
    with columns ["name", "seq"].

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    pl.DataFrame
        - "name": the sequence ID
        - "seq": the full sequence string
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))

    if desc_as_name:
        names = [rec.description for rec in records]
    else:
        names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]

    df = pl.DataFrame({"name": names, "seq": seqs})
    return df
