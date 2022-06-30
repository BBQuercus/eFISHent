"""Common constants used in other methods."""
from typing import List
import luigi

from .config import GeneralConfig
from .config import ProbeConfig
from .config import RunConfig
from .config import SequenceConfig

SAMFILE_COLUMNS = [
    "qname",
    "flag",
    "rname",
    "pos",
    "mapq",
    "cigar",
    "rnext",
    "pnext",
    "tlen",
    "seq",
    "qual",
    "XA",
    "MD",
    "NM",
]

BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend ",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

GTF_COLUMNS = [
    "gene_id",
    "seqname",
    "start",
    "end",
    "frame",
    "exon_number",
    "feature",
]

CONFIG_CLASSES: List[luigi.Config] = [
    GeneralConfig,
    RunConfig,
    SequenceConfig,
    ProbeConfig,
]

FASTA_EXT = (".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa")
