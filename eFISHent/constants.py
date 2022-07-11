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

SAM_FLAG_REVERSE = "16"

CLI_SHORTFORM = {
    "reference_genome": "g",
    "reference_annotation": "a",
    "threads": "t",
    "output_dir": "o",
    "build_indices": "idx",
    "save_intermediates": "int",
    "optimization_method": "om",
    "optimization_time_limit": "otl",
    "sequence_file": "seq",
    "ensembl_id": "id",
    "gene_name": "gen",
    "organism_name": "org",
    "is_plus_strand": "plus",
    "is_endogenous": "endo",
    "min_length": "l",
    "max_length": "L",
    "spacing": "s",
    "min_tm": "tm",
    "max_tm": "TM",
    "min_gc": "gc",
    "max_gc": "GC",
    "formamide_concentration": "f",
    "na_concentration": "na",
    "max_off_targets": "ot",
    "encode_count_table": "ct",
    "max_expression_percentage": "ep",
    "kmer_length": "kl",
    "max_kmers": "km",
    "max_deltag": "dg",
}
