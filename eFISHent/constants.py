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
    "analyze_probeset": "ap",
    "build_indices": "idx",
    "encode_count_table": "ct",
    "ensembl_id": "id",
    "formamide_concentration": "f",
    "gene_name": "gen",
    "is_endogenous": "endo",
    "is_plus_strand": "plus",
    "kmer_length": "kl",
    "max_deltag": "dg",
    "max_expression_percentage": "ep",
    "max_gc": "GC",
    "max_kmers": "km",
    "max_length": "L",
    "max_off_targets": "ot",
    "max_tm": "TM",
    "min_gc": "gc",
    "min_length": "l",
    "min_tm": "tm",
    "na_concentration": "na",
    "optimization_method": "om",
    "optimization_time_limit": "otl",
    "organism_name": "org",
    "output_dir": "o",
    "reference_annotation": "a",
    "reference_genome": "g",
    "save_intermediates": "int",
    "sequence_file": "seq",
    "sequence_similarity": "sim",
    "spacing": "sp",
    "threads": "t",
}
