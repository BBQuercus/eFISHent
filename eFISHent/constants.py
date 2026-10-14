"""Common constants used in other methods."""

import luigi

from .config import GeneralConfig, ProbeConfig, RunConfig, SequenceConfig

SAMFILE_COLUMNS: list[str] = [
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

BLAST_COLUMNS: list[str] = [
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

GTF_COLUMNS: list[str] = [
    "gene_id",
    "seqname",
    "start",
    "end",
    "frame",
    "exon_number",
    "feature",
]

# GTF rows and attributes kept in the prepared parquet annotation. Downstream
# code only reads gene, transcript and exon rows; CDS/UTR/codon rows would
# double the row count for nothing.
GTF_FEATURES: list[str] = ["gene", "transcript", "exon"]
GTF_ATTRIBUTES: list[str] = [
    "gene_id",
    "gene_name",
    "gene_biotype",  # Ensembl
    "gene_type",  # GENCODE
    "transcript_id",
    "exon_number",
]

CONFIG_CLASSES: list[luigi.Config] = [
    GeneralConfig,
    RunConfig,
    SequenceConfig,
    ProbeConfig,
]

FASTA_EXT: tuple[str] = (".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa")

SAM_FLAG_REVERSE: str = "16"

CLI_SHORTFORM: dict[str, str] = {
    "aligner": "al",
    "analyze_probeset": "ap",
    "blast_identity_threshold": "bit",
    "build_indices": "idx",
    "encode_count_table": "ct",
    "ensembl_id": "id",
    "filter_g_quadruplex": "fgq",
    "filter_low_complexity": "flc",
    "formamide_concentration": "f",
    "gene_name": "gen",
    "is_endogenous": "endo",
    "is_plus_strand": "plus",
    "kmer_length": "kl",
    "max_deltag": "dg",
    "max_expression_percentage": "ep",
    "max_gc": "GC",
    "max_homopolymer_length": "hp",
    "max_kmers": "km",
    "max_length": "L",
    "max_off_targets": "ot",
    "max_tm": "TM",
    "max_transcriptome_off_targets": "txot",
    "min_gc": "gc",
    "min_length": "l",
    "min_tm": "tm",
    "na_concentration": "na",
    "no_alternative_loci": "nal",
    "optimization_method": "om",
    "optimization_time_limit": "otl",
    "organism_name": "org",
    "output_dir": "o",
    "reference_annotation": "a",
    "reference_genome": "g",
    "reference_transcriptome": "rt",
    "save_intermediates": "int",
    "sequence_file": "seq",
    "sequence_similarity": "sim",
    "spacing": "sp",
    "threads": "t",
    "filter_rrna": "rrna",
    "intergenic_off_targets": "iot",
    "mask_repeats": "mr",
    "off_target_min_tm": "ottm",
    "min_blast_match_length": "bml",
    "max_probes_per_off_target": "mpot",
    "adaptive_length": "adl",
    "filter_rdna_45s": "rdna45",
    "custom_rdna_fasta": "rdnaf",
    "reject_cross_hybridization": "rxh",
    "allow_no_transcriptome": "ant",
    "max_cpg_fraction": "cpg",
    "accessibility_scoring": "acc",
    "target_regions": "tr",
}
