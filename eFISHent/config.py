"""
Processes:
    - index creation
    - probe generation

General options:
    - reference genome -> pass file, will check if index exists in same folder
    - reference genome annotation
    - ENCODE RNAseq count table for off-target weighting
    - number of threads
    - method of optimization
    - optimization time limit

Probe specific options:
    - (sequence file) or (gene name & organism)
    - strand on genome
    - intronic/exonic/both
    - endogenous or exogenous
    - min and max length
    - min and max melting temperature
    - min and max GC content
    - formamide concentration
    - salt concentration
    - maximum deltaG
"""

import luigi


class GeneralConfig(luigi.Config):
    reference_genome = luigi.Parameter(
        description="Path to reference genome fasta/fa file.",
        significant=True,
        default="",
    )
    reference_annotation = luigi.Parameter(
        description="Path to reference genome annotation (gene transfer format / gtf) file.",
        default="",
    )
    threads = luigi.IntParameter(description="Number of threads to launch.", default=2)
    output_dir = luigi.Parameter(
        description="Path to output directory. "
        "If not specified, will use the current directory.",
        default="",
    )


class RunConfig(luigi.Config):
    build_indices = luigi.BoolParameter(
        description="Build indices for reference genome only.", default=False
    )
    save_intermediates = luigi.BoolParameter(
        description="Save intermediate files?", default=False
    )
    optimization_method = luigi.Parameter(
        description="Optimization method to use. "
        "Greedy will procedurally select the next longest probe. "
        "Optimal will optimize probes for maximum gene coverage - "
        "slow if many overlaps are present (typically with non-stringent parameter settings). "
        "[options: optimal, greedy].",
        default="greedy",
    )
    optimization_time_limit = luigi.IntParameter(
        description="Time limit for optimization in seconds. "
        "Only used if optimization method is set to 'optimal'.",
        default=60,
    )


class SequenceConfig(luigi.Config):
    sequence_file = luigi.Parameter(
        description="Path to the gene's sequence file. "
        "Use this if the sequence can't be easily downloaded or if only certain regions should be targeted.",
        default="",
    )
    ensembl_id = luigi.Parameter(
        description="Ensembl ID of the gene of interest."
        "Can be used instead of gene and organism name to download the gene of interest."
        "Used to filter out the gene of interest from FPKM based filtering - "
        "will do an automated blast-based filtering if not passed.",
        default="",
    )
    gene_name = luigi.Parameter(description="Gene name.", default="")
    organism_name = luigi.Parameter(
        description="Latin name of the organism.", default=""
    )
    is_plus_strand = luigi.BoolParameter(
        description="Is the probe targeting the plus strand?", default=True
    )
    is_endogenous = luigi.BoolParameter(
        description="Is the probe endogenous?", default=True
    )


class ProbeConfig(luigi.Config):
    min_length = luigi.IntParameter(
        description="Minimum probe length in nucleotides.", default=21
    )
    max_length = luigi.IntParameter(
        description="Maximum probe length in nucleotides.", default=25
    )
    spacing = luigi.IntParameter(
        description="Minimum distance in nucleotides between probes.", default=2
    )
    min_tm = luigi.FloatParameter(
        description="Minimum melting temperature in degrees C. "
        "Formamide and Na concentration will affect melting temperature!",
        default=40.0,
    )
    max_tm = luigi.FloatParameter(
        description="Maximum melting temperature in degrees C (see minimum).",
        default=60.0,
    )
    min_gc = luigi.FloatParameter(
        description="Minimum GC content in percentage.", default=20.0
    )
    max_gc = luigi.FloatParameter(
        description="Maximum GC content in percentage.", default=80.0
    )
    formamide_concentration = luigi.FloatParameter(
        description="Formamide concentration as a percentage of the total buffer.",
        default=0.0,
    )
    na_concentration = luigi.FloatParameter(
        description="Na concentration in mM.", default=390.0
    )
    max_off_targets = luigi.IntParameter(
        description="Maximum number of off-targets binding anywhere BUT "
        "the gene of interest in the genome.",
        default=0,
    )
    encode_count_table = luigi.Parameter(
        description="Path to the ENCODE RNAseq count table. "
        "Must contain the columns 'gene_id' and 'FPKM'. "
        "'gene_id' must use ensembl ID's!",
        default="",
    )
    max_off_target_fpkm = luigi.FloatParameter(
        description=(
            "Maximum off-target FPKM (Fragments Per Kilobase of transcript per Million mapped reads) "
            "based on ENCODE RNAseq count table. "
            "Only used if ENCODE RNAseq count table is provided."
        ),
        default=10.0,
    )
    kmer_length = luigi.IntParameter(
        description="Length of k-mer used to filter probe sequences.", default=15
    )
    max_kmers = luigi.IntParameter(
        description="Highest count of sub-k-mers found in reference genome.", default=5
    )
    max_deltag = luigi.FloatParameter(
        description="Maximum predicted deltaG in kcal/mol. "
        "Note: deltaGs range from 0 (no secondary structures) to increasingly negative values!",
        default=-10.0,
    )
