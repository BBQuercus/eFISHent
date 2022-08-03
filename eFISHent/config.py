"""Luigi configurations passed as command line arguments."""

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
        description=(
            "Path to output directory. "
            "If not specified, will use the current directory."
        ),
        default="",
    )


class RunConfig(luigi.Config):
    build_indices = luigi.BoolParameter(
        description="Build indices for reference genome only.", default=False
    )
    analyze_probeset = luigi.Parameter(
        description=(
            "Path to the workflow's output fasta file. "
            "Will analyze the probe set in more detail and create a PDF. "
            "Note: reference genome, target gene of interest, "
            "and (optional) count table will still have to be provided!"
        ),
        default="",
    )
    save_intermediates = luigi.BoolParameter(
        description="Save intermediate files?", default=False
    )
    optimization_method = luigi.Parameter(
        description=(
            "Optimization method to use. "
            "Greedy will procedurally select the next longest probe. "
            "Optimal will optimize probes for maximum gene coverage - "
            "slow if many overlaps are present (typically with non-stringent parameter settings). "
            "[options: optimal, greedy]."
        ),
        default="greedy",
    )
    optimization_time_limit = luigi.IntParameter(
        description=(
            "Time limit for optimization in seconds. "
            "Only used if optimization method is set to 'optimal'."
        ),
        default=60,
    )


class SequenceConfig(luigi.Config):
    sequence_file = luigi.Parameter(
        description=(
            "Path to the gene's sequence file. "
            "Use this if the sequence can't be easily downloaded or if only certain regions should be targeted."
        ),
        default="",
    )
    ensembl_id = luigi.Parameter(
        description=(
            "Ensembl ID of the gene of interest. "
            "Can be used instead of gene and organism name to download the gene of interest. "
            "Used to filter out the gene of interest from FPKM based filtering - "
            "will do an automated blast-based filtering if not passed."
        ),
        default="",
    )
    gene_name = luigi.Parameter(description="Gene name.", default="")
    organism_name = luigi.Parameter(
        description=(
            "Latin name of the organism. "
            "Can be passed together with `ensembl_id` to narrow search."
        ),
        default="",
    )
    is_plus_strand = luigi.BoolParameter(
        description=(
            "Is the probe targeting the plus strand? "
            "Note: Entrez downloads will download the gene 5'-3'. "
            "Specifying the strand is therefore not required."
        ),
        default=True,
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
        description=(
            "Minimum melting temperature in degrees C. "
            "Formamide and Na concentration will affect melting temperature!"
        ),
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
        default=10.0,
    )
    na_concentration = luigi.FloatParameter(
        description="Na concentration in mM.", default=330.0
    )
    max_off_targets = luigi.IntParameter(
        description=(
            "Maximum number of off-targets binding anywhere BUT "
            "the gene of interest in the genome."
        ),
        default=0,
    )
    encode_count_table = luigi.Parameter(
        description=(
            "Path to the ENCODE RNAseq count table provided as TSV, CSV, or TXT format. "
            "The first two columns have to be `gene_id` and `normalized_value` respectively. "
            "Specific column names are not required. "
            "`gene_id` have to be ensemble IDs of the respective genes. "
            "`normalized_value` have to be normalized count values (RPKM, FPKM, TPM, DeSeq2-norm, etc). "
        ),
        default="",
    )
    max_expression_percentage = luigi.FloatParameter(
        description=(
            "Remove any off-target genes expressed more highly than "
            "`max-expression-percentage` percentage of all genes. "
            "Genes to remove are based on ENCODE RNAseq count table provided in `encode-count-table`. "
            "Only used if `encode-count-table` and `reference-annotation` are provided."
        ),
        default=50.0,
    )
    kmer_length = luigi.IntParameter(
        description="Length of k-mer used to filter probe sequences.", default=15
    )
    max_kmers = luigi.IntParameter(
        description="Highest count of sub-k-mers found in reference genome.", default=5
    )
    max_deltag = luigi.FloatParameter(
        description=(
            "Maximum predicted deltaG in kcal/mol. "
            "Note: deltaGs range from 0 (no secondary structures) to increasingly negative values!"
        ),
        default=-10.0,
    )
    sequence_similarity = luigi.IntParameter(
        description=(
            "Maximum percentage probes can be similar. "
            "Applied to their complements and reverse complements only. "
            "Ensures probes don't falsely bind to eachother. "
            "Setting any value above 0 could add multiple minute of runtime! "
            "The lower the value (above 0) the fewer potential probes."
        ),
        default=0,
    )
