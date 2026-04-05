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
    reference_transcriptome = luigi.Parameter(
        description=(
            "Path to reference transcriptome FASTA file. "
            "Enables transcriptome-level off-target filtering via BLAST. "
            "Generate from genome + GTF using gffread."
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
    no_alternative_loci = luigi.BoolParameter(
        description=(
            "When aligning do not consider alternative loci. "
            "These are all _alt labelled sequences in the reference genome."
        ),
        default=False,
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
            "Maximum percentage probes can be similar to each other's "
            "complement/reverse complement. Prevents probe-probe cross-hybridization. "
            "Set to 0 to disable (default). "
            "Recommended: 75-85 for stringent designs."
        ),
        default=0,
    )
    aligner = luigi.Parameter(
        description=(
            "Alignment tool for off-target detection. "
            "'bowtie2' uses sensitive local alignment tuned for short probes "
            "(OligoMiner/Tigerfish parameter set). "
            "'bowtie' uses the legacy Bowtie aligner (max 2 mismatches). "
        ),
        default="bowtie2",
    )
    max_homopolymer_length = luigi.IntParameter(
        description=(
            "Maximum length of homopolymer run allowed in any base (A, T, C, or G). "
            "OligoMiner uses 5 (filters AAAAA, TTTTT, CCCCC, GGGGG). "
            "Set to 0 to disable. Replaces the legacy GGGG-only filter when > 0."
        ),
        default=5,
    )
    filter_low_complexity = luigi.BoolParameter(
        description="Filter probes containing low-complexity regions (dinucleotide repeats, low entropy).",
        default=False,
    )
    max_cpg_fraction = luigi.FloatParameter(
        description=(
            "Maximum CpG dinucleotide fraction allowed per probe (0.0-1.0). "
            "CpG-rich probes cross-hybridize with GC-rich rRNA causing nucleolar "
            "background (p=0.004 in 84-probe-set analysis). "
            "Set to 0 to disable. Recommended: 0.10 for stringent designs."
        ),
        default=0.0,
    )
    max_transcriptome_off_targets = luigi.IntParameter(
        description=(
            "Maximum number of transcriptome off-target hits per probe. "
            "A hit = >= blast_identity_threshold identity over >= min_blast_match_length. "
            "Only used when reference_transcriptome is provided."
        ),
        default=0,
    )
    blast_identity_threshold = luigi.FloatParameter(
        description=(
            "Minimum percent identity for a BLAST hit to count as an off-target. "
            "Only used when reference_transcriptome is provided."
        ),
        default=75.0,
    )
    off_target_min_tm = luigi.FloatParameter(
        description=(
            "Minimum predicted Tm for an off-target alignment to count as significant. "
            "Hits with Tm below this are ignored as thermodynamically unstable — "
            "rescues probes whose off-targets would not bind at hybridization conditions. "
            "Set to 0 to disable (count all alignment hits). "
            "Recommended: set to your hybridization temperature (e.g., 37)."
        ),
        default=0.0,
    )
    mask_repeats = luigi.BoolParameter(
        description=(
            "Ignore off-target hits in low-complexity/repetitive regions. "
            "Uses dustmasker (BLAST+) to identify repeats — "
            "rescues probes whose only off-targets are in repeats. "
            "Requires BLAST+ to be installed."
        ),
        default=False,
    )
    intergenic_off_targets = luigi.BoolParameter(
        description=(
            "Ignore off-target hits in intergenic regions (not overlapping any "
            "annotated gene). Requires --reference-annotation. "
            "Rescues probes whose off-targets fall outside annotated genes."
        ),
        default=False,
    )
    filter_rrna = luigi.BoolParameter(
        description=(
            "Remove probes with off-target hits on ribosomal RNA genes. "
            "rRNA is ~80%% of cellular RNA — even weak off-target binding "
            "causes intense background. Requires --reference-annotation. "
            "Catches 5S and mitochondrial rRNA from standard GTFs; "
            "for 18S/28S use --reference-transcriptome with rRNA sequences."
        ),
        default=False,
    )
    min_blast_match_length = luigi.IntParameter(
        description=(
            "Minimum effective alignment length (alignment_length - gaps) for a "
            "BLAST hit to count as an off-target in transcriptome filtering. "
            "Default is max(18, 0.8 * min_probe_length). Lower values increase "
            "sensitivity but may cause excessive filtering for exogenous genes."
        ),
        default=0,
    )
    max_probes_per_off_target = luigi.IntParameter(
        description=(
            "Maximum number of probes in the final set allowed to hit the same "
            "off-target gene. If exceeded, the weakest probes hitting that gene "
            "are removed. Set to 0 to disable. "
            "Prevents correlated off-target binding that creates false FISH spots."
        ),
        default=0,
    )
    adaptive_length = luigi.BoolParameter(
        description=(
            "Adjust probe length based on local GC content to normalize Tm. "
            "High-GC regions get shorter probes, low-GC regions get longer probes. "
            "Requires min_length < max_length with a range of at least 2nt."
        ),
        default=False,
    )
    filter_rdna_45s = luigi.BoolParameter(
        description=(
            "Screen probes against bundled 45S rDNA (U13369.1) and alpha satellite "
            "consensus to prevent nucleolar/centromeric off-target signal. "
            "The 45S rDNA is absent from standard genome assemblies (GRCh38) but "
            "present at ~300-400 copies — a single near-match probe creates "
            "overwhelming nucleolar background. Requires BLAST+."
        ),
        default=True,
    )
    custom_rdna_fasta = luigi.Parameter(
        description=(
            "Path to a custom rDNA FASTA file for non-human species. "
            "Overrides the bundled human 45S rDNA (U13369.1). "
            "Only used when --filter-rdna-45s is enabled."
        ),
        default="",
    )
    reject_cross_hybridization = luigi.BoolParameter(
        description=(
            "Hard-reject probes with strong cross-hybridization to off-target "
            "transcripts (>=16nt contiguous match at >=95%% identity, no gaps). "
            "Without this, such probes are only flagged as warnings."
        ),
        default=True,
    )
    allow_no_transcriptome = luigi.BoolParameter(
        description=(
            "Allow exogenous probe design without --reference-transcriptome. "
            "Without a transcriptome, exogenous off-target detection relies "
            "almost entirely on the rDNA/satellite filter and genome alignment, "
            "which may miss abundant expressed off-targets. "
            "Only use this if you accept reduced off-target protection."
        ),
        default=False,
    )
