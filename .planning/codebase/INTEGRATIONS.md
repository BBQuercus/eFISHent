# External Integrations

**Analysis Date:** 2026-03-14

## APIs & External Services

**NCBI (Entrez):**
- Service: NCBI Entrez Direct (eutils) - Gene sequence download
- What it's used for: Automatic gene sequence retrieval from NCBI when user provides gene name + organism or Ensembl ID
- SDK/Client: Entrez Direct command-line tools (`esearch`, `efetch`)
  - Implementation: `eFISHent/prepare_sequence.py:DownloadEntrezGeneSequence`
  - Query databases: `gene` (Ensembl ID → NCBI Gene ID resolution), `nuccore` (transcript sequence retrieval)
  - Filters: RefSeq only, mRNA/ncRNA biomol types
- Auth: None required (public API, rate-limited)

## Data Storage

**Databases:**
- None - Application uses user-provided reference files
- Types of inputs accepted:
  - Reference genome: FASTA format (gzip allowed)
  - Genome annotation: GTF/GFF format
  - Transcriptome: FASTA format (generated from genome + GTF)
  - Expression count table: CSV/TSV format (normalized RNA-seq counts)

**File Storage:**
- Local filesystem only
- Default output: Current working directory
- Configurable via `--output-dir` parameter
- Output files:
  - `{GENE}_{HASH}.fasta` - Probe sequences
  - `{GENE}_{HASH}.csv` - Detailed probe metrics (quality scores, off-targets, Tm, etc.)
  - `{GENE}_{HASH}.txt` - Configuration/parameters used
  - `{GENE}_{HASH}_analysis.pdf` - Optional analysis report

**Intermediate Files:**
- Stored in temp/output directory during workflow
- Can be preserved with `--save-intermediates True`
- Deleted automatically at end unless flag set

**Caching:**
- Index files (built once per genome):
  - Bowtie2 index: `.bt2` files in genome directory
  - Bowtie index: `.ebwt` files in genome directory
  - Jellyfish index: `.jf` file in output directory
  - BLAST database: genome.nhr/nin/nsq files (if `--reference-transcriptome` used)

## Authentication & Identity

**Auth Provider:**
- None - Tool is standalone CLI, no user authentication required
- NCBI Entrez API access: Public (no API key required but subject to rate limiting)

**Implementation:**
- No login/session management
- All parameters passed via command-line arguments
- No secrets management (environment variables or config files)

## Monitoring & Observability

**Error Tracking:**
- None - No external error tracking services integrated
- Local logging to stderr with Rich formatting

**Logs:**
- Rich-based terminal output (colors, progress bars, formatted tables)
- Logger: Python logging module with custom handlers
- Levels: DEBUG, INFO, WARNING, ERROR
- Can silence non-critical output with `--silent` flag
- No persistent log files by default (printed to console)

## CI/CD & Deployment

**Hosting:**
- GitHub repository: https://github.com/BBQuercus/eFISHent
- Distribution: PyPI (as Python package) + conda-forge (bioconda)

**CI Pipeline:**
- GitHub Actions: `.github/workflows/` (present)
- Workflow: Tests on Python 3.9+ across platforms
- Badges in README:
  - Github Actions Status (Tests workflow)
  - CodeFactor code quality score
  - codecov coverage report (token: C1SRFYZ5VP)

**Update Mechanism:**
- `efishent --update` command triggers pip upgrade
- Verifies external dependencies still present after Python update

## Environment Configuration

**Required env vars:**
- None explicit - Tool manages PATH internally
- Modified at runtime by CLI:
  - `PATH` - Prepended with `~/.local/efishent/deps/bin` and edirect directory
  - `LD_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` - Set for BLAST library paths
  - `DATAPATH` - Set to thermodynamic data tables directory during secondary structure prediction

**Secrets location:**
- No secrets used - Application is open-source bioinformatics tool
- NCBI Entrez access is public (unauthenticated)

**Install Configuration:**
- Installer: `install.sh`
- Default prefix: `~/.local/efishent`
- Creates wrapper at: `~/.local/bin/efishent`
- Modified shell RC files: `.bashrc`, `.zshrc` (unless `--no-modify-rc` flag)
- Can specify custom prefix: `--prefix /path/to/install`

## Webhooks & Callbacks

**Incoming:**
- None - CLI tool, not a server

**Outgoing:**
- None - No callbacks or webhooks to external services

## External Tool Integration

**Command-line Tools Executed (via subprocess):**

**Alignment:**
- `bowtie2` - Sensitive local alignment (default aligner)
  - Command: `bowtie2 -x <index> -f -U <probes.fasta> [OligoMiner params]`
  - Location: `eFISHent/alignment.py:AlignProbeCandidates`
- `bowtie` - Legacy aligner (when `--aligner bowtie` specified)
  - Command: `bowtie <index> <probes.fasta>`
  - Parameters: 2-mismatch tolerance

**K-mer Filtering:**
- `jellyfish count` - Build k-mer index from reference genome
- `jellyfish query` - Count k-mer occurrences in probes
  - Location: `eFISHent/kmers.py`

**Secondary Structure Prediction:**
- `Fold` - RNAstructure program for ΔG prediction
  - Binary: `eFISHent/Fold_linux` or `Fold_osx` (platform-specific)
  - Output: Free energy values for secondary structure stability
  - Location: `eFISHent/secondary_structure.py:PredictSecondaryStructure`
  - Environment: `DATAPATH` set to thermodynamic tables directory

**Transcriptome BLAST:**
- `makeblastdb` - Build nucleotide BLAST database from transcriptome FASTA
- `blastn` - Alignment with TrueProbes parameters (Neuert Lab, 2025)
  - Command: `blastn -db <transcriptome.db> -query <probes.fasta> -outfmt 6 [parameters]`
  - Parameters tuned for sensitive detection of off-target binding
  - Location: `eFISHent/transcriptome_filter.py`
- `dustmasker` - Low-complexity region masking (from BLAST+)

**GTF/Annotation Processing:**
- `gffread` - Convert GTF + genome to transcriptome FASTA (user-provided, optional)
  - Not called by eFISHent directly; users run: `gffread annotation.gtf -g genome.fa -w transcriptome.fa`

**NCBI Data Download:**
- `esearch` - Search NCBI databases (gene, nuccore)
  - Example: `esearch -db gene -query "ACTB AND homo sapiens[Organism]"`
- `efetch` - Fetch sequences in FASTA format
  - Example: `efetch -db nuccore -id <UID> -format fasta`
  - Location: `eFISHent/prepare_sequence.py:DownloadEntrezGeneSequence`

## Optimization Solver Integration

**GLPK (GNU Linear Programming Kit):**
- When: Only used if `--optimization-method optimal` is specified
- Implementation: Pyomo wrapper around GLPK
- Model: Integer linear program (MILP) to maximize nucleotide coverage
- Location: `eFISHent/optimization.py:OptimizeProbes`

## Data Sources

**Reference Data:**
- User provides:
  - Reference genome (FASTA, e.g., GRCh38, dm6)
  - Gene annotation (GTF, e.g., Ensembl, GENCODE)
  - Optional: Reference transcriptome FASTA
  - Optional: RNA-seq count table for expression weighting
- Can download from:
  - Ensembl FTP: https://ftp.ensembl.org/pub/
  - UCSC Genome Browser: https://hgdownload.soe.ucsc.edu/downloads.html
  - GEO: https://www.ncbi.nlm.nih.gov/gds/ (for expression data)
  - Expression Atlas: https://www.ebi.ac.uk/gxa/home

**Embedded Data:**
- Thermodynamic parameters: `eFISHent/data_tables/*.dat` (~50 files)
  - Used by RNAstructure for nearest-neighbor ΔG calculation
  - Includes DNA/RNA stacking, loop, coaxial, and dangle parameters
  - Source: RNAstructure package (Mathews lab)

## Network Dependencies

**During Installation:**
- Downloads pre-compiled binaries from GitHub releases and mirrored sources
- Downloads Python packages from PyPI (via pip/uv)
- One-time network access; works offline after installation

**During Execution:**
- NCBI Entrez queries (if `--gene-name` or `--ensembl-id` provided)
- No other network calls during probe design workflow

---

*Integration audit: 2026-03-14*
