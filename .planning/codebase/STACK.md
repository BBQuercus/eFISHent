# Technology Stack

**Analysis Date:** 2026-03-14

## Languages

**Primary:**
- Python 3.9+ - Core application language for probe design logic, sequence processing, and analysis

**Secondary:**
- Bash - Installation and environment setup (`install.sh`)
- C/C++ - External tools (Bowtie2, Jellyfish, BLAST+, RNAstructure) distributed as pre-compiled binaries

## Runtime

**Environment:**
- Python 3.9 or later required
- Runs on macOS, Linux, Windows (via WSL), and shared HPC/cluster servers

**Package Manager:**
- uv (UV package manager)
- Setuptools for distribution
- Lockfile: `uv.lock` (present)

## Frameworks

**Core:**
- Luigi 3.x - Workflow orchestration and task dependency management for multi-stage probe design pipeline
- Biopython - Sequence processing, SeqIO, alignment analysis, melting temperature calculations

**Analysis/Data:**
- Pandas - Data table processing and CSV output generation
- NumPy - Numerical computations for optimization and scoring
- Matplotlib - Visualization for probe set analysis reports and distributions
- PyArrow (>=18.0.0) - Parquet/Arrow data format support

**Optimization:**
- Pyomo - Mathematical optimization modeling for optimal probe coverage (MILP solver interface)
- GLPK (5.0) - Linear programming solver for optimization (installed via external dependency)

**UI/Logging:**
- Rich (>=13.0.0) - Enhanced terminal output, progress bars, and formatted logging
- Rich-argparse (>=1.4.0) - Rich formatting for command-line help and argument documentation

**Bioinformatics/GTF:**
- GTFparse - Parsing and handling of genome annotation (GTF/GFF) files
- Pysam - SAM/BAM file handling and sequence alignment processing

## Key Dependencies

**Critical:**
- biopython - DNA/RNA sequence manipulation, melting temperature calculations using nearest-neighbor thermodynamics
- luigi - Manages workflow dependencies between genome indexing, alignment, filtering stages
- pandas - Outputs probe quality reports (CSV format) with per-probe metrics
- pyomo - Mathematical optimization for finding maximum coverage probe sets with constraints

**Infrastructure:**
- gtfparse - Parses genome annotations to identify genes, transcripts, biotypes (rRNA detection)
- pysam - Processes alignment outputs from Bowtie2/Bowtie and transcriptome BLAST results
- pyarrow>=18.0.0 - Efficient data serialization for intermediate results
- matplotlib - Generates PDF analysis reports with probe length/Tm/GC distributions
- numpy - Numerical operations in optimization and sequence scoring
- rich>=13.0.0 - User-facing terminal interface with progress visualization

## External Tools (Pre-compiled Binaries)

**Genome Alignment:**
- Bowtie2 (default) - Sensitive local alignment with OligoMiner/Tigerfish parameters for short probes
- Bowtie (legacy) - Alternative aligner with 2-mismatch tolerance

**K-mer Filtering:**
- Jellyfish 2.3.1 - Fast k-mer frequency counting in reference genome

**Secondary Structure Prediction:**
- RNAstructure Fold (6.4) - Nearest-neighbor thermodynamic model for RNA secondary structure (ΔG prediction)
- Fold_osx, Fold_linux - Platform-specific binaries included in package

**Transcriptome Filtering:**
- BLAST+ 2.17.0 (optional) - blastn for transcriptome BLAST off-target detection
- dustmasker (from BLAST+) - Low-complexity region identification for repeat masking

**Genome Processing:**
- gffread (optional) - Converts GTF + genome FASTA to transcriptome FASTA format

**NCBI Integration:**
- Entrez Direct (eutils) - esearch, efetch for downloading gene sequences from NCBI
- Installed via install.sh (EDirect 2024 release)

**Optimization Solver:**
- GLPK 5.0 - Linear programming solver backend for pyomo (optimal coverage method)

## Configuration

**Build System:**
- `pyproject.toml` - Modern Python packaging configuration (setuptools build)
- `setup.cfg` - Flake8 and Pydocstyle linting rules (max line 120 characters)
- `setup.py` - Legacy compatibility setup script

**CLI/Entry Points:**
- Main: `efishent` → `eFISHent.cli:main`
- Alias: `eFISHent` → `eFISHent.cli:main`

**Package Data:**
- Pre-compiled binaries: `Fold_osx`, `Fold_linux` (included in package)
- Thermodynamic data tables: `data_tables/*` (RNAstructure nearest-neighbor parameters, ~50 .dat files)

**Installation:**
- Script: `install.sh` - One-command installation that:
  - Downloads pre-compiled binaries (Bowtie2, Jellyfish, BLAST+, GLPK, Entrez Direct)
  - Creates isolated venv at `~/.local/efishent/venv`
  - Installs Python package via pip
  - Creates wrapper command in `~/.local/bin/efishent`
  - No sudo, Docker, or conda required

**Dependencies Directory Structure (after install):**
```
~/.local/efishent/
├── deps/
│   ├── bin/           # Bowtie2, Bowtie, Jellyfish, BLAST+, gffread, GLPK
│   ├── lib/           # Shared libraries
│   └── edirect/       # Entrez Direct (esearch, efetch)
├── venv/              # Python virtual environment
└── (Python package installed in venv/lib/pythonX.Y/site-packages/)
```

## Platform Requirements

**Development:**
- Python 3.9+
- Git
- curl or wget for downloads
- C/C++ compiler (for building some external tools)
- Supports macOS (Intel/Apple Silicon), Linux (x86_64/ARM), Windows (WSL2)

**Production:**
- Python 3.9+
- ~500MB disk space for binary dependencies and venv
- Multi-core CPU recommended (supports `--threads` parameter)
- Works on shared HPC/cluster servers via SSH (no root required)

**Type Checking:**
- MyPy (strict mode configured in pyproject.toml)
- Check untyped defs, disallow generics, strict equality

---

*Stack analysis: 2026-03-14*
