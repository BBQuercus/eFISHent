# Codebase Structure

**Analysis Date:** 2026-03-14

## Directory Layout

```
efishent/
├── eFISHent/                       # Main package (20 Python modules)
│   ├── __init__.py                 # Version export
│   ├── __main__.py                 # Entry point for python -m eFISHent
│   ├── cli.py                      # Command-line interface, argument parsing, dependency checks
│   ├── config.py                   # Luigi configuration classes (4 config groups)
│   ├── constants.py                # Global constants (SAM flags, CLI shortforms, column definitions)
│   ├── presets.py                  # Parameter presets (smfish, merfish, dna-fish, strict, relaxed)
│   ├── console.py                  # Rich terminal UI, progress bars, formatted output
│   ├── util.py                     # Utility functions (paths, logging, stage tracking, data tables)
│   ├── generate_probes.py          # GenerateAllProbes task (candidate generation)
│   ├── basic_filtering.py          # BasicFiltering task (TM/GC/complexity filters)
│   ├── alignment.py                # AlignProbeCandidates task (bowtie/bowtie2 alignment)
│   ├── kmers.py                    # BuildJellyfishIndex, KMerFiltering tasks
│   ├── secondary_structure.py      # SecondaryStructureFiltering task (RNAstructure Fold)
│   ├── optimization.py             # OptimizeProbeCoverage task (greedy/optimal algorithms)
│   ├── cleanup.py                  # CleanUpOutput task (final output formatting)
│   ├── prepare_sequence.py         # PrepareSequence, DownloadEntrezGeneSequence tasks
│   ├── indexing.py                 # BuildBowtieIndex, BuildBowtie2Index tasks
│   ├── transcriptome_filter.py     # TranscriptomeFiltering, BuildTranscriptomeBlastDB tasks
│   ├── analyze.py                  # AnalyzeProbeset task (probe set analysis + PDF)
│   ├── gene_annotation.py          # Gene annotation lookup (ENSEMBL/NCBI mapping)
│   ├── data_tables/                # RNAstructure thermodynamic tables (included)
│   ├── Fold_osx                    # RNAstructure Fold binary (macOS)
│   ├── Fold_linux                  # RNAstructure Fold binary (Linux)
│   └── luigi.cfg                   # Luigi configuration template
│
├── tests/                          # Test suite
│   ├── __init__.py
│   ├── conftest.py                 # pytest configuration (currently commented)
│   ├── data/                       # Test fixtures (reference genomes, sequences)
│   ├── test_unit_*.py              # Unit tests (11 files)
│   ├── test_integration.py         # Integration tests
│   ├── test_adaptive_probes.py    # Tests for adaptive probe length
│   ├── test_gene_annotation.py     # Tests for gene annotation mapping
│   └── test_quality_scores.py      # Tests for probe quality metrics
│
├── scripts/                        # Utility scripts
│
├── pyproject.toml                  # Project metadata, build config, tool settings
├── setup.py                        # Setup script (legacy, setuptools configured via pyproject.toml)
├── setup.cfg                       # setuptools configuration
├── README.md                       # User documentation
│
├── .planning/                      # GSD planning documents
│   └── codebase/                   # Codebase analysis documents
│       ├── ARCHITECTURE.md         # (This file group)
│       └── STRUCTURE.md
│
├── .github/                        # GitHub configuration
├── refs/                           # Reference files (unused in active codebase)
├── TrueProbes-1/                   # Legacy TrueProbes project (research reference)
└── test_output/                    # Test execution output
```

## Directory Purposes

**eFISHent/ (Main Package):**
- Purpose: Core probe design library and CLI
- Contains: 20 Python modules organized by function (tasks, utilities, UI, configuration)
- Key files: cli.py (entry point), config.py (parameters), 14 luigi Task classes

**tests/:**
- Purpose: Test suite covering units, integration, and probe quality
- Contains: 15 test files, fixture data in tests/data/
- Key files: test_unit_*.py (individual module tests), test_integration.py (full pipeline), test_quality_scores.py (statistical validation)

**scripts/:**
- Purpose: Utility scripts for development/analysis
- Contains: Supporting tools (currently empty/minimal)

**data_tables/ (in eFISHent/):**
- Purpose: RNAstructure thermodynamic reference tables
- Contains: Precomputed folding coefficients required by Fold binary
- Generated: No, distributed with package
- Committed: Yes

## Key File Locations

**Entry Points:**
- `eFISHent/__main__.py`: Module entry point (python -m eFISHent)
- `eFISHent/cli.py`: main() function, CLI orchestration
- `pyproject.toml` [project.scripts]: `efishent` and `eFISHent` console commands map to cli:main()

**Configuration:**
- `eFISHent/cli.py`: Argument parsing, parameter validation, Luigi config creation
- `eFISHent/config.py`: 4 luigi.Config subclasses (GeneralConfig, RunConfig, SequenceConfig, ProbeConfig)
- `eFISHent/luigi.cfg`: Template configuration file (overridden at runtime)
- `eFISHent/presets.py`: Named parameter sets (smfish, merfish, dna-fish, strict, relaxed)
- `eFISHent/constants.py`: Global constants (CLI shortforms, column definitions, flags)

**Core Logic:**
- **Sequence Preparation**: `eFISHent/prepare_sequence.py` (DownloadEntrezGeneSequence, PrepareSequence)
- **Candidate Generation**: `eFISHent/generate_probes.py` (GenerateAllProbes, create_candidate_probes_generator)
- **Basic Filtering**: `eFISHent/basic_filtering.py` (BasicFiltering, TM/GC/complexity metrics)
- **Genome Alignment**: `eFISHent/alignment.py` (AlignProbeCandidates)
- **K-mer Filtering**: `eFISHent/kmers.py` (KMerFiltering, BuildJellyfishIndex)
- **Secondary Structure**: `eFISHent/secondary_structure.py` (SecondaryStructureFiltering, get_free_energy)
- **Transcriptome Filtering**: `eFISHent/transcriptome_filter.py` (TranscriptomeFiltering, BuildTranscriptomeBlastDB)
- **Optimization**: `eFISHent/optimization.py` (OptimizeProbeCoverage, greedy_model, optimal_model)
- **Output Finalization**: `eFISHent/cleanup.py` (CleanUpOutput, prettify_table)
- **Analysis**: `eFISHent/analyze.py` (AnalyzeProbeset)
- **Indexing**: `eFISHent/indexing.py` (BuildBowtieIndex, BuildBowtie2Index)

**Testing:**
- `tests/conftest.py`: pytest configuration
- `tests/data/`: Test fixtures (reference genome, sequences, annotation files)
- `tests/test_unit_cli.py`: CLI parsing and validation tests
- `tests/test_unit_generate_probes.py`: Candidate generation tests
- `tests/test_unit_basic_filtering.py`: Filtering metric tests
- `tests/test_unit_alignment.py`: Alignment task tests
- `tests/test_quality_scores.py`: Probe quality scoring validation
- `tests/test_gene_annotation.py`: Gene annotation mapping tests

**Utilities & Support:**
- `eFISHent/util.py`: Path/file management, stage tracking (PIPELINE_STAGES), data table creation
- `eFISHent/console.py`: Rich terminal UI (progress bars, tables, formatted output)
- `eFISHent/gene_annotation.py`: Gene name/ID mapping via ENSEMBL/NCBI

## Naming Conventions

**Files:**
- Task modules: Lowercase with underscores (e.g., `generate_probes.py`, `basic_filtering.py`)
- Utility modules: Lowercase (e.g., `util.py`, `console.py`, `constants.py`)
- Binary executables: CamelCase with platform suffix (e.g., `Fold_osx`, `Fold_linux`)

**Directories:**
- Package directory: CamelCase (`eFISHent/`)
- Test directory: Lowercase (`tests/`)
- Support directories: Lowercase (`scripts/`, `refs/`)

**Functions:**
- Public functions: snake_case (e.g., `create_candidate_probes_generator()`, `get_free_energy()`)
- Private functions: Prefix with underscore (e.g., `_add_groups()`, `_parse_args()`, `_ensure_deps_on_path()`)
- Type validators: snake_case describing validation (e.g., `existing_file()`, `positive_int()`, `percentage()`)

**Classes:**
- Task classes: PascalCase, descriptive (e.g., `GenerateAllProbes`, `BasicFiltering`, `OptimizeProbeCoverage`)
- Configuration classes: PascalCase with "Config" suffix (e.g., `GeneralConfig`, `ProbeConfig`)

**Constants:**
- Global constants: UPPER_SNAKE_CASE (e.g., `SAMFILE_COLUMNS`, `BLAST_COLUMNS`, `PIPELINE_STAGES`)
- Private constants: Prefix with underscore (e.g., `_BOOL_METAVAR`)

**Variables:**
- Local variables: snake_case (e.g., `probe_sequences`, `alignment_path`, `filtered_candidates`)
- Configuration keys: snake_case matching luigi.Parameter names (e.g., `min_length`, `max_tm`, `reference_genome`)

## Where to Add New Code

**New Feature (Filter Type):**
- Primary code: Create new file `eFISHent/feature_name.py` with a luigi.Task class
  - Follow pattern: inherit from luigi.Task, implement requires() → run() → output()
  - Use `util.log_stage_start()` for progress tracking
  - Read input FASTA via `Bio.SeqIO.parse()`, write output via `Bio.SeqIO.write()`
  - Store metrics in task attributes for cleanup.py aggregation
- Tests: Create `tests/test_unit_feature_name.py` with pytest tests
- Integration: Add task to cleanup.py requires() dict, update PIPELINE_STAGES in util.py

**New Component/Module:**
- Utility function: Add to `eFISHent/util.py` if generic, otherwise to relevant module
- Task class: Follow new feature pattern above
- Configuration: Add new luigi.Parameter to appropriate Config class in `eFISHent/config.py`
- CLI argument: Add to CONFIG_CLASSES list in config.py, parameter automatically surfaces in CLI

**Utilities/Helpers:**
- Shared helpers: `eFISHent/util.py` for path management, logging, generic operations
- Filtering metrics: `eFISHent/basic_filtering.py` for reusable scoring functions (get_melting_temp, get_gc_content, etc.)
- Console formatting: `eFISHent/console.py` for Rich terminal output
- Annotation lookups: `eFISHent/gene_annotation.py` for ENSEMBL/NCBI mapping

**Validators/Parsers:**
- CLI validators: Add to PARAM_VALIDATORS dict in `eFISHent/cli.py`
- Type coercion functions: Add below parameter validators in cli.py (e.g., string_to_bool, positive_int)

## Special Directories

**data_tables/ (eFISHent/data_tables/):**
- Purpose: RNAstructure thermodynamic reference coefficients
- Generated: No (distributed with package)
- Committed: Yes
- Usage: Set as DATAPATH environment variable when invoking Fold binary

**Fold_osx, Fold_linux (eFISHent/):**
- Purpose: Pre-compiled RNAstructure Fold binary for secondary structure prediction
- Generated: No (compiled separately, included in distribution)
- Committed: Yes (essential dependency)
- Usage: Invoked directly by secondary_structure.py get_free_energy() function

**tests/data/ (tests/data/):**
- Purpose: Test fixtures (reference genomes, GTF annotations, sequence files)
- Generated: No (checked in)
- Committed: Yes
- Usage: Loaded by test files to avoid external dependencies

**test_output/ (test_output/):**
- Purpose: Temporary directory for pipeline test outputs
- Generated: Yes (created during testing)
- Committed: No (git-ignored)
- Cleanup: Can be safely deleted

**TrueProbes-1/ (TrueProbes-1/):**
- Purpose: Legacy predecessor project reference
- Generated: No (historical archive)
- Committed: Yes (for context)
- Usage: None in active codebase

