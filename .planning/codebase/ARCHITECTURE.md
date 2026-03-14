# Architecture

**Analysis Date:** 2026-03-14

## Pattern Overview

**Overall:** Luigi task-based pipeline orchestration for RNA FISH probe design.

**Key Characteristics:**
- Task dependency graph managed by Luigi workflow engine
- Each processing stage is a discrete luigi.Task with defined inputs/outputs
- Configuration-driven via CLI-to-luigi parameter binding
- Lazy imports of heavy dependencies (bioinformatics tools)
- External tool invocation via subprocess (bowtie/bowtie2, jellyfish, BLAST, RNAstructure Fold)

## Layers

**Configuration & CLI Layer:**
- Purpose: Parse user arguments, validate parameters, manage run configuration
- Location: `eFISHent/cli.py`, `eFISHent/config.py`, `eFISHent/presets.py`, `eFISHent/constants.py`
- Contains: Argument parsing, type validators, parameter constraints, configuration classes
- Depends on: argparse, rich-argparse, luigi.Config
- Used by: main() entry point, all tasks

**Sequence Preparation Layer:**
- Purpose: Acquire target gene sequence (from file, NCBI Entrez, or Ensembl)
- Location: `eFISHent/prepare_sequence.py`
- Contains: DownloadEntrezGeneSequence, PrepareSequence tasks
- Depends on: Bio.SeqIO, subprocess (esearch/efetch), config
- Used by: GenerateAllProbes

**Candidate Generation Layer:**
- Purpose: Generate all possible probe subsequences with GC-adaptive length selection
- Location: `eFISHent/generate_probes.py`
- Contains: GenerateAllProbes task, create_candidate_probes_generator()
- Depends on: Bio.SeqRecord, configuration
- Used by: BasicFiltering

**Filtering & Validation Layers:**
- Purpose: Sequentially filter probes by biophysical and genomic criteria
- Location: Multiple files, each is a luigi.Task
  - `eFISHent/basic_filtering.py` - TM/GC content, homopolymer, complexity
  - `eFISHent/alignment.py` - Genome uniqueness via bowtie/bowtie2
  - `eFISHent/transcriptome_filter.py` - Off-target binding via BLAST (optional)
  - `eFISHent/kmers.py` - K-mer frequency via jellyfish index
  - `eFISHent/secondary_structure.py` - Folding energy via RNAstructure Fold
- Contains: Task classes for each filter type, helper functions for metrics
- Depends on: Bio utilities, external alignment/indexing tools, biopython
- Used by: Optimization layer

**Optimization Layer:**
- Purpose: Select maximum-coverage probe set from filtered candidates
- Location: `eFISHent/optimization.py`
- Contains: OptimizeProbeCoverage task with greedy and optimal (GLPK) algorithms
- Depends on: pandas, pyomo, numpy
- Used by: CleanUpOutput

**Output Finalization Layer:**
- Purpose: Aggregate results, compute statistics, format deliverables
- Location: `eFISHent/cleanup.py`
- Contains: CleanUpOutput task, prettify_table(), probe statistics computation
- Depends on: All filtering layer outputs, alignment data
- Used by: CLI/main()

**Index Building Layer (Parallel):**
- Purpose: Build reference indices used by filtering tasks
- Location: `eFISHent/indexing.py`, `eFISHent/kmers.py`
- Contains: BuildBowtieIndex, BuildBowtie2Index, BuildJellyfishIndex tasks
- Depends on: External tools (bowtie, bowtie2, jellyfish)
- Used by: Filtering layer

**Analysis Layer (Alternative Pipeline):**
- Purpose: Detailed statistical analysis of a provided probe set (--analyze-probeset mode)
- Location: `eFISHent/analyze.py`
- Contains: AnalyzeProbeset task with matplotlib visualization
- Depends on: All metrics functions, PDF generation
- Used by: CLI/main()

**Utility & Support Layer:**
- Purpose: Cross-cutting helper functions and logging
- Location: `eFISHent/util.py`, `eFISHent/console.py`, `eFISHent/gene_annotation.py`
- Contains: Path management, stage progress tracking, rich terminal UI, gene annotation lookup
- Used by: All tasks

## Data Flow

**Main Probe Design Pipeline:**

1. **Initialization**: User provides reference genome + target gene sequence
2. **Sequence Preparation**: DownloadEntrezGeneSequence → PrepareSequence (FASTA output)
3. **Candidate Generation**: GenerateAllProbes reads sequence, creates all subsequences
4. **Filtering Cascade**:
   - BasicFiltering: TM/GC/complexity → candidates_filtered.fasta
   - AlignProbeCandidates: Genome alignment via bowtie/bowtie2 → alignment.csv
   - TranscriptomeFiltering (optional): BLAST against transcriptome → txome_filtered.fasta
   - KMerFiltering: Jellyfish index lookup → kmers_filtered.fasta
   - SecondaryStructureFiltering: Fold predictions → sstruct_filtered.fasta
5. **Optimization**: OptimizeProbeCoverage selects final probe set via greedy/optimal algorithm
6. **Output**: CleanUpOutput aggregates all data into FASTA + CSV table + config file

**Index Building (Parallel):**
- BuildJellyfishIndex: Reference genome → jellyfish index (used by KMerFiltering)
- BuildBowtieIndex or BuildBowtie2Index: Reference genome → aligner index (used by AlignProbeCandidates)
- BuildTranscriptomeBlastDB: Transcriptome FASTA → BLAST database (used by TranscriptomeFiltering)

**Analysis Pipeline (Alternative):**
- AnalyzeProbeset reads user-provided FASTA file
- Computes metrics for all probes + gene alignment
- Generates multi-page PDF with histograms, boxplots, heatmaps

**State Management:**
- Config: Immutable throughout run, stored in luigi.cfg (created in temp directory)
- Intermediate files: Stored in output directory, optionally kept with --save-intermediates
- Task dependencies: Luigi automatically reruns tasks if inputs are missing/stale
- Progress tracking: Captured via PIPELINE_STAGES dict in util.py, displayed via console.py

## Key Abstractions

**Luigi Task:**
- Purpose: Represents a single processing stage with requires() → run() → output() contract
- Examples: `eFISHent/generate_probes.py` GenerateAllProbes, `eFISHent/optimization.py` OptimizeProbeCoverage
- Pattern: Each task depends on previous tasks via requires(), processes, writes to output()

**SeqRecord:**
- Purpose: Biopython object representing a DNA sequence with metadata
- Examples: Probe candidates, filtered probes, final selected probes
- Pattern: Created by GenerateAllProbes, modified by each filtering task, output to FASTA files

**Configuration Classes:**
- Purpose: Group related parameters with validation (GeneralConfig, RunConfig, SequenceConfig, ProbeConfig)
- Examples: `eFISHent/config.py` luigi.Config subclasses
- Pattern: Instantiated fresh per task, accessed via config().param_name

**Validator Functions:**
- Purpose: Type coercion + validation for CLI arguments
- Examples: `eFISHent/cli.py` string_to_bool(), positive_int(), existing_file()
- Pattern: Applied to each argument via PARAM_VALIDATORS dict during parsing

## Entry Points

**Main CLI Entry Point:**
- Location: `eFISHent/__main__.py` → `eFISHent/cli.py` main()
- Triggers: User runs `efishent` or `eFISHent` command
- Responsibilities:
  1. Parse CLI arguments via _parse_args()
  2. Validate cross-parameter constraints
  3. Handle special modes: --check (dependency verification), --update (self-update), --preset (configuration templates)
  4. Create temp luigi.cfg with parsed parameters
  5. Instantiate root task (CleanUpOutput, AnalyzeProbeset, or index-building tasks)
  6. Build task dependency graph via luigi.build()
  7. Display progress and completion summary

**Index Building Entry Point:**
- Triggered: --build-indices flag
- Tasks: BuildJellyfishIndex, BuildBowtieIndex/BuildBowtie2Index
- Output: Indices stored in output directory for reuse

**Probeset Analysis Entry Point:**
- Triggered: --analyze-probeset FILE flag
- Task: AnalyzeProbeset
- Output: PDF file with statistical analysis

## Error Handling

**Strategy:** Fail-fast validation with user-friendly messages

**Patterns:**
- **CLI Validation**: Cross-parameter constraints checked in validate_args() before task execution
- **File Validation**: existing_file(), existing_fasta_file(), existing_gtf_file() type validators prevent invalid paths
- **Dependency Checking**: check_required_dependencies() verifies external tools before pipeline start
- **Task Failures**: Luigi catches subprocess errors; logs written to efishent.log in debug mode
- **Parameter Warnings**: validate_parameter_warnings() identifies risky combinations (narrow TM window, high formamide, etc.) and surfaces as non-fatal warnings
- **Logging Levels**: Debug (file output), info (console progress), warning (non-fatal issues), error (fatal)

## Cross-Cutting Concerns

**Logging:**
- Framework: Python logging module
- Configuration: Centralized via set_logging_level() in cli.py
- Patterns: Custom logger "custom-logger" used by all tasks, Luigi loggers suppressed to WARNING
- Files: Debug mode writes to efishent.log, silent mode uses standard handlers

**Validation:**
- CLI arguments validated by type validators + cross-parameter checks in validate_args()
- External tool availability checked via check_required_dependencies()
- Parameter constraints validated post-parsing (min ≤ max for TM/GC/length ranges)

**Authentication:**
- Not applicable (no external services requiring auth)

**Multiprocessing:**
- Used in: Secondary structure prediction (get_free_energy), cleanup statistics (prettify_table)
- Pattern: multiprocessing.Pool(GeneralConfig().threads)
- Configuration: Thread count set via --threads parameter

**Progress Tracking:**
- Framework: Rich terminal UI via rich library
- Patterns: Pipeline stages defined in PIPELINE_STAGES dict (util.py), displayed via pipeline_progress() context manager (console.py)
- Modes: Full pipeline (8 stages), index building (2 stages), analysis (9 steps)

