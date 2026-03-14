# Codebase Concerns

**Analysis Date:** 2026-03-14

## Tech Debt

**Subprocess Management Without Consistency:**
- Issue: The codebase uses `subprocess.check_call()` for bowtie/bowtie2 execution, but `subprocess.run()` with `check=True` for external tools. No timeout specified for long-running alignment jobs.
- Files: `eFISHent/alignment.py` (lines 167, 210), `eFISHent/indexing.py` (lines 41, 73), `eFISHent/kmers.py` (line 174), `eFISHent/transcriptome_filter.py` (lines 45, 113)
- Impact: Alignment jobs can hang indefinitely. If bowtie/bowtie2 process stalls, the entire pipeline blocks without recovery.
- Fix approach: Add `timeout` parameter to all `subprocess` calls. Consider creating a helper function `run_external_tool()` that enforces timeout, error handling, and logging consistently across all external commands.

**Duplicate Binding Affinity Logic:**
- Issue: The `is_binding()` function in `eFISHent/optimization.py` is duplicated with custom implementation in `eFISHent/analyze.py` (lines 152-161 has TODO comment noting this).
- Files: `eFISHent/optimization.py` (lines 250-268), `eFISHent/analyze.py` (lines 153-161)
- Impact: Maintenance burden. Bug fixes must be applied to two locations. Inconsistent behavior if one is updated.
- Fix approach: Extract `is_binding()` to shared utility module (`eFISHent/util.py` or new `sequence_utils.py`) and import from both locations.

**Undocumented Mathematical Model:**
- Issue: The optimal model in `eFISHent/optimization.py` (line 292) has TODO comment: "add proper mathematical description"
- Files: `eFISHent/optimization.py` (lines 293-345)
- Impact: Complex optimization logic is hard to verify or maintain. New developers cannot understand the constraints and objective function without reverse-engineering the Pyomo model definition.
- Fix approach: Add docstring with mathematical notation explaining the coverage problem, constraint definitions (C1: probe assignment, C2: no over-coverage), and objective.

## Known Bugs

**Soft Exception Handling in File Reading:**
- Issue: Multiple functions return `None` or empty lists without proper error propagation when files fail to read or parse.
- Files: `eFISHent/cleanup.py` (lines 145-146, 152, 155, 172, 204, 209, 278, 494, 498, 577), `eFISHent/kmers.py` (line 83)
- Trigger: When BLAST CSV file is corrupted, transcriptome FASTA is unreadable, or temporary files are removed prematurely
- Workaround: Pipeline continues silently; final output contains incomplete data (e.g., missing `txome_off_targets` column)
- Fix approach: Use explicit logging before early returns. Raise exceptions instead of returning None for critical paths. Allow callers to decide if missing data is recoverable.

**Bowtie2 Output Parsing Fragility:**
- Issue: FASTA description field characters (angle brackets, pipes) cause bowtie2 to corrupt read names in SAM output. Workaround in `eFISHent/alignment.py` (lines 181-184) strips BioPython formatting by rewriting FASTA, but original sequence metadata is lost.
- Files: `eFISHent/alignment.py` (lines 178-184)
- Trigger: When BioPython creates sequence records with default "<unknown description>" template
- Impact: Probe metadata (e.g., original probe position) cannot be recovered from SAM output
- Fix approach: Use `Bio.SeqIO.FastaIO` with custom header formatting function. Validate probe IDs are unique and contain no special characters before alignment.

**Unsafe Hardcoded Platform Detection:**
- Issue: Secondary structure prediction selects binary based on `sys.platform` check (`eFISHent/secondary_structure.py`, lines 32-40). No validation that binary exists or is executable.
- Files: `eFISHent/secondary_structure.py` (lines 32-40)
- Trigger: If `Fold_linux` binary is missing or permissions wrong, subprocess.run fails at runtime with unhelpful error
- Workaround: Manual PATH setup; silent failures if binary missing
- Fix approach: Call a validation function at startup (already exists: `check_all_dependencies()` in cli.py). Store Fold executable path in config and validate presence before first use.

**Off-Target Melting Temp Rescue Logic Catches KeyError But Continues:**
- Issue: In `eFISHent/alignment.py` (lines 527-529), when genome.fetch() raises ValueError or KeyError (e.g., chromosome not found), the code increments `significant_count` and continues, treating missing reference as equivalent to a hit.
- Files: `eFISHent/alignment.py` (lines 527-529)
- Impact: Probes against non-existent genome regions are incorrectly "rescued" because missing fetch is counted as a significant off-target match
- Fix approach: Log why fetch failed (missing chrom? out of bounds?). Only increment count if it's a real alignment issue, not a reference lookup failure.

## Security Considerations

**No Shell Injection Protection for NCBI Queries:**
- Risk: User provides `ensembl_id` or `organism_name` that are embedded into NCBI `esearch` command via list (safe) but the query string building is not validated.
- Files: `eFISHent/prepare_sequence.py` (lines 40-45, 86-90)
- Current mitigation: subprocess uses list args (not shell=True), so command injection is prevented at OS level. No special characters are passed.
- Recommendations: Add validation in CLI to reject ensembl_id/organism_name with unexpected characters (non-alphanumeric, hyphens, underscores only). Document this limit.

**BLAST Output Parsing Could Be Exploited:**
- Risk: BLAST results are parsed as CSV, then gene names extracted via regex. Malformed BLAST output or crafted input could cause parsing errors or unexpected data in gene lists.
- Files: `eFISHent/cleanup.py` (lines 149-151), `eFISHent/gene_annotation.py`
- Current mitigation: DataFrame operations ignore malformed lines; regex is bounded.
- Recommendations: Validate BLAST output format before parsing. Add schema validation for result dataframes.

**Temporary Files Not Explicitly Cleaned:**
- Risk: Temp files created in jellyfish k-mer counting (tempfile.NamedTemporaryFile in `eFISHent/kmers.py`) and BLAST operations are cleaned with `os.unlink()`, but if exception occurs, cleanup may be skipped.
- Files: `eFISHent/kmers.py` (lines 46, 57), `eFISHent/cleanup.py` (line 271-272)
- Impact: Disk space leaks if pipeline is interrupted or exceptions occur
- Recommendations: Use context managers (`TemporaryDirectory()`) or try/finally blocks explicitly.

## Performance Bottlenecks

**O(n²) Binding Affinity Matrix in Analysis:**
- Problem: When analyzing probes, `eFISHent/analyze.py` (lines 171-176) builds a full sequence similarity matrix comparing all probes against all other probes using `Bio.pairwise2.align.globalxx()`.
- Files: `eFISHent/analyze.py` (lines 163-178)
- Cause: For 100 probes, this is 10,000 alignments; for 500 probes, it's 250,000 alignments. Each alignment is 20ms+.
- Improvement path: Warn users about runtime (already at line 167-170). For large probe sets (>50 probes), use faster heuristic (k-mer overlap) instead of exact alignment. Implement incremental updates if probes are sorted.

**Redundant Transcriptome BLAST Rescoring in Cleanup:**
- Problem: `eFISHent/cleanup.py` (lines 143-146) re-BLASTs final probes against transcriptome, even though transcriptome filtering already ran during main pipeline.
- Files: `eFISHent/cleanup.py` (lines 143-146, 199-278)
- Cause: Final probes may not be in intermediate transcriptome hit CSV; re-BLAST ensures accurate final counts.
- Improvement path: Cache/reuse transcriptome hits from filtering stage if probes haven't changed. Alternatively, document why re-run is necessary and consider lazy evaluation (only compute on request).

**Melting Temp Uniformity Refinement O(n²) Loop:**
- Problem: `eFISHent/optimization.py` (lines 139-242) iterates through assigned probes and for each one, checks all unassigned probes for overlap via nested loop over full dataframe.
- Files: `eFISHent/optimization.py` (lines 139-242)
- Cause: For 100 assigned probes with 200 unassigned candidates, this is 20,000 overlap checks per alternative evaluation.
- Improvement path: Build spatial index (interval tree) of assigned probes once, then query it for overlaps. Reduce overlap checks from O(n²) to O(n log n).

## Fragile Areas

**Luigi Task Dependency Graph Complexity:**
- Files: `eFISHent/alignment.py` (lines 59-68), `eFISHent/cleanup.py` (lines 50-56)
- Why fragile: Multiple tasks declare overlapping dependencies (e.g., both alignment and cleanup need BuildJellyfishIndex). If a task's requires() is modified, downstream tasks may break silently or re-execute unnecessarily.
- Safe modification: Document task dependency graph (create ASCII diagram in ARCHITECTURE.md). Use type hints on requires() return values. Add integration tests that verify complete pipeline execution.
- Test coverage: `tests/test_integration.py` exists but doesn't verify incremental task re-execution or cache invalidation.

**Global State in Console Module:**
- Files: `eFISHent/console.py` (lines 38-46, 66, 434, 627, 669) - global variables `_silent_mode`, `_current_stage`, `_total_stages`, `_pipeline_progress`
- Why fragile: Multiple threads or parallel task execution could corrupt progress state. If pipeline runs twice in same process (unlikely but possible in tests), state is not reset.
- Safe modification: Refactor globals into a `PipelineProgressContext` class; use context managers to ensure cleanup. Add state reset in test fixtures.
- Test coverage: Unit tests for console don't verify concurrent access or state leakage between runs.

**Config Class Instantiation via `config()` Calls:**
- Files: Throughout codebase, e.g. `eFISHent/optimization.py` (lines 125-126, 154), `eFISHent/cleanup.py` (lines 86-87)
- Why fragile: `ProbeConfig()` instantiates a new config object each time. If Luigi cached config is not thread-safe or if config is modified mid-pipeline, different tasks may see different values.
- Safe modification: Cache config objects at pipeline start. Use immutable config frozen objects. Document that configs must not be modified after pipeline begins.
- Test coverage: No tests verify config consistency across multiple instantiations.

**Transcriptome Filtering Conditional Logic:**
- Files: `eFISHent/alignment.py` (lines 42-57), `eFISHent/cleanup.py` (lines 119-198)
- Why fragile: Transcriptome filtering is enabled only if both `reference_transcriptome` AND `reference_annotation` are provided. Logic is scattered across two modules. Missing either parameter silently disables transcriptome filtering with no warning.
- Safe modification: Add explicit validation in CLI to warn if user provided only one but not the other. Consolidate transcriptome path validation in a single function. Add debug logs showing which filters are active.

## Scaling Limits

**Bowtie/Bowtie2 Alignment Memory:**
- Current capacity: Tested with reference genomes ~3 GB (human, mouse)
- Limit: For genomes >5 GB, bowtie2 index creation may exceed available RAM or disk space. No checks in code.
- Scaling path: Warn users about memory requirements in docs. Detect available RAM at startup and bail out early if insufficient.

**Jellyfish K-mer Index Size:**
- Current capacity: Default k-mer length 31 requires index ~500 MB for human genome
- Limit: Very large genomes or custom k-mer lengths could consume all disk space
- Scaling path: Pre-calculate index size; warn user. Allow customization of k-mer length in config.

**DataFrame Merging with Large Expression Tables:**
- Current capacity: Expression tables tested up to 20k genes
- Limit: If expression table has millions of rows, pandas merge in `eFISHent/cleanup.py` (lines 97-100) will be slow
- Scaling path: Use indexed merge (ensure gene ID column is indexed). Consider chunking for very large tables.

## Dependencies at Risk

**Pyomo Solver Dependency (GLPK):**
- Risk: `eFISHent/optimization.py` (lines 340-342) uses `glpk` solver via Pyomo. If solver is not installed or version mismatch occurs, optimization fails at runtime.
- Impact: `--optimization-method optimal` fails silently or with confusing error if GLPK not found
- Migration plan: Add solver availability check in CLI (alongside dependency checks). Fall back to greedy method with warning if optimal unavailable. Or, use default solver via Pyomo (CBC if installed).

**BioPython Pairwise2 Performance:**
- Risk: `Bio.pairwise2.align.globalxx()` is slow for large probe counts. No alternative algorithm implemented.
- Impact: Analysis of large probe sets (>200 probes) takes hours
- Migration plan: Implement faster k-mer based similarity as alternative. Vendorize ssw (Smith-Waterman) library for speed.

## Missing Critical Features

**No Progress Checkpoint/Resume:**
- Problem: If pipeline crashes halfway (e.g., alignment timeout), user must restart from scratch. Luigi creates `.luigistate` file but it's not reliably checkpoint-able.
- Blocks: Cannot restart long-running jobs; forces re-alignment of probes already computed
- Recommendation: Implement explicit checkpoint files at each pipeline stage. Allow `--resume` flag to skip completed stages.

**No Quality Control Report:**
- Problem: Output CSV has per-probe quality scores (`eFISHent/cleanup.py`, line 110) but no overall pass/fail summary.
- Blocks: Users cannot quickly assess if probe set is usable without manual inspection
- Recommendation: Add summary stats (e.g., % probes with quality >= 80, coverage %, off-target count distribution) to final report.

**No Dry-Run Mode:**
- Problem: No way to validate inputs and configs without running full pipeline
- Blocks: Users cannot quickly test parameter combinations; must commit to full run
- Recommendation: Add `--dry-run` flag that validates inputs, checks dependencies, and prints what would be executed without running alignments.

## Test Coverage Gaps

**No Integration Test for Transcriptome Filtering:**
- What's not tested: Full pipeline with `--reference-transcriptome` flag. Downstream effects on cleanup stage when transcriptome hits are present.
- Files: `tests/test_integration.py` has no test marked with transcriptome param
- Risk: Refactoring transcriptome logic could break without detection
- Priority: Medium - feature is used but not end-to-end tested

**No Test for Bowtie2 Alignment Output Parsing:**
- What's not tested: SAM output from bowtie2 with various FASTA header formats (angle brackets, pipes, special chars). Only bowtie (v1) is implicitly tested.
- Files: `tests/test_unit_alignment.py` doesn't mock bowtie2 output
- Risk: Bowtie2 adoption could silently produce incorrect results
- Priority: High - bowtie2 is now primary aligner

**No Test for Concurrent Multiprocessing Scenarios:**
- What's not tested: `multiprocessing.Pool()` behavior in cleanup (secondary structure, k-mer counting) and optimization (Tm refinement). Tests run sequentially.
- Files: `tests/test_integration.py` runs with default threads=2 but doesn't vary thread count
- Risk: Race conditions, deadlocks could emerge with high thread counts
- Priority: Low - most users have 4-8 threads; unlikely to hit edge cases

**No Test for Out-of-Memory Conditions:**
- What's not tested: Behavior when k-mer index is corrupted, alignment index missing, or temp files run out of disk
- Files: No error handling tests in `tests/`
- Risk: Unclear error messages if resources exhausted
- Priority: Medium - affects production stability

**No Test for Invalid GTF/GFF Annotations:**
- What's not tested: Behavior when GTF is malformed, missing required columns, or has non-standard format
- Files: `eFISHent/gene_annotation.py` parses GTF via `gtfparse`, no validation of expected structure
- Risk: Silent data loss if GTF parsing fails
- Priority: Low - gtfparse is robust, but edge cases exist

---

*Concerns audit: 2026-03-14*
