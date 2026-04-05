# Changelog

## v0.0.12

### Bug Fixes
- Fix rDNA/satellite off-target blind spots for exogenous probes: probes are now screened against 45S rDNA (U13369.1) and centromeric alpha satellite consensus sequences, which are absent from the GRCh38 primary assembly but can cause intense nucleolar/centromeric background signal
- Bundle `Fold_linux` binary in the package and sync `uv.lock`
- Remove unused imports flagged by ruff

## v0.0.11

### Bug Fixes
- Fix secondary structure filtering crash on Linux: use bundled `Fold_linux` binary instead of expecting `Fold` on system PATH

## v0.0.10

### New Features

#### Automated Probe Validation Report (B1/B2/B3)
- Output CSV now includes transcriptome off-target details:
  - `txome_off_targets`: count of unique off-target transcripts (≥90% identity / ≥18bp)
  - `off_target_genes`: comma-separated `GENE(count)` list mapping transcript IDs to gene names via GTF
  - `worst_match`: best off-target match quality string (e.g., `95%/20bp/0mm`)
  - `recommendation`: automated PASS / FLAG / FAIL recommendation per probe
  - `expression_risk`: expression-level risk annotation when count table is provided
- New `eFISHent/gene_annotation.py` module for transcript-to-gene mapping supporting Ensembl, RefSeq, gffread, and pipe-delimited ID formats
- Final probes are re-BLASTed against the transcriptome for detailed reporting (even if transcriptome filtering was used during candidate selection)

#### Adaptive Probe Length (C1)
- New `--adaptive-length` flag adjusts probe length based on local GC content:
  - High-GC regions (>55%): prefer shorter probes to lower Tm
  - Low-GC regions (<45%): prefer longer probes to raise Tm
  - Balanced GC (45-55%): use middle length
- Enabled by default in the `smfish` preset (now generates 18-22nt probes instead of fixed 20nt)

#### Tm Uniformity Optimization (C2)
- After probe selection, a refinement step swaps Tm outliers (>5°C from set median) for better alternatives at nearby positions
- Produces more uniform Tm across the final probe set

#### Cumulative Off-Target Cap (C3)
- New `--max-probes-per-off-target` parameter (default: 0 = disabled)
- When enabled, iteratively removes lowest-quality probes when too many in the final set hit the same off-target gene
- Prevents correlated off-target binding that creates false FISH spots (Stellaris-style safety net)

#### 16-nt Cross-Hybridization Warning (C4)
- Probes with ≥16nt contiguous match at ≥95% identity to off-target transcripts are flagged in debug logs
- Warning-level only (not hard rejection), complements the existing BLAST filtering

#### Exogenous Gene Support
- **K-mer filter bypass**: When `--is-endogenous False`, k-mer filtering is skipped entirely — short k-mer matches in the host genome are random chance for foreign sequences, and the transcriptome BLAST filter catches real off-targets
- **New `exogenous` preset**: `--preset exogenous` sets optimal defaults for GFP, Renilla, mCherry, and other reporter genes (19-22nt, adaptive length, strict BLAST, no k-mer filter)
- **Smart warnings**: CLI warns when exogenous mode is used without `--reference-transcriptome`
- This change alone increases probe yield for exogenous genes from ~11 to ~30+ probes (previously k-mer filter eliminated 97% of candidates)

#### Improved Greedy Probe Selection
- Greedy optimizer now sorts by start position with length as tiebreaker, maximizing probe count while preferring longer probes at each position
- For Renilla: increased from 29 to 35 probes (better than Stellaris's 34)

### Improvements

#### Quality Score Reweighting (C5)
- Updated to Stellaris-informed weights: Tm=25%, GC=15%, ΔG=20%, off-target=25%, k-mer=15%
- Tm score now uses the **actual median Tm of the selected set** instead of config midpoint — penalizes outliers within the set
- Off-target score combines genome and transcriptome off-target counts

#### Configurable BLAST Match Length (A1)
- New `--min-blast-match-length` CLI parameter replaces the hardcoded `min_match_len = 15`
- Default: `max(18, 0.8 * min_probe_length)` — less aggressive for exogenous genes
- Fixes the issue where the transcriptome filter rejected nearly all probes for exogenous genes due to chance 15bp partial matches

#### Updated smFISH Preset
- Now uses 18-22nt probe range with adaptive length enabled (was fixed 20nt)

### Infrastructure Fixes

#### makeblastdb Installation (A2)
- `install.sh --with-blast` now installs `makeblastdb` alongside `blastn` and `dustmasker`
- Added to dependency checks and required-deps list when `--reference-transcriptome` is used
- Fixes silent failure of transcriptome filter when `makeblastdb` was missing

#### gffread Installation (A5)
- `install.sh --with-blast` now also installs `gffread` for building transcriptome FASTA from genome + GTF
- Downloaded from GitHub releases as a single static binary

#### Summary Stats
- Tm standard deviation now included in completion summary

---

## v0.0.9

- Fix off-target filtering bugs: exogenous rRNA gate, thermo rescue pre-filter, substring matching
- Add human quick-start guide with Ensembl download instructions
- Add probe quality scoring, rRNA filtering, and updated presets
- Fix installer issues: BLAST+, Fold, PATH, gffread, version string
