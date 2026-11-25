# Changelog

## Unreleased

- **Checksums for pre-built genomes**: `--genome` / `--download-genome` now check every downloaded file against the SHA256 in the genome's `metadata.json`. On a mismatch the download stops with an error and the bad files are deleted, so the next run downloads them again. A file without a checksum entry is kept with a warning. A supplemental `genome_<k>.jf` that fails its checksum is deleted and the index is built locally.

## v0.0.19

**Probe sets can change in this release.** The same command can give a different probe set than v0.0.18; on the benchmark below the new sets are more specific. No command-line option name, command-line default, column of the final probe table or Python package dependency changed. The internal BLAST settings and the prepared annotation `.parq` did change (see below).

**To get the new probe set for an earlier design, use a fresh `--output-dir`.** The output name hash does not include the eFISHent version, so the same command in the same directory reuses the v0.0.18 result and its cached intermediate files.

### Off-target Filtering and Probe Selection
- **More sensitive transcriptome BLAST**: with `-max_target_seqs 50 -evalue 100`, short perfect matches crowded out the 18+ nt one-mismatch hits that the off-target rule counts. The filter and the final report now use evalue 1000 and 1000 targets (TrueProbes settings). On 1,392 Renilla candidates, the old settings missed qualifying off-targets on 257 probes
- **Strand-aware hits**: only plus/minus hits (probe antisense to the transcript) count as off-targets. Plus/plus hits, e.g. to an overlapping antisense RNA, cannot bind and no longer reject probes. This assumes a transcriptome FASTA in transcript orientation (gffread `-w`, Ensembl cDNA)
- **Same-locus rule counts only perfect hits of at least 15 nt**: a transcript that matches at least half of the candidate probes perfectly is treated as the target's own locus (overlapping gene model, read-through) and not as an off-target. Shorter perfect matches no longer count. With the more sensitive BLAST they would whitelist other transcripts on a small transcriptome (yeast AAD4: the paralog AAD6 and six long unrelated genes)
- **Off-target burden**: every non-self hit ≥15 nt gets a duplex Tm and a binding weight. The per-probe burden ranks candidates during selection and is written to `GENE_HASH_txome_burden.csv` when `--save-intermediates` is set
- **`--off-target-min-tm` now also applies to the transcriptome filter**: before, it only changed the genome count gate, which is skipped for runs with a transcriptome and `--max-off-targets 0`. The default stays 0 (disabled). The `relaxed` preset already sets it to 37, so `relaxed` runs with a transcriptome now leave weak hits out of the transcriptome off-target count without an explicit flag. The separate cross-hybridization rejection still sees those hits. On the four Drosophila benchmark targets, 37 added 1 to 6 probes but let 2 to 8 off-target probes back in; the three exogenous targets did not change. It remains a coverage/specificity trade-off
- **Genome screen before selection (exogenous targets)**: probes with a near-full-length genomic BLAST hit were removed after optimization and nothing replaced them (Renilla: 36 → 20 probes). The same check, with the `--max-off-targets` allowance, now also runs on all candidates before selection, so the optimizer can backfill. The final verification is kept
- **Exact probe selection**: `--optimization-method greedy` (the default) now makes the initial selection by exact weighted interval scheduling over a combined weight (Tm refinement can still swap probes afterwards): every probe counts 1, plus up to 0.5 for a low off-target burden (rank among the candidates) and up to 0.25 for sequence quality. The old left-to-right pass took whichever probe started first. The new selection favours more probes, but it does not strictly maximize the probe count: it can prefer a slightly smaller set of cleaner probes. `optimal` (GLPK) now also discounts each probe's quality weight by its burden rank

Benchmark on 7 targets (Renilla, EGFP, mCherry against human; dpn, Oli, CycB, Syp in Drosophila), every set scored by one independent BLAST evaluator (antisense hits only, hits to the target gene excluded). The table compares the probe-selection code of v0.0.18 with the reworked filtering and selection (mean per target). These runs were made before the GTF streaming and the internal module split landed; Renilla and dpn were re-run on the code after those two changes and gave identical probe sets. The same-locus rule above was tightened after that; for Renilla, mCherry and dpn it whitelists the same transcripts as before, the other four targets were not re-checked:

| Metric | Before | After |
|--------|--------|-------|
| Probes with an off-target (≥90% identity over ≥18 nt) | 4.0 | 0 |
| Target coverage | 51.9% | 62.8% |
| Probes | 47.4 | 52.1 |

Yeast (test genome with transcriptome and annotation, same evaluator), v0.0.18 → v0.0.19:

| Target | Probes | Coverage | Probes with an off-target |
|--------|--------|----------|---------------------------|
| AAD4 | 6 → 6 | 14.1% → 13.8% | 0 → 0 |
| ACT1 | 39 → 38 | 78.3% → 75.4% | 1 → 0 |
| CDC19 | 56 → 56 | 83.6% → 81.1% | 2 → 0 |
| HSP104 | 95 → 90 | 80.0% → 74.8% | 3 → 0 |
| TDH3 | 37 → 21 | 83.2% → 45.4% | 33 → 16 |
| Renilla (exogenous) | 39 → 39 | 83.3% → 81.1% | 0 → 0 |

### Behaviour Changes
- A target with paralogs can get clearly fewer probes: yeast TDH3 37 → 21. In v0.0.18, 33 of the 37 probes also matched TDH2 and 24 matched TDH1 (final report). In v0.0.19, none of the 21 probes matches TDH1, and 16 still match TDH2: a very close paralog can pass the same-locus rule, and the candidate filter then does not reject probes for it. The final probe table still counts these hits in `txome_off_targets` and lists the genes in `off_target_genes`
- Two of the seven benchmark targets get slightly fewer probes: CycB 86 → 84 (coverage 73.2% → 69.5%) and Syp 100 → 97 (71.1% → 67.6%). Their v0.0.18 sets contained 3 and 10 probes with an off-target; the new sets contain none
- Runs with `--reference-transcriptome` take longer because of the more sensitive BLAST: 17 to 41 s more end-to-end in the same seven benchmark runs (71–185 s before, 88–225 s after, 6 threads). On the yeast test genome the run needs about 1 GB more memory (0.5–0.9 GB before, 1.4–2.0 GB after)
- `--max-expression-percentage` can pick a different cutoff, and so remove different probes, once the annotation `.parq` is regenerated: the cutoff is a quantile over annotation rows, and the prepared annotation no longer contains CDS/UTR/codon rows (see below)

### Performance
- **GTF annotations are streamed into a slim parquet** (gene, transcript and exon rows; only the attributes eFISHent reads). Peak memory for the human Ensembl 115 GTF (3.3 GB) is now 558 MB; before, it exceeded a 20 GB machine. C. elegans: 1.96 GB → 262 MB and 9.9 s → 3.4 s. A `.parq` file prepared by an older version is still read as is

### Installer
- `install.sh` works on Linux without `unzip`, and on Linux x86-64 without a C compiler: zip archives fall back to python3, and jellyfish uses the pre-built 2.3.0 binary (other platforms still build from source). A missing `libgomp` for BLAST+ is fetched without root where `apt-get` and `dpkg-deb` exist; otherwise the installer prints an install hint
- The BLAST+ verify step now runs `blastn -version` instead of only checking that the file exists

### Other
- The PyPI update check only follows https URLs, redirects included
- Internal: large modules were split and lint rules pinned, with identical output

## v0.0.13 – v0.0.18

See the [GitHub releases](https://github.com/BBQuercus/eFISHent/releases).

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
