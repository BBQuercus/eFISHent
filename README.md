<!-- [![Conda download statistics](https://anaconda.org/bioconda/efishent/badges/downloads.svg)]() -->
[![GitHub code licence is MIT](https://anaconda.org/bioconda/efishent/badges/license.svg)]()
[![Anaconda package version number](https://anaconda.org/bioconda/efishent/badges/version.svg)]()
[![Github Actions Status](https://github.com/bbquercus/eFISHent/workflows/Tests/badge.svg)]()
[![CodeFactor](https://www.codefactor.io/repository/github/bbquercus/efishent/badge)](https://www.codefactor.io/repository/github/bbquercus/efishent)
[![codecov](https://codecov.io/gh/BBQuercus/eFISHent/branch/main/graph/badge.svg?token=C1SRFYZ5VP)](https://codecov.io/gh/BBQuercus/eFISHent)
[![Maintainability](https://api.codeclimate.com/v1/badges/7470189fdd927276f80e/maintainability)](https://codeclimate.com/github/BBQuercus/eFISHent/maintainability)
[![DOI](https://zenodo.org/badge/501129295.svg)](https://zenodo.org/badge/latestdoi/501129295)

<img src="https://raw.githubusercontent.com/BBQuercus/eFISHent/main/logo.png" width="200px" align="right" alt="Logo of eFISHent.">

# eFISHent

A command-line based tool to facilitate the creation of eFISHent single-molecule RNA fluorescence in-situ hybridization (RNA smFISH) oligonucleotide probes.

## Contents

- [Description](#description)
- [Installation](#installation)
- [Getting Genomes, Annotations, Count Tables](#getting-genomes-annotations-count-tables)
- [Workflow](#workflow)
- [Usage](#usage)
- [Output](#output)
- [Full Examples](#full-examples)
- [FAQ](#faq)

## Description

eFISHent is a tool to facilitate the creation of eFISHent RNA smFISH oligonucleotide probes. Some of the key features of eFISHent are:

* One-line installation using conda (available through bioconda)
* Automatic gene sequence download from NCBI when providing a gene and species name (or pass a FASTA file)
* Filtering steps to remove low-quality probes including off-targets, frequently occurring short-mers, secondary structures, etc.
* Mathematical or greedy optimization to ensure highest coverage

## Installation

eFISHent is tested on macOS and Linux with Python 3.9+. Windows is not supported due to bioinformatics dependencies. For Windows users, we recommend [WSL](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support#1-overview) or a [Virtual Machine](https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview).

### Quick Install (Recommended)

Use the provided installation scripts to install all external dependencies:

```bash
# macOS
curl -sL https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install_macos.sh | bash

# Linux
curl -sL https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install_linux.sh | bash
```

Then add the activation script to your shell profile and install eFISHent:

```bash
# Add to ~/.zshrc (macOS) or ~/.bashrc (Linux)
echo 'source ~/.local/efishent-deps/activate.sh' >> ~/.zshrc  # or ~/.bashrc

# Reload shell
source ~/.zshrc  # or ~/.bashrc

# Install eFISHent
pip install efishent
```

### Custom Installation Path

To install dependencies to a custom location:

```bash
# Download and run with custom path
curl -sL https://raw.githubusercontent.com/BBQuercus/eFISHent/main/install_macos.sh -o install.sh
chmod +x install.sh
./install.sh /path/to/custom/location

# Then add the custom activation script to your profile
echo 'source /path/to/custom/location/activate.sh' >> ~/.zshrc
```

### Using Conda

Alternatively, use [conda](https://conda.io/) to manage dependencies:

```bash
# Create environment with dependencies
conda env create bbquercus/efishent

# Activate and install
conda activate efishent
pip install efishent
```

Updates can be done via pip: `pip install --upgrade efishent`

### Development Installation

To get the nightly build you can directly install eFISHent from this repository:

```bash
# Install other conda based dependencies
conda install -c bioconda bowtie blast jellyfish entrez-direct

# ONLY on linux!
conda install -c bioconda rnastructure

# Clone the repository
git clone https://github.com/BBQuercus/eFISHent.git
cd eFISHent/

# Install development version
pip install -e .
```

## Getting Genomes, Annotations, Count Tables

For some of the steps described below, you'll need to provide a few key resources unique to your model organism/cell line:

### Genome Sequence and Annotations

1. Go to the [UCSC genome browser](https://hgdownload.soe.ucsc.edu/downloads.html)
2. Find your organism, typically select the `Genome sequence files and select annotations` option
3. For the actual genome, download the file ending in `.fa.gz` (which will have to be unzipped e.g. using `gunzip`)
4. For the annotation, it's typically in the `genes/` directory ending in `.gtf.gz` (will also have to be unzipped)

### Count Table

1. Go to the [GEO dataset search](https://www.ncbi.nlm.nih.gov/gds/)
2. Search for your organism/cell line followed by `RNA-seq`
3. Select a sample that you think will represent your data the best (make sure it's RNA-seq - sometimes the search isn't the best...)
4. Scroll down and, if available, under `Supplementary file` select any files ending in `/FPKM/TPM/RPKMs.txt.gz` or similar. Do not download raw counts!
5. Alternatively you can also visit the [Expression Atlas](https://www.ebi.ac.uk/gxa/home) - count tables there might need some minor editing beforehand to ensure the required format of Ensembl ID and count values in columns 1 and 2 respectively

## Workflow

eFISHent works by iteratively selecting probes passing various filtering steps as outlined below:

1. A list of all candidate probes is generated from an input FASTA file containing the gene sequence. This sequence file can be passed manually or downloaded automatically from NCBI when providing a gene and species name.
2. The first round of filtering removes any probes not passing basic sequence-specific criteria including melting temperature as given formamide and salt concentrations, GC content, and G-quadruplets.
3. Probes are aligned to the reference genome using bowtie and candidates with off-targets are removed. In case of shorter genes or if off-targets are unavoidable, off-targets can be weighted using an encode count table to remove highly expressed genes.
4. The targets are divided into short k-mers and discarded if they appear above a determined threshold in the reference genome using Jellyfish.
5. The secondary structure of each candidate is predicted using a nearest neighbor thermodynamic model and filtered if the free energy is too high which could result in motifs hindering hybridization.
6. This gives the set of all viable candidates which are still overlapping. The final step is to use mathematical or greedy optimization to maximize probe non-overlapping coverage across the gene sequence.

<p align="center">
  <img src="https://github.com/BBQuercus/eFISHent/raw/main/workflow.png" width="500" />
</p>

## Usage

### Quick Start

```bash
eFISHent --reference-genome <reference-genome> --gene-name <gene> --organism-name <organism>
```

### Index Building

While there is only one main workflow, the slightly more time-intensive index creation step can be run ahead of time. Indexes are unique to each reference genome and can be created using:

```bash
eFISHent \
    --reference-genome <path to genome fasta file> \
    --build-indices True
```

### Passing Gene Sequence

The actual probe-creating workflow will then not only require the reference genome but also the sequence against which probes should be designed. Probes can be passed in one of three ways:

* `--sequence-file` - Path to a fasta file containing the gene sequence
* `--ensembl-id` (& `--organism-name`) - Ensembl ID of gene of interest. Will be downloaded from Entrez. The organism name can also be passed to avoid some wonky organism genes that have similar names but isn't required.
* `--gene-name` & `--organism-name` - Instead of ensembl ID, both gene and organism name can be provided. The sequence will also be downloaded from Entrez.

### Optimization Options

There are two ways in which the final set of probe candidates that passed filtering can be assigned/selected using `--optimization-method`:

* **greedy** - Uses the next best possibility in line starting with the first probe. This has a time complexity of `O(n)` with n being the number of candidates. Therefore, even with very loosely set parameters and a lot of candidates, this will still be very fast. This is the default option.
* **optimal** - Uses a mathematical optimization model to yield highest coverage (number of nucleotides bound to a gene). This has a time complexity of `O(n**2)` meaning the more probes there are the exponentially slower it will get. Despite breaking the problem into chunks, this might be restrictively slow. However, you can set a time limit (`--optimization-time-limit`) to stop the optimization process after a given amount of seconds. The resultant probes will be the best ones found so far.

### Off-Target Handling

To minimize the effect of off-targets, you can employ one of two strategies:

**Off-target minimization** - Using the maximum off-target flag, you can specify the maximum number of off-target bindings in the genome. By default this is set to zero meaning there aren't any known off-targets. However, for shorter or more repetitive genes, this might pose an issue which is why you can also use...

**Off-target weighting** - If off-targets are unavoidable, you can provide three parameters to select how high their expression is allowed to get:

* `--reference-annotation` - A GTF genome annotation file to know which genes correspond to which genomic loci
* `--encode-count-table` - A `csv` or `tsv` file with any normalized RNA-seq count table format (FPKM, FPKM, TPM, etc.) as well as the encode ID matching the entries in the GTF file
* `--max-expression-percentage` - The percentage of genes to be excluded sorted based on expression level (using the provided count table)

If you don't have your own RNA-seq dataset, you can download available datasets (make sure you're not using raw, but only normalized counts!) at [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/). Search for `RNA-seq` and the name of your organism/cell line.

### General Filtering Parameters

There are a bunch of parameters that can be set to adjust filtering steps:

| Parameter | Description |
|-----------|-------------|
| `--min-length`, `--max-length` | Probe lengths in nucleotides |
| `--spacing` | Minimum distances between probes |
| `--min-tm`, `--max-tm` | Minimum and maximum possible melting temperature (affected by length, GC content, and formamide/Na concentration) |
| `--min-gc`, `--max-gc` | GC content in percentage |
| `--formamide-concentration` | Percentage of formamide in buffer |
| `--na-concentration` | Sodium ion concentration in mM |
| `--kmer-length`, `--max-kmers` | Jellyfish-based short-mer filtering. If candidate probes contain kmers of length `--kmer-length` that are found more than `--max-kmers` in the reference genome, the candidate will get discarded |
| `--max-deltag` | Predicted secondary structure threshold |
| `--sequence-similarity` | Will remove probes that might potentially bind to each other (set the similarity to the highest allowed binding percentage). Will reduce the number of probes because (due to otherwise exceptionally high runtimes) has to be run after optimization |

### Remaining Options

There are a few additional options:

| Parameter | Description |
|-----------|-------------|
| `--is-plus-strand` | Set true/false depending on gene of interest |
| `--is-endogenous` | Set true/false depending on gene of interest |
| `--threads` | Wherever multiprocessing is available, spawn that many threads. Set this to as many cores as you have available |
| `--save-intermediates` | Save all intermediary files. Can be used to gauge which filtering steps are set too aggressively |
| `--verbose` | Set to get some more information on progress |

### Probe Set Analysis

After generating a probe set, you can analyze it in more detail using the `--analyze-probeset` option. This creates a PDF report with visualizations of various probe characteristics:

```bash
eFISHent \
    --reference-genome <path to genome fasta file> \
    --sequence-file <path to gene fasta file> \
    --analyze-probeset <path to probe set fasta file>
```

The analysis includes:

| Plot | Description |
|------|-------------|
| Lengths | Distribution of probe lengths |
| Melting temperatures | Boxplot of calculated Tm values |
| GC Content | Boxplot of GC percentages |
| G quadruplet | Count of G-quadruplet motifs per probe |
| K-mer count | Maximum k-mer frequency in genome |
| Free energy | Predicted secondary structure stability (ΔG) |
| Off target count | Number of off-target binding sites per probe |
| Binding affinity | Probe-to-probe similarity matrix (potential cross-hybridization) |
| Gene coverage | Visual map of probe positions along the target sequence |

The output is saved as `<probeset_name>_analysis.pdf` in the current directory.

## Output

By default eFISHent will output three unique files:

* `GENE_HASH.fasta` - All probes in FASTA format for subsequent usage
* `GENE_HASH.csv` - A table containing all probes as well as basic parameters (such as melting temperature)
* `GENE_HASH.txt` - A configuration file to check which parameters were used during the run as well as the command to start it

`GENE` is a reinterpreted gene name dependent on the options passed but should be immediately clear as to where it's from. The `HASH` is a unique set of characters that identifies the parameters passed for the run. This way, if and only if the same parameters are passed again eFISHent doesn't have to rerun anything. All intermediary files during the run will be saved in the same format but will get deleted at the end unless `--save-intermediates` is set to true.

## Full Examples

First, the indexes for the respective genome have to be built:

```bash
eFISHent \
    --reference-genome ./hg-38.fa \
    --build-indices True
```

An example to create 45 to 50-mers for a gene of interest downloaded from Entrez:

```bash
eFISHent \
    --reference-genome ./hg-38.fa \
    --gene-name "norad" \
    --organism-name "homo sapiens" \
    --is-plus-strand True \
    --optimization-method optimal \
    --min-length 45 \
    --max-length 50 \
    --formamide-concentration 45 \
    --threads 8
```

Another example using a custom sequence:

```bash
eFISHent \
    --reference-genome ./dm6.fa \
    --sequence-file "./renilla.fasta" \
    --is-endogenous False \
    --threads 8
```

An example with off-target weighting:

```bash
eFISHent \
    --reference-genome ./hg-38.fa \
    --reference-annotation ./hg-38.gtf \
    --ensembl-id ENSG00000128272 \
    --organism-name "homo sapiens" \
    --is-plus-strand False \
    --max-off-targets 5 \
    --encode-count-table ./count_table.tsv \
    --max-expression-percentage 20 \
    --threads 8
```

Lastly, an example to analyze an existing probe set:

```bash
eFISHent \
    --reference-genome ./hg-38.fa \
    --sequence-file ./my_gene.fasta \
    --analyze-probeset ./my_gene_probes.fasta
```

## FAQ

Have questions? Open an issue on [GitHub](https://github.com/BBQuercus/eFISHent/issues).
