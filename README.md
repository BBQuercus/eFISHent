<!-- [![Github Actions Status](https://github.com/bbquercus/efishent/workflows/main/badge.svg)]() -->
<!-- [![Conda download statistics](https://anaconda.org/bioconda/efishent/badges/downloads.svg)]() -->
[![GitHub code licence is MIT](https://anaconda.org/bioconda/efishent/badges/license.svg)]()
[![Anaconda package version number](https://anaconda.org/bioconda/efishent/badges/version.svg)]()
[![DOI](https://zenodo.org/badge/501129295.svg)](https://zenodo.org/badge/latestdoi/501129295)

<img src="https://github.com/bbquercus/efishent/raw/master/logo.png" width="200px" align="right" alt="Logo of deepBlink.">

# eFISHent

A tool to facilitate the creation of eFISHent RNA FISH oligonucleotide probes.

## Description

eFISHent is a tool to facilitate the creation of eFISHent RNA FISH oligonucleotide probes. Some of the key features of eFISHent are:

* One-line installation using conda (available through bioconda)
* Automatic gene sequence download from NCBI when providing a gene and species name (or pass a FASTA file)
* Filtering steps to remove low-quality probes:
  * Basics such as melting temperature, GC content, and length
  * Off-target bindings (based on raw count or weighted through RNAseq count tables)
  * Frequently occuring short-mers in the genome
  * Predicted free energy in the secondary structure
* Mathematical or greedy optimization to ensure highest coverage

## Installation

eFISHent can be installed using the [conda](https://conda.io/) package manager.

```bash
conda install -c bioconda efishent
```

## Usage

A detailed usage guide can be found on the [GitHub wiki](https://github.com/bbquercus/eFISHent/wiki) but here is a quick example:

```bash
eFISHent --reference-genome <reference-genome> --gene-name <gene> --organism-name <organism>
```

## Component overview

eFISHent is built up modularly using the following components...

Index creation workflow:

* Bowtie index (`alignment.py`)
* Jellyfish indices (`kmers.py`)

Probe filtering workflow:

* Download / prepare sequences (`prepare_sequences.py`)
* Generate candidate probes (`generate_probes.py`)
* Filter with basic filters (`basic_filtering.py`)
* Align probes to reference genome (`alignment.py`)
* Filter based on alignment score and uniqueness (`alignment.py`)
* Filter reoccuring k-mers (`kmers.py`)
* Filter based on secondary structure prediction (`secondary_structure.py`)
* Create final list of probes (`optimization.py`)
* Write final list of probes to file with report (`cleanup.py`)

## TODO

* **General**
  * [x] Find a name and description for the tool
  * [x] Add basic documentation in CLI and in README.md
  * [x] Add thread limiting check with `os.cpu_count()`
  * [x] Create bioconda recipe
  * [x] Add logo to repository
  * [ ] Add more detailed documentation as wiki page(s)
    * [ ] Add benchmarks for deltaG, FPKM
    * [ ] Add links to genomes and RNAseq databases
    * [x] Add examples from multiple sources
  * [ ] Add mathematical description for model (in wiki?)
  * [x] Add citation
  * [ ] Set up CI pipeline
