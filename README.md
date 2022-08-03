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

## Description

eFISHent is a tool to facilitate the creation of eFISHent RNA smFISH oligonucleotide probes. Some of the key features of eFISHent are:

* One-line installation using conda (available through bioconda)
* Automatic gene sequence download from NCBI when providing a gene and species name (or pass a FASTA file)
* Filtering steps to remove low-quality probes:
  * Basics such as melting temperature, GC content, and length
  * Off-target bindings (based on raw count or weighted through RNAseq count tables)
  * Frequently occuring short-mers in the genome
  * Predicted free energy in the secondary structure
* Mathematical or greedy optimization to ensure highest coverage

## Installation

eFISHent is being tested on MacOS and Linux. Unfortunately, due to the bioinformatics dependencies Windows is not supported. For Windows users, we reccommend installing "Windows Subsystem for Linux (WSL)" ([Windows 10](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#1-overview), [Windows 11](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support#1-overview)) or using a fully fledged [Virtual Machine](https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview). eFISHent can be installed using the [conda](https://conda.io/) package manager.

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

* Bowtie index (`indexing.py`)
* Blast database (`indexing.py`)
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

Probe set analysis plotting:

* Create a simple overview over the key parameters (`analyze.py`)
* (Using inputs from most other files)

## TODO

* [ ] Add more detailed documentation as wiki page(s)
  * [x] Add links to genomes and RNAseq databases
  * [x] Add examples from multiple sources
  * [ ] Add benchmarks for deltaG, counts
* [ ] Add mathematical description for model (in wiki?)
* [x] Add filtering step to prevent off chance of probes hybridizing with themselves
  * [x] Create function to check if probes could bind
  * [x] Greedy model - skip if probe rev complement is in set
  * [x] Optimal model - add constraint between similar probes to not get assigned
* [x] Add probe set analysis visualization
  * [x] Input probe set fasta
  * [x] All filtering step components as graphs (tm, gc, length, # off targets)
  * [x] Save output as pdf?
* [ ] Add probe set analysis txt file with off-target locations / potentially harmful probes
