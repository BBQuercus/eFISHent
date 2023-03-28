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

* One-line installation using conda (available through bioconda*)
* Automatic gene sequence download from NCBI when providing a gene and species name (or pass a FASTA file)
* Filtering steps to remove low-quality probes including off-targets, frequently occuring short-mers, secondary structures, etc.
* Mathematical or greedy optimization to ensure highest coverage

\* The release on bioconda is always associated with waiting times. Therefore, the easiest approach is to install conda dependencies and install eFISHent using pip.

## Installation

eFISHent is being tested on MacOS and Linux with Python versions 3.8 - 3.10. Unfortunately, due to the bioinformatics dependencies Windows is not supported. For Windows users, we reccommend installing "Windows Subsystem for Linux (WSL)" ([Windows 10](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#1-overview), [Windows 11](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support#1-overview)) or using a fully fledged [Virtual Machine](https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview). Using [conda](https://conda.io/) environment, install eFISHent as follows:

```bash
# Create an environment and install all dependencies (e.g. python)
conda env create bbquercus/efishent  

# Activate environment
conda activate efishent

# Install efishent via pypi
pip install efishent
```

Any updates can then simply be done via pypi (`pip install --upgrade efishent`).

## Usage

A detailed usage guide can be found on the [GitHub wiki](https://github.com/bbquercus/eFISHent/wiki) but here is a quick example:

```bash
eFISHent --reference-genome <reference-genome> --gene-name <gene> --organism-name <organism>
```

## Component overview

eFISHent is built up modularly using the following components...

Index creation workflow:

* Bowtie index
* Jellyfish indices

Probe filtering workflow:

* Download / prepare sequences
* Generate candidate probes
* Filter with basic filters
* Align probes to reference genome
* Filter based on alignment score and uniqueness
* Filter reoccuring k-mers
* Filter based on secondary structure prediction
* Create final list of probes
* Write final list of probes to file with report

Probe set analysis plotting:

* Create a simple overview over the key parameters

## TODO

* [ ] Add more detailed documentation as wiki page(s)
  * [x] Add links to genomes and RNAseq databases
  * [x] Add examples from multiple sources
  * [ ] Add benchmarks for deltaG, counts
* [ ] Add mathematical description for model (in wiki?)
* [ ] Add probe set analysis txt file with off-target locations / potentially harmful probes
