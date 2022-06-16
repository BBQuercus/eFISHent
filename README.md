# eFISHent

A tool to facilitate the creation of eFISHent RNA FISH oligonucleotide probes.

## Description

eFISHent is a tool to facilitate the creation of eFISHent RNA FISH oligonucleotide probes. Some of the key features of eFISHent are:

* Single-line installation using conda
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

A detailed usage guide can be found on the [GitHub wiki]() but here is a quick example:

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
  * [ ] Add mathematical description for model
  * [x] Add basic documentation in CLI and in README.md
  * [x] Add more detailed documentation as wiki page(s)
  * [ ] Create bioconda recipe

* **Components**
  * [x] Add spacing option (length +- nt in optimization step)
  * [ ] Allow selection of isoform in entrez query
  * [ ] Select intron/exon based on blast from sequence
  * [ ] Verify bowtie parameters run on endo and exo targets
  * [ ] Remove gene of interest when using rna-seq data
  * [ ] Save alignment and kmer scores to csv output

* **Interface**
  * [x] Clean up output files
  * [x] Clean up logging and error handling
  * [x] Decide on way to handle temporary files (tempdir?)
  * [x] Find way to handle rerunning same gene with different parameters (unique name hash?)
  * [x] Find way to make CLI alter config (luigi.configuration.add_config_path)

* **Testing**
  * [ ] Add unit tests for core components
    * Sliding window probe generation
    * Entrez query error handling
    * Parsing alignment output
    * Parsing blast output
    * Selecting intron/exon
    * RNA-seq data loading
    * RNA-seq data gene of interest removal
    * RNA-seq data sorting
    * Secondary structure prediction
    * Optimization (small mock cases) - mathematical model
    * Optimization - greedy model
  * [ ] Add integration tests for the following examples
    * dm6/hr38 (download from entrez, endo, long, multiple isoforms)
    * dm6/kayak (provided, select introns, endo)
    * dm6/Renilla (provided, exo)
    * hg38/RPL12 (provided, endo, repetitive, short, minus strand)
    * hg38/VIM (provided, endo, plus strand)
