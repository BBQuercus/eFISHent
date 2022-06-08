# eFISHent

A tool to facilitate the creation of eFISHent RNA FISH oligonucleotide probes.

## TODO

- **General**
  - [x] Find a name and description for the tool
  - [ ] Add mathematical description for model
  - [ ] Add documentation in CLI and in README.md
  - [ ] Create bioconda recipe

- **Components**
  - [ ] Allow selection of isoform in entrez query
  - [ ] Select intron/exon based on blast from sequence
  - [ ] Verify bowtie parameters run on endo and exo targets
  - [ ] Remove gene of interest when using rna-seq data

- **Interface**
  - [ ] Clean up output files
  - [ ] Clean up logging and error handling
  - [ ] Decide on way to handle temporary files (tempdir?)
  - [ ] Find way to handle rerunning same gene with different parameters (unique name hash?)
  - [x] Find way to make CLI alter config (luigi.configuration.add_config_path)

- **Testing**
  - [ ] Add unit tests for core components
    - Sliding window probe generation
    - Entrez query error handling
    - Parsing alignment output
    - Parsing blast output
    - Selecting intron/exon
    - RNA-seq data loading
    - RNA-seq data gene of interest removal
    - RNA-seq data sorting
    - Secondary structure prediction
    - Optimization (small mock cases) - mathematical model
    - Optimization - greedy model
  - [ ] Add integration tests for the following examples
    - dm6/hr38 (download from entrez, endo, long, multiple isoforms)
    - dm6/kayak (provided, select introns, endo)
    - dm6/Renilla (provided, exo)
    - hg38/RPL12 (provided, endo, repetitive, short, minus strand)
    - hg38/VIM (provided, endo, plus strand)
