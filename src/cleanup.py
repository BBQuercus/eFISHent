"""
Clean up the output directory and remove the files that are not needed.
Prettify the output files.

Output files:
    - .fasta with the probe sequences
    - .csv with probe information (Tm, GC, deltaG, etc.)
"""

import os

import luigi
import pandas as pd
import Bio.SeqIO

class CleanUpOutput(luigi.Task):
    """"""
    