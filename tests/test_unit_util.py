import configparser
import os
import shutil
import tempfile

import Bio.SeqIO
import luigi
import pytest

from eFISHent.util import create_config_hash
from eFISHent.util import create_data_table
from eFISHent.util import get_gene_name
from eFISHent.util import get_genome_name
from eFISHent.util import get_output_dir


def test_create_config_hash():
    # Original
    luigi.configuration.add_config_path("./eFISHent/luigi.cfg")
    hash1 = create_config_hash(luigi.configuration.get_config())

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Identical in different place
        config_file = os.path.join(tmp_dir, "luigi.cfg")
        shutil.copy("./eFISHent/luigi.cfg", config_file)
        luigi.configuration.add_config_path(config_file)
        hash2 = create_config_hash(luigi.configuration.get_config())

        # Single modification
        config = configparser.ConfigParser()
        config.read(config_file)
        config.set(
            "GeneralConfig", "threads", str(int(config["GeneralConfig"]["threads"]) + 1)
        )
        with open(config_file, "w") as f:
            config.write(f)
        luigi.configuration.add_config_path(config_file)
        hash3 = create_config_hash(luigi.configuration.get_config())

    assert hash1 == hash2
    assert hash1 != hash3


def test_create_data_table():
    sequences = list(Bio.SeqIO.parse("./tests/data/renilla_data_table.fasta", "fasta"))
    df = create_data_table(sequences)

    for col in ["name", "length", "sequence", "start", "end"]:
        assert col in df.columns

    raw_sequences = [seq.seq for seq in sequences]
    for seq in df["sequence"]:
        assert seq in raw_sequences


def test_gene_name_genefile():
    class Config(luigi.Config):
        sequence_file = luigi.Parameter("./tests/data/renilla.fasta")
        ensembl_id = luigi.Parameter("")
        gene_name = luigi.Parameter("")
        organism_name = luigi.Parameter("")

    assert get_gene_name(config=Config).startswith("renilla")
    assert get_gene_name(hashed=False, config=Config) == "renilla"


def test_gene_name_ensembl():
    class Config(luigi.Config):
        sequence_file = luigi.Parameter("")
        ensembl_id = luigi.Parameter("ENSG00000026025")
        gene_name = luigi.Parameter("")
        organism_name = luigi.Parameter("")

    assert get_gene_name(config=Config).startswith("ENSG00000026025")


def test_gene_name_genename():
    class Config(luigi.Config):
        sequence_file = luigi.Parameter("")
        ensembl_id = luigi.Parameter("")
        gene_name = luigi.Parameter("ACOOLGene123")
        organism_name = luigi.Parameter("Latin name")

    assert get_gene_name(config=Config).startswith("Latin_name_ACOOLGene123_")


def test_gene_name_genefile_genename():
    class Config(luigi.Config):
        sequence_file = luigi.Parameter("./tests/data/renilla.fasta")
        ensembl_id = luigi.Parameter("")
        gene_name = luigi.Parameter("ACOOLGene123")
        organism_name = luigi.Parameter("Latin name")

    assert get_gene_name(config=Config).startswith("renilla")


def test_genome_name_notpassed():
    class Config(luigi.Config):
        reference_genome = luigi.Parameter("")

    with pytest.raises(ValueError):
        get_genome_name(config=Config)


def test_genome_name_invalid():
    class Config(luigi.Config):
        reference_genome = luigi.Parameter("./some_random_file.notfasta")

    with pytest.raises(ValueError):
        get_genome_name(config=Config)


def test_genome_name_valid():
    class Config(luigi.Config):
        reference_genome = luigi.Parameter("./tests/data/sacCer3.fa")

    assert get_genome_name(config=Config) == os.path.abspath("./tests/data/sacCer3")


def test_output_dir_passed():
    # Passed
    class Config(luigi.Config):
        output_dir = luigi.Parameter("./tests/data/")

    assert get_output_dir(config=Config) == os.path.abspath("./tests/data/")


def test_output_dir_unpassed():
    class Config(luigi.Config):
        output_dir = luigi.Parameter("")

    assert get_output_dir(config=Config) == os.getcwd()
