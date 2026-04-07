import os
from unittest.mock import patch, MagicMock

import luigi
import pandas as pd
import pytest

from eFISHent.indexing import (
    BuildBowtieIndex,
    BuildBowtie2Index,
    PrepareAnnotationFile,
)


@pytest.fixture
def task_prepare():
    return PrepareAnnotationFile()


@pytest.mark.skip(reason="gtfparse uses deprecated polars API (toggle_string_cache)")
def test_prepare_annotations(task_prepare):
    fname_output = "./mylittle.gtf.parq"
    task_prepare.prepare_gtf_file("./tests/data/sacCer3.gtf", fname_output)

    assert os.path.isfile(fname_output)

    df = pd.read_parquet(fname_output)
    for col in ["seqname", "source", "feature", "start", "end"]:
        assert col in df.columns
    os.remove(fname_output)


class TestBuildBowtieIndexOutput:
    """Test BuildBowtieIndex output targets."""

    def test_output_returns_six_targets(self):
        """Bowtie index should produce 6 ebwt files."""
        with patch("eFISHent.util.get_genome_name", return_value="/tmp/test_genome"):
            task = BuildBowtieIndex()
            outputs = task.output()
            assert len(outputs) == 6

    def test_output_extensions(self):
        """All 6 bowtie index extensions should be present."""
        with patch("eFISHent.util.get_genome_name", return_value="/tmp/test_genome"):
            task = BuildBowtieIndex()
            outputs = task.output()
            paths = [o.path for o in outputs]
            expected_exts = [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt",
                             ".rev.1.ebwt", ".rev.2.ebwt"]
            for ext in expected_exts:
                assert any(p.endswith(ext) for p in paths), f"Missing {ext}"

    def test_output_uses_genome_name(self):
        """Output paths should be based on get_genome_name."""
        with patch("eFISHent.util.get_genome_name", return_value="/data/my_genome"):
            task = BuildBowtieIndex()
            outputs = task.output()
            for o in outputs:
                assert o.path.startswith("/data/my_genome")


class TestBuildBowtie2IndexOutput:
    """Test BuildBowtie2Index output target."""

    def test_output_returns_single_target(self):
        """Bowtie2 index output should be a single .1.bt2 file."""
        with patch("eFISHent.util.get_genome_name", return_value="/tmp/test_genome"):
            task = BuildBowtie2Index()
            output = task.output()
            assert isinstance(output, luigi.LocalTarget)
            assert output.path == "/tmp/test_genome.1.bt2"


class TestPrepareAnnotationFileOutput:
    """Test PrepareAnnotationFile output target."""

    def test_output_is_parquet(self):
        """Output should be a .parq file based on reference_annotation."""
        with patch("eFISHent.indexing.GeneralConfig") as mock_config:
            mock_config.return_value.reference_annotation = "/data/genome.gtf"
            task = PrepareAnnotationFile()
            output = task.output()
            assert isinstance(output, luigi.LocalTarget)
            assert output.path == "/data/genome.gtf.parq"


class TestBuildBowtieIndexRun:
    """Test BuildBowtieIndex run method without actually calling bowtie."""

    def test_build_bowtie_index_calls_subprocess(self):
        """build_bowtie_index should call bowtie-build with correct args."""
        task = BuildBowtieIndex()
        with patch("eFISHent.indexing.subprocess.check_call") as mock_call, \
             patch("eFISHent.indexing.BuildBowtieIndex.logger"), \
             patch("eFISHent.console.spinner", return_value=MagicMock(__enter__=MagicMock(), __exit__=MagicMock())):
            task.build_bowtie_index("/path/to/genome.fa", "/path/to/genome")
            mock_call.assert_called_once()
            args = mock_call.call_args[0][0]
            assert args[0] == "bowtie-build"
            assert args[1] == "/path/to/genome.fa"
            assert args[2] == "/path/to/genome"

    def test_run_calls_build_bowtie_index(self):
        """run() should call build_bowtie_index with correct arguments."""
        task = BuildBowtieIndex()
        with patch.object(task, "build_bowtie_index") as mock_build, \
             patch("eFISHent.util.log_stage_start"), \
             patch("eFISHent.util.get_genome_name", return_value="/tmp/test_genome"), \
             patch("eFISHent.indexing.GeneralConfig") as mock_config:
            mock_config.return_value.reference_genome = "/tmp/genome.fa"
            task.run()
            mock_build.assert_called_once()
            call_kwargs = mock_build.call_args[1]
            assert call_kwargs["genome"] == "/tmp/test_genome"


class TestBuildBowtie2IndexRun:
    """Test BuildBowtie2Index run method."""

    def test_run_calls_bowtie2_build(self):
        """run() should call bowtie2-build via subprocess."""
        task = BuildBowtie2Index()
        with patch("eFISHent.indexing.subprocess.check_call") as mock_call, \
             patch("eFISHent.util.log_stage_start"), \
             patch("eFISHent.util.get_genome_name", return_value="/tmp/test_genome"), \
             patch("eFISHent.indexing.GeneralConfig") as mock_config, \
             patch("eFISHent.console.spinner", return_value=MagicMock(__enter__=MagicMock(), __exit__=MagicMock())):
            mock_config.return_value.reference_genome = "/tmp/genome.fa"
            task.run()
            mock_call.assert_called_once()
            args = mock_call.call_args[0][0]
            assert args[0] == "bowtie2-build"
            assert "/tmp/genome.fa" in args[1]
            assert args[2] == "/tmp/test_genome"


class TestPrepareAnnotationFileRun:
    """Test PrepareAnnotationFile run method."""

    def test_run_calls_prepare_gtf_file(self):
        """run() should call prepare_gtf_file with correct paths."""
        task = PrepareAnnotationFile()
        with patch.object(task, "prepare_gtf_file") as mock_prep, \
             patch("eFISHent.indexing.GeneralConfig") as mock_config, \
             patch.object(task, "output") as mock_output:
            mock_config.return_value.reference_annotation = "/data/genome.gtf"
            mock_target = MagicMock()
            mock_target.path = "/data/genome.gtf.parq"
            mock_output.return_value = mock_target
            task.run()
            mock_prep.assert_called_once_with("/data/genome.gtf", "/data/genome.gtf.parq")
