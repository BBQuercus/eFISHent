from unittest.mock import MagicMock, patch

import luigi
import pandas as pd
import pytest

from eFISHent.indexing import (
    BuildBowtie2Index,
    BuildBowtieIndex,
    PrepareAnnotationFile,
    read_annotation,
)


@pytest.fixture
def task_prepare():
    return PrepareAnnotationFile()


def test_prepare_annotations(task_prepare, tmp_path):
    fname_output = str(tmp_path / "sacCer3.gtf.parq")
    task_prepare.prepare_gtf_file("./tests/data/sacCer3.gtf", fname_output)

    df = pd.read_parquet(fname_output)
    for col in ["seqname", "feature", "start", "end", "frame", "exon_number"]:
        assert col in df.columns
    assert set(df["feature"]) == {"gene", "transcript", "exon"}
    assert len(df) == 7127 + 7127 + 7507
    assert df["start"].dtype == "int64"
    assert isinstance(df["seqname"].dtype, pd.CategoricalDtype)
    # Only attributes downstream code reads; gene_type is absent from this GTF
    assert "gene_biotype" in df.columns
    assert "gene_type" not in df.columns
    assert "gene_source" not in df.columns
    assert "source" not in df.columns


GTF_LINES = [
    "#!genome-build test",
    'chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "ONE"; gene_type "rRNA";',
    'chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ONE";',
    'chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";',
    'chr1\tsrc\tCDS\t120\t200\t.\t+\t0\tgene_id "G1"; transcript_id "T1";',
    'chr2\tsrc\tgene\t10\t50\t.\t-\t.\tgene_id "G2"; gene_type "lncRNA";',
    "malformed line",
]


def _write_gtf(tmp_path):
    path = tmp_path / "small.gtf"
    path.write_text("\n".join(GTF_LINES) + "\n")
    return str(path)


def test_prepare_annotations_keeps_features_and_attributes(task_prepare, tmp_path):
    out = str(tmp_path / "small.gtf.parq")
    task_prepare.prepare_gtf_file(_write_gtf(tmp_path), out)
    df = pd.read_parquet(out)

    assert list(df["feature"]) == ["gene", "transcript", "exon", "gene"]
    assert list(df["start"]) == [100, 100, 100, 10]
    assert df["gene_type"].tolist()[::3] == ["rRNA", "lncRNA"]
    assert df["gene_type"].isna().tolist() == [False, True, True, False]
    assert df["gene_name"].isna().tolist() == [False, False, True, True]
    assert df.loc[2, "exon_number"] == "1"
    # gene_biotype never appears, so the column is dropped like a full parse would
    assert "gene_biotype" not in df.columns


def test_prepare_annotations_chunking_is_transparent(task_prepare, tmp_path):
    gtf = _write_gtf(tmp_path)
    whole = str(tmp_path / "whole.parq")
    chunked = str(tmp_path / "chunked.parq")
    task_prepare.prepare_gtf_file(gtf, whole)
    task_prepare.prepare_gtf_file(gtf, chunked, chunk_size=1)
    pd.testing.assert_frame_equal(pd.read_parquet(whole), pd.read_parquet(chunked))


def test_prepare_annotations_empty_gtf(task_prepare, tmp_path):
    gtf = tmp_path / "empty.gtf"
    gtf.write_text("#!only a header\n")
    out = str(tmp_path / "empty.parq")
    task_prepare.prepare_gtf_file(str(gtf), out)
    df = pd.read_parquet(out)
    assert df.empty
    assert "feature" in df.columns


def test_read_annotation_skips_missing_columns(tmp_path):
    path = str(tmp_path / "a.parq")
    pd.DataFrame({"feature": ["gene"], "gene_biotype": ["rRNA"], "extra": [1]}).to_parquet(path)

    df = read_annotation(path, ["feature", "gene_biotype", "gene_type"])
    assert list(df.columns) == ["feature", "gene_biotype"]
    assert list(read_annotation(path).columns) == ["feature", "gene_biotype", "extra"]


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
