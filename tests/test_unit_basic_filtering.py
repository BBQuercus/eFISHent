import Bio.Seq
import Bio.SeqRecord
import luigi
import pytest

from eFISHent.basic_filtering import BasicFiltering
from eFISHent.basic_filtering import get_g_quadruplet_count
from eFISHent.basic_filtering import get_gc_content
from eFISHent.basic_filtering import get_melting_temp


@pytest.mark.parametrize(
    "seq",
    [
        Bio.Seq.Seq("ATGCATGC"),
        Bio.Seq.Seq("GCGCGCGC"),
        Bio.Seq.Seq("ATATATATTATAATATAT"),
    ],
)
def test_melting_temperature(seq):
    # More formamide -> lower melting temp
    assert (
        get_melting_temp(seq, 100, 0)
        > get_melting_temp(seq, 100, 10)
        > get_melting_temp(seq, 100, 20)
    )

    # More salt -> higher melting temp
    assert (
        get_melting_temp(seq, 10, 0)
        < get_melting_temp(seq, 100, 0)
        < get_melting_temp(seq, 200, 0)
    )


@pytest.mark.parametrize(
    "seq,output",
    [
        (Bio.Seq.Seq("ATATAT"), 0),
        (Bio.Seq.Seq("GCGCGC"), 100),
        (Bio.Seq.Seq("AGTCAGTC"), 50),
    ],
)
def test_gc_content(seq, output):
    assert get_gc_content(seq) == output


@pytest.mark.parametrize(
    "seq,output",
    [
        (Bio.Seq.Seq("ATTTTTC"), 0),
        (Bio.Seq.Seq("AGGGGGGA"), 1),
        (Bio.Seq.Seq("AGGGGGGA"), 1),
        (Bio.Seq.Seq("AGGGGCGGGGT"), 2),
    ],
)
def test_quadruplet_count(seq, output):
    assert get_g_quadruplet_count(seq) == output


@pytest.mark.parametrize(
    "seq,tm,gc,valid",
    [
        ("AGAGAGAGA", (10, 20), (0, 100), True),
        ("AGAGAGGGA", (10, 20), (0, 100), False),
        ("ATCGATGCAC", (20, 30), (40, 60), True),
        ("ATCGATGCAC", (0, 10), (0, 100), False),
        ("ATCGATGCAC", (20, 30), (60, 70), False),
    ],
)
def test_is_valid(seq, tm, gc, valid):
    class Config(luigi.Config):
        min_tm = luigi.FloatParameter(tm[0])
        max_tm = luigi.FloatParameter(tm[1])
        min_gc = luigi.FloatParameter(gc[0])
        max_gc = luigi.FloatParameter(gc[1])
        na_concentration = luigi.IntParameter(390)
        formamide_concentration = luigi.IntParameter(10)

    assert (
        BasicFiltering().is_candidate_valid(
            Bio.SeqRecord.SeqRecord(seq, id="sequence"), Config()
        )
        == valid
    )
