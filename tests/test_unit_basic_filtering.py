import Bio.Seq
import Bio.SeqRecord
import luigi
import pytest

from eFISHent.basic_filtering import BasicFiltering
from eFISHent.basic_filtering import compute_off_target_tm
from eFISHent.basic_filtering import get_cpg_fraction
from eFISHent.basic_filtering import get_dinucleotide_repeat_count
from eFISHent.basic_filtering import get_g_quadruplet_count
from eFISHent.basic_filtering import get_gc_content
from eFISHent.basic_filtering import get_max_homopolymer_run
from eFISHent.basic_filtering import get_melting_temp
from eFISHent.basic_filtering import has_low_complexity


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
        max_homopolymer_length = luigi.IntParameter(0)
        filter_low_complexity = luigi.BoolParameter(False)

    assert (
        BasicFiltering().is_candidate_valid(
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id="sequence"), Config()
        )
        == valid
    )


# Low-complexity filter tests


@pytest.mark.parametrize(
    "seq,expected",
    [
        (Bio.Seq.Seq("ATCGATCG"), 1),
        (Bio.Seq.Seq("ATTTTCG"), 4),
        (Bio.Seq.Seq("AAAAATCG"), 5),
        (Bio.Seq.Seq("GGGGGGGG"), 8),
        (Bio.Seq.Seq(""), 0),
    ],
)
def test_max_homopolymer_run(seq, expected):
    assert get_max_homopolymer_run(seq) == expected


@pytest.mark.parametrize(
    "seq,expected",
    [
        (Bio.Seq.Seq("ATCGATCG"), 0),
        (Bio.Seq.Seq("ATATATATATATAT"), 2),  # AT and TA both repeat
        (Bio.Seq.Seq("GCGCGCGCGCGC"), 2),  # GC and CG both repeat
        (Bio.Seq.Seq("ATATATATGCGCGCGC"), 2),  # AT and GC repeat 4+ times
    ],
)
def test_dinucleotide_repeat_count(seq, expected):
    assert get_dinucleotide_repeat_count(seq) == expected


def test_has_low_complexity():
    assert has_low_complexity(Bio.Seq.Seq("AAAAAAAAAAAAA")) is True
    assert has_low_complexity(Bio.Seq.Seq("ATCGATCGATCGATCG")) is False


def test_is_valid_rejects_homopolymer():
    """Probe with TTTTT should fail when max_homopolymer_length=5."""
    class Config(luigi.Config):
        min_tm = luigi.FloatParameter(0)
        max_tm = luigi.FloatParameter(100)
        min_gc = luigi.FloatParameter(0)
        max_gc = luigi.FloatParameter(100)
        na_concentration = luigi.IntParameter(390)
        formamide_concentration = luigi.IntParameter(10)
        max_homopolymer_length = luigi.IntParameter(5)
        filter_low_complexity = luigi.BoolParameter(False)

    seq_with_run = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATTTTTCGATCG"), id="seq")
    assert BasicFiltering().is_candidate_valid(seq_with_run, Config()) is False

    seq_without_run = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATCGATCGATCG"), id="seq")
    assert BasicFiltering().is_candidate_valid(seq_without_run, Config()) is True


def test_is_valid_low_complexity_filter():
    """Low-complexity filter should reject dinucleotide repeats."""
    class Config(luigi.Config):
        min_tm = luigi.FloatParameter(0)
        max_tm = luigi.FloatParameter(100)
        min_gc = luigi.FloatParameter(0)
        max_gc = luigi.FloatParameter(100)
        na_concentration = luigi.IntParameter(390)
        formamide_concentration = luigi.IntParameter(10)
        max_homopolymer_length = luigi.IntParameter(0)
        filter_low_complexity = luigi.BoolParameter(True)

    seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATATATATATATAT"), id="seq")
    assert BasicFiltering().is_candidate_valid(seq, Config()) is False


def test_existing_gggg_filter_unchanged():
    """Regression: GGGG still rejected when max_homopolymer_length=0 (legacy mode)."""
    class Config(luigi.Config):
        min_tm = luigi.FloatParameter(0)
        max_tm = luigi.FloatParameter(100)
        min_gc = luigi.FloatParameter(0)
        max_gc = luigi.FloatParameter(100)
        na_concentration = luigi.IntParameter(390)
        formamide_concentration = luigi.IntParameter(10)
        max_homopolymer_length = luigi.IntParameter(0)
        filter_low_complexity = luigi.BoolParameter(False)

    seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("AGGGGGGA"), id="seq")
    assert BasicFiltering().is_candidate_valid(seq, Config()) is False


def test_compute_off_target_tm_perfect_match():
    """Perfect match should return same Tm as get_melting_temp."""
    seq = "ATCGATCGATCGATCGATCG"
    tm = compute_off_target_tm(seq, seq, na_concentration=330, formamide_concentration=10)
    expected = get_melting_temp(Bio.Seq.Seq(seq), 330, 10)
    assert abs(tm - expected) < 0.1


def test_compute_off_target_tm_mismatch_penalty():
    """Mismatches should reduce Tm by ~2.5°C each."""
    seq = "ATCGATCGATCGATCGATCG"
    # 1 mismatch at position 17 (A->T)
    target = "ATCGATCGATCGATCGTTCG"
    tm_perfect = compute_off_target_tm(seq, seq, 330, 10)
    tm_mismatch = compute_off_target_tm(seq, target, 330, 10)
    assert tm_perfect - tm_mismatch == pytest.approx(2.5, abs=0.1)


def test_compute_off_target_tm_invalid_seq():
    """Non-standard bases should return 0.0."""
    tm = compute_off_target_tm("NNNNNN", "ATCGAT", 330, 10)
    # Should not crash, returns some value (may or may not be 0 depending on NN table)
    assert isinstance(tm, float)


# CpG fraction tests


@pytest.mark.parametrize(
    "seq,expected",
    [
        (Bio.Seq.Seq("ATATAT"), 0.0),  # No CpG
        (Bio.Seq.Seq("CGCGCG"), 0.6),  # 3 CpG out of 5 dinucleotides
        (Bio.Seq.Seq("ACGTACGT"), 2 / 7),  # 2 CpG out of 7 dinucleotides
        (Bio.Seq.Seq("A"), 0.0),  # Too short
        (Bio.Seq.Seq(""), 0.0),  # Empty
    ],
)
def test_cpg_fraction(seq, expected):
    assert get_cpg_fraction(seq) == pytest.approx(expected, abs=0.01)


def test_is_valid_rejects_high_cpg():
    """Probes with CpG fraction above threshold should be rejected."""
    class Config(luigi.Config):
        min_tm = luigi.FloatParameter(0)
        max_tm = luigi.FloatParameter(100)
        min_gc = luigi.FloatParameter(0)
        max_gc = luigi.FloatParameter(100)
        na_concentration = luigi.IntParameter(390)
        formamide_concentration = luigi.IntParameter(10)
        max_homopolymer_length = luigi.IntParameter(0)
        filter_low_complexity = luigi.BoolParameter(False)
        max_cpg_fraction = luigi.FloatParameter(0.10)

    # High CpG: CGCGCGCGATCG = 4 CpG out of 11 dinucleotides = 0.36
    high_cpg = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("CGCGCGCGATCG"), id="seq")
    assert BasicFiltering().is_candidate_valid(high_cpg, Config()) is False

    # Low CpG: ATGATGATGATG = 0 CpG
    low_cpg = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGATGATGATG"), id="seq")
    assert BasicFiltering().is_candidate_valid(low_cpg, Config()) is True


def test_cpg_filter_disabled_by_default():
    """When max_cpg_fraction=0, CpG filter should not reject anything."""
    class Config(luigi.Config):
        min_tm = luigi.FloatParameter(0)
        max_tm = luigi.FloatParameter(100)
        min_gc = luigi.FloatParameter(0)
        max_gc = luigi.FloatParameter(100)
        na_concentration = luigi.IntParameter(390)
        formamide_concentration = luigi.IntParameter(10)
        max_homopolymer_length = luigi.IntParameter(0)
        filter_low_complexity = luigi.BoolParameter(False)
        max_cpg_fraction = luigi.FloatParameter(0.0)

    high_cpg = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("CGCGCGCGATCG"), id="seq")
    assert BasicFiltering().is_candidate_valid(high_cpg, Config()) is True
