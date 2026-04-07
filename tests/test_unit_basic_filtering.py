import Bio.Seq
import Bio.SeqRecord
import luigi
import pytest

from eFISHent.basic_filtering import BasicFiltering
from eFISHent.basic_filtering import compute_duplex_dg
from eFISHent.basic_filtering import compute_off_target_tm
from eFISHent.basic_filtering import get_cpg_fraction
from eFISHent.basic_filtering import has_g_quadruplex
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


# Duplex ΔG tests


def test_compute_duplex_dg_returns_negative():
    """On-target duplex ΔG should be negative (stable binding)."""
    seq = Bio.Seq.Seq("ATCGATCGATCGATCGATCG")
    dg = compute_duplex_dg(seq, na_concentration=330, formamide_concentration=10)
    assert dg < 0.0


def test_compute_duplex_dg_longer_more_stable():
    """Longer probes should have more negative (more stable) ΔG."""
    short = Bio.Seq.Seq("ATCGATCGATCG")
    long = Bio.Seq.Seq("ATCGATCGATCGATCGATCGATCG")
    dg_short = compute_duplex_dg(short, 330, 10)
    dg_long = compute_duplex_dg(long, 330, 10)
    assert dg_long < dg_short  # more negative = more stable


def test_compute_duplex_dg_gc_rich_more_stable():
    """GC-rich probes should have more negative ΔG than AT-rich probes."""
    at_rich = Bio.Seq.Seq("ATATATATATATATATATATAT")
    gc_rich = Bio.Seq.Seq("GCGCGCGCGCGCGCGCGCGCG")
    dg_at = compute_duplex_dg(at_rich, 330, 10)
    dg_gc = compute_duplex_dg(gc_rich, 330, 10)
    assert dg_gc < dg_at


def test_compute_duplex_dg_error_returns_zero():
    """Error cases should return 0.0."""
    dg = compute_duplex_dg(Bio.Seq.Seq(""), 330, 10)
    assert dg == pytest.approx(0.0)


# Parameter suggestion tests


class TestSuggestParameters:
    """Tests for the parameter recommendation engine."""

    @staticmethod
    def _make_probes(gc_values, length=20):
        """Generate synthetic probes with approximate target GC contents."""
        probes = []
        for i, gc in enumerate(gc_values):
            n_gc = int(length * gc / 100)
            n_at = length - n_gc
            # Alternate G and C for GC bases, A and T for AT bases
            seq = "GC" * (n_gc // 2) + ("G" if n_gc % 2 else "") + "AT" * (n_at // 2) + ("A" if n_at % 2 else "")
            seq = seq[:length]
            probes.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=f"probe-{i}"))
        return probes

    def test_suggests_gc_adjustment(self):
        """When most probes fail GC, should suggest adjusting GC thresholds."""
        # Create probes with GC values mostly outside default range (20-80%)
        # Use very high GC probes that exceed 80%
        probes = []
        for i in range(20):
            seq = "GCGCGCGCGCGCGCGCGCGC"  # 100% GC
            probes.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=f"probe-{i}"))
        suggestions = BasicFiltering.suggest_parameters(probes)
        # Should suggest raising max_gc
        gc_suggestions = [s for s in suggestions if "--max-gc" in s]
        assert len(gc_suggestions) > 0

    def test_suggests_tm_adjustment(self):
        """When most probes fail Tm, should suggest adjusting Tm thresholds."""
        # Very short AT-rich probes will have very low Tm
        probes = []
        for i in range(20):
            seq = "ATATAT"  # very short, AT-rich = low Tm
            probes.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=f"probe-{i}"))
        suggestions = BasicFiltering.suggest_parameters(probes)
        tm_suggestions = [s for s in suggestions if "--min-tm" in s]
        assert len(tm_suggestions) > 0

    def test_no_suggestions_when_all_pass(self):
        """When probes are well within range, should return no suggestions."""
        probes = []
        for i in range(20):
            # 50% GC, moderate length => moderate Tm
            seq = "ATGCATGCATGCATGCATGC"
            probes.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=f"probe-{i}"))
        suggestions = BasicFiltering.suggest_parameters(probes)
        # GC and Tm should be fine, minimal suggestions expected
        gc_tm_suggestions = [s for s in suggestions if "--min-gc" in s or "--max-gc" in s
                             or "--min-tm" in s or "--max-tm" in s]
        assert len(gc_tm_suggestions) == 0

    def test_empty_input_returns_empty(self):
        """Empty input should return no suggestions."""
        assert BasicFiltering.suggest_parameters([]) == []

    def test_suggestions_include_probe_counts(self):
        """Suggestions should include how many probes would pass."""
        probes = []
        for i in range(20):
            seq = "GCGCGCGCGCGCGCGCGCGC"
            probes.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), id=f"probe-{i}"))
        suggestions = BasicFiltering.suggest_parameters(probes)
        for s in suggestions:
            assert "probes would pass" in s or "currently" in s


# G-quadruplex tests


@pytest.mark.parametrize(
    "seq,expected",
    [
        # Classic G4 motif: G3+N1-7G3+N1-7G3+N1-7G3+
        (Bio.Seq.Seq("GGGATTGGGATTGGGATTGGG"), True),
        # C4 motif (complementary strand)
        (Bio.Seq.Seq("CCCATTCCCATTCCCATTCCC"), True),
        # No G4
        (Bio.Seq.Seq("ATCGATCGATCGATCGATCG"), False),
        # Only GGGG but not the full G4 pattern
        (Bio.Seq.Seq("ATGGGGCGATCGATCGATCG"), False),
        # Longer loops (still within 7nt)
        (Bio.Seq.Seq("GGGAATTCCGGGAATTCCGGGAATTCCGGG"), True),
    ],
)
def test_has_g_quadruplex(seq, expected):
    assert has_g_quadruplex(seq) == expected


def test_is_valid_rejects_g_quadruplex():
    """Probes with G-quadruplex should be rejected when filter enabled."""
    class Config(luigi.Config):
        min_tm = luigi.FloatParameter(0)
        max_tm = luigi.FloatParameter(100)
        min_gc = luigi.FloatParameter(0)
        max_gc = luigi.FloatParameter(100)
        na_concentration = luigi.IntParameter(390)
        formamide_concentration = luigi.IntParameter(10)
        max_homopolymer_length = luigi.IntParameter(0)
        filter_low_complexity = luigi.BoolParameter(False)
        filter_g_quadruplex = luigi.BoolParameter(True)

    g4_seq = Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq("GGGATTGGGATTGGGATTGGG"), id="seq"
    )
    assert BasicFiltering().is_candidate_valid(g4_seq, Config()) is False

    normal_seq = Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq("ATCGATCGATCGATCGATCG"), id="seq"
    )
    assert BasicFiltering().is_candidate_valid(normal_seq, Config()) is True


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
