import types

import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.generate_probes import GenerateAllProbes
from eFISHent.generate_probes import _compute_gc
from eFISHent.generate_probes import _preferred_length
from eFISHent.generate_probes import create_candidate_probes_generator


@pytest.fixture
def task_generate():
    return GenerateAllProbes()


@pytest.fixture
def sequence():
    return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGC" * 4), id="my_id")


@pytest.mark.parametrize("subseq_length", [4, 6, 8])
def test_generate_probes_count(task_generate, sequence, subseq_length):
    output = task_generate.create_candidate_probes(
        sequence, min_length=subseq_length, max_length=subseq_length
    )
    expected_count = len(sequence) - subseq_length + 1
    assert len(output) == expected_count
    for candidate in output:
        assert len(candidate) == subseq_length


@pytest.mark.parametrize("min_length,max_length", [(2, 4), (1, 5), (5, 8)])
def test_generate_probes_length(task_generate, sequence, min_length, max_length):
    output = task_generate.create_candidate_probes(
        sequence, min_length=min_length, max_length=max_length
    )
    for candidate in output:
        assert min_length <= len(candidate) <= max_length


def test_generate_probes_error(task_generate, sequence):
    with pytest.raises(ValueError):
        task_generate.create_candidate_probes(sequence, min_length=5, max_length=3)

    with pytest.raises(ValueError):
        task_generate.create_candidate_probes(sequence, min_length=18, max_length=20)


class TestProbeGenerator:
    """Tests for the memory-efficient generator implementation."""

    def test_generator_returns_generator(self, sequence):
        """Test that create_candidate_probes_generator returns a generator."""
        result = create_candidate_probes_generator(sequence, 4, 6)
        assert isinstance(result, types.GeneratorType)

    def test_generator_produces_same_probes_as_list(self, task_generate, sequence):
        """Test that generator produces identical probes to list version."""
        list_output = task_generate.create_candidate_probes(sequence, 4, 6)
        generator_output = list(create_candidate_probes_generator(sequence, 4, 6))

        assert len(list_output) == len(generator_output)
        for list_probe, gen_probe in zip(list_output, generator_output):
            assert str(list_probe.seq) == str(gen_probe.seq)
            assert list_probe.id == gen_probe.id

    @pytest.mark.parametrize("subseq_length", [4, 6, 8])
    def test_generator_correct_count(self, sequence, subseq_length):
        """Test generator produces correct number of probes."""
        output = list(
            create_candidate_probes_generator(
                sequence, min_length=subseq_length, max_length=subseq_length
            )
        )
        expected_count = len(sequence) - subseq_length + 1
        assert len(output) == expected_count

    def test_generator_error_min_greater_than_max(self, sequence):
        """Test generator raises error when min_length > max_length."""
        with pytest.raises(ValueError):
            # Must consume the generator to trigger the error
            list(create_candidate_probes_generator(sequence, min_length=5, max_length=3))

    def test_generator_error_min_exceeds_sequence(self, sequence):
        """Test generator raises error when min_length >= sequence length."""
        with pytest.raises(ValueError):
            list(create_candidate_probes_generator(sequence, min_length=18, max_length=20))


# ── _compute_gc correctness ──────────────────────────────────────────────


class TestComputeGC:
    """Tests for the _compute_gc helper."""

    def test_pure_at(self):
        assert _compute_gc("AATTAATT") == 0.0

    def test_pure_gc(self):
        assert _compute_gc("GGCCGGCC") == 1.0

    def test_mixed(self):
        assert _compute_gc("ATGC") == 0.5

    def test_empty_string(self):
        assert _compute_gc("") == 0.0

    def test_case_insensitive(self):
        assert _compute_gc("atgc") == _compute_gc("ATGC")


# ── _preferred_length correctness ────────────────────────────────────────


class TestPreferredLength:
    """Tests for the _preferred_length helper."""

    def test_high_gc_returns_min(self):
        assert _preferred_length(0.60, 20, 30) == 20

    def test_low_gc_returns_max(self):
        assert _preferred_length(0.40, 20, 30) == 30

    def test_balanced_returns_middle(self):
        result = _preferred_length(0.50, 20, 30)
        assert result == 25  # 20 + 10 // 2

    def test_min_equals_max(self):
        assert _preferred_length(0.50, 25, 25) == 25

    def test_range_less_than_two(self):
        assert _preferred_length(0.50, 25, 26) == 25


# ── Generator vs list equivalence ────────────────────────────────────────


class TestGeneratorListEquivalence:
    """Verify that generator and list methods produce identical results."""

    @pytest.fixture
    def long_sequence(self):
        return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGC" * 10), id="long")

    def test_same_sequences(self, long_sequence):
        task = GenerateAllProbes()
        list_output = task.create_candidate_probes(long_sequence, 8, 12)
        gen_output = list(create_candidate_probes_generator(long_sequence, 8, 12))

        assert len(list_output) == len(gen_output)
        for lp, gp in zip(list_output, gen_output):
            assert str(lp.seq) == str(gp.seq)
            assert lp.id == gp.id

    def test_probe_ids_contain_position(self, long_sequence):
        probes = list(create_candidate_probes_generator(long_sequence, 8, 8))
        for probe in probes:
            # ID format: candidate-{idx}-{start_pos}
            parts = probe.id.split("-")
            assert len(parts) == 3
            assert parts[0] == "candidate"
            # start_pos should be a valid integer
            int(parts[2])


# ── Adaptive mode ────────────────────────────────────────────────────────


class TestAdaptiveMode:
    """Tests for the adaptive probe generation mode."""

    def test_high_gc_produces_shorter_probes(self):
        # Pure GC sequence -> adaptive should prefer shorter probes
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("GC" * 30), id="highgc")
        probes = list(create_candidate_probes_generator(seq, 10, 20, adaptive=True))
        lengths = [len(p) for p in probes]
        avg_length = sum(lengths) / len(lengths)
        # With high GC, preferred is min_length (10), so average should be near 10-11
        assert avg_length < 15

    def test_low_gc_produces_longer_probes(self):
        # Pure AT sequence -> adaptive should prefer longer probes
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("AT" * 30), id="lowgc")
        probes = list(create_candidate_probes_generator(seq, 10, 20, adaptive=True))
        lengths = [len(p) for p in probes]
        avg_length = sum(lengths) / len(lengths)
        # With low GC, preferred is max_length (20), so average should be near 19-20
        assert avg_length > 15

    def test_adaptive_vs_standard_probe_count(self):
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGC" * 15), id="mixed")
        standard = list(create_candidate_probes_generator(seq, 10, 14, adaptive=False))
        adaptive = list(create_candidate_probes_generator(seq, 10, 14, adaptive=True))
        # Adaptive generates fewer lengths per position (preferred ±1 vs all)
        assert len(adaptive) < len(standard)


# ── Edge cases ───────────────────────────────────────────────────────────


class TestEdgeCases:
    """Edge case tests for probe generation."""

    def test_sequence_one_longer_than_min_length(self):
        """When sequence is min_length + 1, exactly one probe of min_length fits."""
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGCATGCA"), id="exact")  # len=9
        probes = list(create_candidate_probes_generator(seq, 8, 10))
        # length 8: 2 probes, length 9: 1 probe, length 10: 0 probes
        assert len(probes) == 3
        eight_mers = [p for p in probes if len(p) == 8]
        assert len(eight_mers) == 2

    def test_sequence_exactly_min_length_raises(self):
        """When sequence length == min_length, a ValueError is raised."""
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGCATGC"), id="exact")  # len=8
        with pytest.raises(ValueError):
            list(create_candidate_probes_generator(seq, 8, 10))

    def test_single_base_min_length_one(self):
        seq = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("A"), id="single")
        # min_length must be < len(sequence), so min_length=1 with len=1 should fail
        with pytest.raises(ValueError):
            list(create_candidate_probes_generator(seq, 1, 1))
