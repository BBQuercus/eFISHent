"""Tests for adaptive probe generation."""

import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.generate_probes import (
    _compute_gc,
    _preferred_length,
    create_candidate_probes_generator,
)


@pytest.fixture
def short_sequence():
    """A 50nt sequence with mixed GC content."""
    return Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"),
        id="test_seq",
    )


@pytest.fixture
def high_gc_sequence():
    """A 60nt sequence with high GC content (~70%)."""
    return Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATGCGCGCGCGCGCGCGCGCGCGCGCGCG"),
        id="high_gc",
    )


@pytest.fixture
def low_gc_sequence():
    """A 60nt sequence with low GC content (~30%)."""
    return Bio.SeqRecord.SeqRecord(
        Bio.Seq.Seq("ATATATATATATATATATATGCATATATATATATATATATATATATATATATATATATATATAT"),
        id="low_gc",
    )


class TestComputeGC:
    def test_all_gc(self):
        assert _compute_gc("GCGCGCGC") == 1.0

    def test_no_gc(self):
        assert _compute_gc("ATATATAT") == 0.0

    def test_half_gc(self):
        assert _compute_gc("ATGC") == pytest.approx(0.5)

    def test_empty_sequence(self):
        assert _compute_gc("") == 0.0

    def test_case_insensitive(self):
        assert _compute_gc("atgc") == pytest.approx(0.5)

    def test_known_sequence(self):
        # 3 G/C out of 10
        assert _compute_gc("AAAGCGATTT") == pytest.approx(0.3)


class TestPreferredLength:
    def test_high_gc_prefers_short(self):
        """High GC (>55%) should prefer minimum length."""
        result = _preferred_length(0.70, min_length=18, max_length=25)
        assert result == 18

    def test_low_gc_prefers_long(self):
        """Low GC (<45%) should prefer maximum length."""
        result = _preferred_length(0.30, min_length=18, max_length=25)
        assert result == 25

    def test_balanced_gc_prefers_middle(self):
        """Balanced GC (45-55%) should prefer middle length."""
        result = _preferred_length(0.50, min_length=18, max_length=25)
        expected_middle = 18 + (25 - 18) // 2  # 21
        assert result == expected_middle

    def test_narrow_range(self):
        """When range < 2, should return min_length."""
        result = _preferred_length(0.50, min_length=20, max_length=21)
        assert result == 20

    def test_same_min_max(self):
        """When min == max, should return min_length."""
        result = _preferred_length(0.50, min_length=20, max_length=20)
        assert result == 20

    def test_boundary_55(self):
        """GC exactly 0.55 should use balanced (middle) range."""
        result = _preferred_length(0.55, min_length=18, max_length=25)
        expected_middle = 18 + (25 - 18) // 2
        assert result == expected_middle

    def test_boundary_45(self):
        """GC exactly 0.45 should use balanced (middle) range."""
        result = _preferred_length(0.45, min_length=18, max_length=25)
        expected_middle = 18 + (25 - 18) // 2
        assert result == expected_middle


class TestAdaptiveProbes:
    def test_adaptive_fewer_probes(self, short_sequence):
        """Adaptive mode should produce fewer probes than standard mode."""
        min_len, max_len = 4, 8

        standard_probes = list(
            create_candidate_probes_generator(
                short_sequence, min_len, max_len, adaptive=False
            )
        )
        adaptive_probes = list(
            create_candidate_probes_generator(
                short_sequence, min_len, max_len, adaptive=True
            )
        )

        assert len(adaptive_probes) < len(standard_probes)

    def test_all_probes_within_bounds(self, short_sequence):
        """All adaptive probes should respect min/max length."""
        min_len, max_len = 4, 8
        probes = list(
            create_candidate_probes_generator(
                short_sequence, min_len, max_len, adaptive=True
            )
        )
        for probe in probes:
            assert min_len <= len(probe.seq) <= max_len

    def test_probe_ids_have_start_positions(self, short_sequence):
        """Probe IDs should contain correct start positions."""
        probes = list(
            create_candidate_probes_generator(
                short_sequence, 4, 8, adaptive=True
            )
        )
        for probe in probes:
            # ID format: candidate-{idx}-{start_pos}
            parts = probe.id.split("-")
            assert len(parts) == 3
            assert parts[0] == "candidate"
            start_pos = int(parts[2])
            assert 0 <= start_pos < len(short_sequence)

    def test_standard_mode_all_lengths(self, short_sequence):
        """Standard mode should produce probes at every length in range."""
        min_len, max_len = 4, 6
        probes = list(
            create_candidate_probes_generator(
                short_sequence, min_len, max_len, adaptive=False
            )
        )
        lengths = {len(p.seq) for p in probes}
        assert lengths == {4, 5, 6}

    def test_adaptive_with_same_min_max(self, short_sequence):
        """Adaptive mode with min==max should behave like standard."""
        standard = list(
            create_candidate_probes_generator(
                short_sequence, 6, 6, adaptive=False
            )
        )
        adaptive = list(
            create_candidate_probes_generator(
                short_sequence, 6, 6, adaptive=True
            )
        )
        # When min==max, adaptive can't adjust, should produce same count
        assert len(standard) == len(adaptive)

    def test_adaptive_high_gc_prefers_shorter(self, high_gc_sequence):
        """For high-GC sequences, adaptive should bias toward shorter probes."""
        min_len, max_len = 18, 25
        probes = list(
            create_candidate_probes_generator(
                high_gc_sequence, min_len, max_len, adaptive=True
            )
        )
        lengths = [len(p.seq) for p in probes]
        avg_length = sum(lengths) / len(lengths)
        midpoint = (min_len + max_len) / 2
        # Average should be below the midpoint for high-GC
        assert avg_length < midpoint

    def test_adaptive_low_gc_prefers_longer(self, low_gc_sequence):
        """For low-GC sequences, adaptive should bias toward longer probes."""
        min_len, max_len = 18, 25
        probes = list(
            create_candidate_probes_generator(
                low_gc_sequence, min_len, max_len, adaptive=True
            )
        )
        lengths = [len(p.seq) for p in probes]
        avg_length = sum(lengths) / len(lengths)
        midpoint = (min_len + max_len) / 2
        # Average should be above the midpoint for low-GC
        assert avg_length > midpoint

    def test_adaptive_error_handling(self, short_sequence):
        """Adaptive mode should still raise errors for invalid params."""
        with pytest.raises(ValueError):
            list(
                create_candidate_probes_generator(
                    short_sequence, 10, 5, adaptive=True
                )
            )
