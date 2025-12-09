import types

import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.generate_probes import GenerateAllProbes
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
