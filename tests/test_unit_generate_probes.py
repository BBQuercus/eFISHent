import Bio.Seq
import Bio.SeqRecord
import pytest

from eFISHent.generate_probes import GenerateAllProbes


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
