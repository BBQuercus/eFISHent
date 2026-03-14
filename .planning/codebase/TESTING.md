# Testing Patterns

**Analysis Date:** 2026-03-14

## Test Framework

**Runner:**
- pytest (>=7.0)
- Config: `pyproject.toml` under `[dependency-groups] dev`
- No `pytest.ini` or `setup.cfg` pytest configuration file detected

**Assertion Library:**
- pytest assertions (standard `assert` statements)
- `pytest.approx()` for floating-point comparisons
- `pytest.raises()` for exception testing

**Run Commands:**
```bash
pytest                           # Run all tests
pytest tests/test_unit_*.py      # Run unit tests only
pytest tests/test_integration.py # Run integration tests
pytest -v                        # Verbose output
pytest -m skipif                 # Run tests with external tool requirements
```

## Test File Organization

**Location:**
- Co-located in `tests/` directory (separate from `eFISHent/`)
- Maintains parallel module structure for clarity

**Naming:**
- Unit tests: `test_unit_<module>.py` (e.g., `test_unit_cli.py`, `test_unit_optimization.py`)
- Feature/integration tests: `test_<feature>.py` (e.g., `test_adaptive_probes.py`, `test_integration.py`)
- Data files: `tests/data/` directory

**Structure:**
```
tests/
├── conftest.py                    # Shared fixtures (mostly empty/commented)
├── test_unit_cli.py
├── test_unit_optimization.py
├── test_unit_basic_filtering.py
├── test_unit_kmers.py
├── test_unit_presets.py
├── test_unit_cleanup.py
├── test_integration.py
├── test_quality_scores.py
├── test_adaptive_probes.py
├── test_gene_annotation.py
└── data/
    ├── sacCer3.fa
    ├── sacCer3_15.jf              # Jellyfish k-mer index
    ├── renilla.fasta
    └── renilla_data_table.fasta
```

## Test Structure

**Suite Organization:**
```python
class TestStringToBool:
    """Tests for string_to_bool validator."""

    @pytest.mark.parametrize(
        "value", ["true", "True", "TRUE", "yes", "Yes", "y", "1", "t"]
    )
    def test_truthy_values(self, value):
        assert string_to_bool(value) is True

    @pytest.mark.parametrize(
        "value", ["false", "False", "FALSE", "no", "No", "n", "0", "f"]
    )
    def test_falsy_values(self, value):
        assert string_to_bool(value) is False
```

**Patterns:**
- Test classes group related test functions (one class per tested unit)
- Test function names start with `test_` followed by what is being tested: `test_truthy_values()`, `test_gc_content()`
- Docstrings on classes explain test purpose
- `@pytest.mark.parametrize()` heavily used for testing multiple inputs/outputs

**Setup/Teardown:**
- Fixtures defined with `@pytest.fixture` decorator
- `autouse=True` for fixtures that run before every test:
  ```python
  @pytest.fixture(autouse=True)
  def _reset_state():
      """Reset console module state between tests."""
      reset_funnel_data()
      set_silent_mode(False)
      yield
      set_silent_mode(False)
      reset_funnel_data()
  ```
- Context managers for temporary resources (`tempfile.TemporaryDirectory()`, `tempfile.NamedTemporaryFile()`)
- Manual cleanup in `finally` blocks when needed

**Assertion Pattern:**
- Direct assertions: `assert result == expected`
- Equality assertions: `assert get_gc_content(seq) == output`
- Comparative assertions: `assert get_melting_temp(seq, 100, 0) > get_melting_temp(seq, 100, 10)`
- Collection assertions: `assert col in df.columns`, `assert len(output) == len(df)`

## Mocking

**Framework:**
- Standard library `unittest.mock` via pytest
- Mock objects via inline class definitions in test parameters

**Patterns:**
```python
# Inline Config mock for parameter validation
class Config(luigi.Config):
    min_tm = luigi.FloatParameter(10)
    max_tm = luigi.FloatParameter(20)
    min_gc = luigi.FloatParameter(0)
    max_gc = luigi.FloatParameter(100)
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
```

**What to Mock:**
- Config objects - create lightweight Config subclasses for testing different parameter combinations
- External tools - use `shutil.which()` to check availability, skip tests if missing

**What NOT to Mock:**
- Core business logic (filters, validators, optimization functions)
- Data structures (pandas DataFrames, Bio.SeqRecord objects)
- Actual computations (melting temp, GC content)

## Fixtures and Factories

**Test Data:**
```python
@pytest.fixture
def df():
    return pd.DataFrame(
        {
            "name": ["seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"],
            "sequence": [
                "GTAATTACAAAATAAGCAACG",
                "GCTTGCTTTGAGATTTTGTTC",
                "TCAATTCTACTGTCTCAGT",
                "GTAATTACAAAATAAGCAACG",
                "GGGACGAATTCTTTGTCTATTC",
                "TGGATTTCATAATGTTTATTTCAC",
                "GGGATCTTACGACATAAATCG",
            ],
            "start": [0, 2, 4, 5, 8, 10, 20],
            "end": [3, 5, 7, 9, 11, 14, 25],
        }
    )

@pytest.fixture
def sequential_solution():
    return ["seq1", "seq3", "seq5", "seq7"]
```

**Location:**
- Fixtures defined at top of test files
- Shared fixtures would go in `conftest.py` (currently mostly empty)
- Test data files in `tests/data/` directory

## Coverage

**Requirements:**
- No coverage enforcement configured (no coverage thresholds in `pyproject.toml`)
- Not enforced but recommended practice

**View Coverage:**
```bash
pytest --cov=eFISHent --cov-report=html    # Generate HTML coverage report
pytest --cov=eFISHent --cov-report=term    # Print coverage to terminal
```

## Test Types

**Unit Tests:**
- Scope: Individual functions and methods in isolation
- Approach: Parametrized tests with multiple input/output pairs
- Examples: `test_unit_cli.py` tests validators, `test_unit_basic_filtering.py` tests filtering functions
- Pattern: Fast, no external dependencies (unless marked with skipif)
- File pattern: `test_unit_<module>.py`

**Integration Tests:**
- Scope: Full pipeline execution with external tools (bowtie2, jellyfish, BLAST)
- Approach: Conditional skips if tools not available
- Examples: `test_integration.py` runs full eFISHent workflow with actual tools
- Conditional execution:
  ```python
  BOWTIE2_AVAILABLE = shutil.which("bowtie2") is not None
  ALL_TOOLS_AVAILABLE = all([
      BOWTIE2_AVAILABLE,
      BOWTIE2_BUILD_AVAILABLE,
      JELLYFISH_AVAILABLE,
  ])

  @pytest.mark.skipif(not ALL_TOOLS_AVAILABLE, reason="Required tools not installed")
  def test_full_pipeline(tmp_path):
      # Full integration test
  ```

**E2E Tests:**
- Not explicitly used
- Integration tests serve as end-to-end validation

## Common Patterns

**Async Testing:**
- Not applicable (non-async codebase)

**Error Testing:**
```python
@pytest.mark.parametrize("value", ["abc", "1.5", "", "one"])
def test_non_integer(self, value):
    with pytest.raises(argparse.ArgumentTypeError) as exc_info:
        positive_int(value)
    assert "Invalid integer" in str(exc_info.value)
```

**Parameter Validation Testing:**
- Heavy use of `@pytest.mark.parametrize()` for testing validator functions
- Test matrix approach: multiple values, multiple expected outcomes
- Example from `test_unit_cli.py`:
  ```python
  @pytest.mark.parametrize("value,expected", [("1", 1), ("5", 5), ("100", 100)])
  def test_valid_values(self, value, expected):
      assert positive_int(value) == expected
  ```

**External Tool Dependency Handling:**
- Check tool availability with `shutil.which()`
- Mark tests with `@pytest.mark.skipif()` if tool not available
- Examples:
  - `@pytest.mark.skipif(not JELLYFISH_AVAILABLE, reason="jellyfish not installed")`
  - `@pytest.mark.skipif(not GLPK_AVAILABLE, reason="GLPK solver (glpsol) not installed")`
- Integration tests check multiple tools and skip gracefully

**Fixture Cleanup:**
- Manual file deletion in `finally` blocks:
  ```python
  try:
      output = task_cleanup.prettify_table(...)
      assert output is not None
  finally:
      os.unlink(alignment_path)
  ```
- `tempfile` context managers for automatic cleanup:
  ```python
  with tempfile.TemporaryDirectory() as tmp_dir:
      # Test code using tmp_dir
      # Auto-cleaned up on exit
  ```
- `pytest_sessionstart()` hook to clean up leftover files from previous runs

**Bio.SeqRecord Testing:**
- Tests create SeqRecord objects inline:
  ```python
  Bio.SeqRecord.SeqRecord(Bio.Seq.Seq("ATGCATGC"), id="sequence")
  ```
- Used to test sequence processing functions
- Real FASTA files in `tests/data/` for integration tests

---

*Testing analysis: 2026-03-14*
