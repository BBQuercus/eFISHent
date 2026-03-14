# Coding Conventions

**Analysis Date:** 2026-03-14

## Naming Patterns

**Files:**
- Module files use snake_case: `basic_filtering.py`, `alignment.py`, `optimization.py`
- Test files follow pattern `test_unit_<module>.py` or `test_<feature>.py`: `test_unit_cli.py`, `test_integration.py`
- Package directory: `eFISHent/` (mixed case to match brand name)

**Functions:**
- Private helper functions use leading underscore: `_compute_gc()`, `_preferred_length()`, `_parse_parquet_gtf()`
- Public functions use snake_case: `get_melting_temp()`, `is_overlapping()`, `string_to_bool()`
- Validator functions follow pattern `<noun>()` or `<adjective>_<noun>()`: `positive_int()`, `existing_file()`, `existing_fasta_file()`

**Variables:**
- snake_case for local variables and module-level constants: `na_concentration`, `sequences`, `PIPELINE_STAGES`
- ALL_CAPS for module-level constants: `SAMFILE_COLUMNS`, `FASTA_EXT`, `SAM_FLAG_REVERSE`
- Config parameters match Luigi config attribute names: `min_length`, `max_length`, `spacing`

**Types:**
- PascalCase for classes: `GeneralConfig`, `ProbeConfig`, `AlignProbeCandidates`, `CleanUpOutput`
- Classes inherit from `luigi.Task`, `luigi.Config`, or standard Python types
- Type hints used throughout: `List[str]`, `Dict[str, dict]`, `Tuple[int, int]`, `Optional[Progress]`

## Code Style

**Formatting:**
- No explicit formatter configuration detected (no `.black`, `.flake8` files)
- Ruff linting enabled with per-file configuration
- Per-file ignores in `pyproject.toml`:
  - `eFISHent/cli.py`: E402 (module level import not at top), F401 (unused imports)

**Linting:**
- Ruff linter configured in `pyproject.toml`
- Tool: `ruff` (modern Python linter)
- No explicit `ruff.toml` configuration file

## Import Organization

**Order:**
1. Standard library imports (pathlib, typing, logging, os, etc.)
2. Third-party imports (Bio, luigi, pandas, matplotlib, numpy, etc.)
3. Local imports from current package (from . import util, from .config import GeneralConfig)
4. TYPE_CHECKING block for optional imports: `if TYPE_CHECKING: from .analyze import AnalyzeProbeset`

**Example from `cli.py`:**
```python
from pathlib import Path
from typing import Any, Dict, List, TYPE_CHECKING
import argparse
import configparser
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

import luigi

from . import __version__
from .constants import CLI_SHORTFORM
from .util import UniCode

if TYPE_CHECKING:
    from .analyze import AnalyzeProbeset
    from .cleanup import CleanUpOutput
```

**Path Aliases:**
- No path aliases (`jsconfig.json` or `tsconfig.json` equivalents)
- Relative imports use dot notation: `from .config import`, `from . import util`

## Error Handling

**Patterns:**
- Validators raise `argparse.ArgumentTypeError` with descriptive messages:
  ```python
  if not 0 <= fvalue <= 100:
      raise argparse.ArgumentTypeError(f"Percentage must be 0-100, got {fvalue}")
  ```
- Custom exceptions not used; relies on standard Python exceptions
- File existence and type validation in CLI validators: `existing_file()`, `existing_fasta_file()`
- Logger assertions via `log_and_check_candidates()` utility for filtering stage results
- `ValueError` raised for invalid configuration: `raise ValueError(f"Invalid optimization method: {RunConfig().optimization_method}")`

## Logging

**Framework:** `logging` (standard library)

**Patterns:**
- All classes that are Luigi tasks have class-level logger: `logger = logging.getLogger("custom-logger")`
- Logger name hardcoded to `"custom-logger"` across all task classes
- Logging configured via `set_logging_level()` in `cli.py`
- Rich logging handler via `get_rich_handler()` in `console.py`
- Two-tier output: Rich console for user display + log handler for file/debug output
- Suppressed loggers during CLI startup:
  ```python
  logging.getLogger("luigi-interface").setLevel(level=logging.CRITICAL)
  logging.getLogger("matplotlib.font_manager").setLevel(level=logging.CRITICAL)
  ```

## Comments

**When to Comment:**
- Module docstrings present on all files describing purpose
- Function docstrings describe complex logic or unclear intent
- Inline comments explain non-obvious algorithms (e.g., entropy calculation, mismatch penalty in `basic_filtering.py`)
- No trailing explanatory comments; docstrings preferred

**Docstring Style:**
- Google-style docstrings with description, then sections:
  ```python
  def get_free_energy(sequence: Bio.SeqRecord.SeqRecord) -> float:
      """Get free energy of sequence using RNAfold.

      Returns 0.0 on subprocess error (e.g., RNAfold not installed).
      """
  ```

**JSDoc/TSDoc:**
- Not applicable (Python project)

## Function Design

**Size:**
- Generally short, focused functions (15-40 lines typical)
- Longer methods in Luigi tasks (run methods 50-100+ lines) due to workflow orchestration
- Helper functions extracted for reusability (e.g., validation, filtering logic)

**Parameters:**
- Explicit parameter passing preferred over global state (except Config classes)
- Type hints required for function parameters and returns
- Config objects passed as parameters when needed (e.g., `config=Config()` in tests)
- Default parameter values used sparingly, mostly for optional arguments

**Return Values:**
- Functions return data structures or computed values, not None unless side-effect only
- Tuple returns for multiple values: `Tuple[int, int]`, `List[str]`
- Dict returns for structured data: `Dict[str, Dict]` (presets)
- Optional returns for potentially missing values: `Optional[Progress]`

## Module Design

**Exports:**
- Implicit exports (all public functions available via import)
- No `__all__` declarations observed
- Private functions use leading underscore to indicate internal use

**Barrel Files:**
- `__init__.py` contains version: `__version__ = "0.0.9"`
- No barrel exports in `__init__.py`

## Code Organization Patterns

**Luigi Task Hierarchy:**
- Tasks define `requires()`, `output()`, `run()` methods as standard workflow pattern
- `run()` method handles:
  1. Logging via `log_stage_start()`
  2. Reading input sequences with `Bio.SeqIO.parse()`
  3. Creating pandas DataFrames for intermediate data
  4. Writing output via `Bio.SeqIO.write()` or `df.to_csv()`

**Config Pattern:**
- Configuration is Luigi Config subclass with typed Parameters
- Four main config classes: `GeneralConfig`, `RunConfig`, `SequenceConfig`, `ProbeConfig`
- Config classes stored in `constants.CONFIG_CLASSES`
- Configs accessed via `ConfigClass()` to get singleton instance
- Defaults specified in Parameter definitions, not hardcoded

**Data Table Pattern:**
- Intermediate results stored in pandas DataFrames
- Utility function `create_data_table()` converts FASTA to DataFrame with columns: name, length, sequence, start, end
- Filtering and metadata added to DataFrame columns
- Final table prettified before output via `prettify_table()`

---

*Convention analysis: 2026-03-14*
