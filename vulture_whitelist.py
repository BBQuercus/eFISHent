# Vulture whitelist: names vulture reports as unused that are in fact used, or kept on purpose.
# Not imported by anything. Run the dead-code check with `scripts/lint_dead.sh`.

# Library settings assigned for their side effect (rich-argparse, luigi).
_.group_name_formatter  # RichHelpFormatter class setting (eFISHent/arguments.py)
_.log_level  # luigi.interface.core setting (eFISHent/cli.py)

# Pyomo model components: assigning them to the model registers them with the solver.
_.objective  # eFISHent/selection.py
_.C1  # eFISHent/selection.py
_.C2  # eFISHent/selection.py

# Kept on purpose: backwards-compatible alias for callers outside this repo (eFISHent/util.py).
PIPELINE_STAGES

# Tests: pytest hook signature, and MagicMock magic-method setup.
pytest_sessionstart  # tests/test_integration.py – pytest only calls hooks from conftest.py; kept for review
session  # pytest_sessionstart argument
_.__getitem__  # MagicMock setup in tests/test_unit_alignment.py and tests/test_unit_cleanup.py
