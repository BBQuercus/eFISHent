#!/usr/bin/env bash
# Static analysis that CI runs in the "static-analysis" job of .github/workflows/tests.yml.
#
# - vulture: dead code (config in pyproject.toml, known false positives in vulture_whitelist.py)
# - deptry: unused, missing and transitive dependencies (config in pyproject.toml)
# - jscpd: duplicated code above the baseline threshold in .jscpd.json (needs Node/npx)
#
# Usage:
#   bash scripts/lint_dead.sh

set -euo pipefail

cd "$(dirname "$0")/.."

uv run --frozen vulture
uv run --frozen deptry .
npx --yes jscpd@5.4.3
