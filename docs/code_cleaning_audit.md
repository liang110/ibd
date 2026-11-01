# Code cleaning audit

The public package is a copy of the local scaffold's reusable package and
primary runners. It was cleaned without changing the mathematical estimand or
the data-processing order.

- Removed redundant inline comments and retained behavior-defining docstrings.
- Removed unused imports and duplicate module docstrings.
- Applied Ruff formatting and Ruff lint checks.
- Removed generated `__pycache__`, `.pytest_cache` and `.ruff_cache` directories.
- Replaced the original absolute data catalogue with a relative-path example.
- Corrected the Nielsen cohort configuration to retain both Denmark and Spain.
- Kept all data access explicit through command-line paths; no input is bundled.

Validation performed on the package:

- AST parse: 49 Python files passed.
- Unit and contract tests: 17 passed with `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1`.
- Ruff check: passed for `src`, `scripts` and `tests`.
- Token-level comment scan: only executable shebangs remain.

The local scaffold itself was not modified by this cleaning pass.
