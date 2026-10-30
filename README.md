# IBD microbiome transferability benchmark

This repository contains the minimal, data-free implementation package for the
IBD cross-cohort transferability benchmark. It implements the frozen analytical
contract used by the manuscript without shipping participant-level data,
manuscript files, derived result tables, or local workstation paths.

## Scope

The primary estimand is a fixed ordered 126-species HMP2–PRISM transductive
contract, source-only preprocessing and model selection, six directed
source-to-target edges, ten source-side perturbations, and participant-cluster
bootstrap uncertainty. The package does not estimate strict source-specific
feature-universe deployment, clinical utility, real-world confounding control,
broad population generalization, or biological invariance.

The 20% source reserve is used for source perturbation and is not an outer test
set. Model selection is performed by grouped inner cross-validation within the
source training partition. Target labels are read only by scoring scripts.

## Repository map

```text
src/ibd_microbiome_framework/   reusable contract, harmonization, transforms,
                                splits, models, benchmark and synthetic controls
scripts/                        primary benchmark, audit and uncertainty entry points
configs/                        sanitized experiment configuration templates
tests/                          unit and contract tests
data/                           input-schema documentation; no raw data
docs/                           reproducibility and reconstruction notes
```

## Installation

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows: .venv\\Scripts\\activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

The standard model benchmark additionally requires `xgboost`. Nielsen export
requires the optional `third-cohort` dependencies. The package intentionally
does not bundle those public cohort files; obtain them under their original
licenses and record their SHA-256 values in the run manifest.

## Input boundary

Commands accept explicit paths to user-supplied matrices and manifests. Use the
schemas in `data/input_schema.md` and the configuration templates under
`configs/`. Never place participant-level inputs or access tokens in this
repository. The source feature contract must be generated and hashed before a
locked target is projected.

## Primary commands

```bash
python scripts/build_data_contract.py \
  --legacy-root /path/to/frozen_inputs \
  --output-dir results/contract

python scripts/run_standard_model_benchmark.py \
  --contract-dir results/contract \
  --legacy-root /path/to/frozen_inputs \
  --nielsen-projection /path/to/nielsen_projection.csv \
  --nielsen-manifest /path/to/nielsen_manifest.csv \
  --output-dir results/standard_benchmark

python scripts/run_reciprocal_nielsen_standard_benchmark.py \
  --contract-dir results/contract \
  --legacy-root /path/to/frozen_inputs \
  --nielsen-projection /path/to/nielsen_projection.csv \
  --nielsen-manifest /path/to/nielsen_manifest.csv \
  --output-dir results/reciprocal_nielsen

python scripts/bootstrap_grouped_predictions.py \
  --predictions results/standard_benchmark/standard_model_benchmark_predictions.csv \
  --output-dir results/cluster_bootstrap \
  --group-columns direction,model,split \
  --aggregate-columns direction,model
```

The commands write machine-readable metrics, predictions, tuning records and
audit summaries. Figures and manuscript documents are built outside this code
package from a frozen result registry.

The SIAMCAT input/export and native-reference runner are optional and are not
merged into the CLR model ranking. They require a separately versioned R/SIAMCAT
environment and should be reported as a protocol-level reference only.

## Tests and checks

```bash
python -m compileall -q src scripts tests
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -p no:cacheprovider -q
ruff check src scripts tests
```

## Reproducibility and disclosure

Every public run should archive the exact input checksums, package revision,
configuration, software versions, seeds, contract hash and output manifest.
The release does not include the source manuscript, author metadata, funding
statements, code-availability URLs, or DOI archive identifiers.

The default license is intentionally not assumed. Replace
`LICENSE_PLACEHOLDER.txt` with the authors' approved license before public
release.
