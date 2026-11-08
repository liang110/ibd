# IBD microbiome transferability benchmark

This repository contains the minimal, data-free implementation package for the
IBD cross-cohort transferability benchmark. It implements the registered
analytical contract without shipping participant-level data, manuscript files,
derived result tables, or local workstation paths.

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

The primary microbiome contract is independent of metabolome files. Optional
metabolome matrices can be included in the contract audit only with
`--include-metabolome`; they are not used by the IBD transfer benchmark.

The public default disables scikit-learn MLP sample-level early stopping. The
manuscript comparator was run with `--allow-ungrouped-early-stopping`, which sets
`early_stopping=True`, `validation_fraction=0.15`, and `n_iter_no_change=30` for
the shallow MLP. This internal validation split is sample-level and is reported as a
descriptive comparator; it is not participant-aware validation. The flag must not be
used to support participant-aware validation claims, and its value is written to the
run manifest.

## Input boundary

Analysis commands accept explicit paths to user-supplied cohort matrices, manifests,
and (where applicable) a Nielsen projection; the repository does not download or
embed participant-level inputs. Use the schemas in `data/input_schema.md` and the
configuration templates under `configs/`. Never place participant-level inputs or
access tokens in this repository. The source feature contract must be generated and
hashed before a locked target is projected.

## Primary commands

```bash
python scripts/build_data_contract.py \
  --legacy-root /path/to/cohort_inputs \
  --output-dir results/contract

python scripts/run_standard_model_benchmark.py \
  --contract-dir results/contract \
  --legacy-root /path/to/cohort_inputs \
  --nielsen-projection /path/to/nielsen_projection.csv \
  --nielsen-manifest /path/to/nielsen_manifest.csv \
  --output-dir results/standard_benchmark

python scripts/run_reciprocal_nielsen_standard_benchmark.py \
  --contract-dir results/contract \
  --legacy-root /path/to/cohort_inputs \
  --nielsen-projection /path/to/nielsen_projection.csv \
  --nielsen-manifest /path/to/nielsen_manifest.csv \
  --output-dir results/reciprocal_nielsen

python scripts/bootstrap_grouped_predictions.py \
  --predictions results/standard_benchmark/standard_model_benchmark_predictions.csv \
  --output-dir results/cluster_bootstrap \
  --group-columns direction,model,split \
  --aggregate-columns direction
```

For the manuscript's common-target protocol, `aggregate-columns direction` is
required: it creates one participant-weight vector per direction and replicate,
then applies it to every source perturbation, split and model in that direction.
Using `aggregate-columns direction,model` is a different marginal-per-model
protocol and must not be used for cross-model paired contrasts. The bootstrap
command writes the vectors to `shared_target_bootstrap_weights.csv`. This is
the common-target weighting protocol used by the manuscript. The commands write
machine-readable metrics,
predictions, tuning records and audit summaries. Figures and manuscript
documents are built outside this code package from the recorded benchmark outputs.

The six standard comparator implementations and the prevalence baseline in
`src/ibd_microbiome_framework/models/standard.py` are the executable models for
the primary benchmark. They are not SIAMCAT reimplementations and must not be
removed merely because the SIAMCAT paper-analysis repository delegates model
fitting to the upstream SIAMCAT/mlr stack.

The SIAMCAT input/export and native-reference runner are optional adapters and
are not merged into the CLR model ranking. They require a separately versioned
R/SIAMCAT environment and an external reference script. The repository does
not contain the SIAMCAT R/Bioconductor implementation; obtain it from the
upstream SIAMCAT project and report the resulting analysis as a protocol-level
reference only. These adapters must not be described as a self-contained
SIAMCAT baseline or as an independent-data replication.

The primary common-target result path is the only path for manuscript Table 2,
Supplementary Tables 10, 16 and 17, and Figures 3–5. The Nielsen origin-composition
audit is Supplementary Figure S3. Marginal or prior
bootstrap scripts are retained for audit compatibility and must not be used to
regenerate the paired manuscript results.

## Tests and checks

```bash
python -m compileall -q src scripts tests
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -p no:cacheprovider -q
ruff check src scripts tests
```

## Reproducibility and disclosure

Every public run should archive the exact input checksums, package revision,
configuration, software versions, seeds, contract hash, shared bootstrap weight
manifest and output manifest. Use `docs/runtime_manifest_template.json` to
record the exact runtime rather than relying on lower-bound dependency ranges.
The code-only package does not include the source
manuscript, author metadata, funding statements, figures, derived result CSVs,
or DOI archive identifiers; those artifacts must be linked in the versioned
submission archive.

The default license is intentionally not assumed. Replace
`LICENSE_PLACEHOLDER.txt` with the authors' approved license before public
release.


## Explicit label-contract rules

HMP2 and PRISM use explicit CD/UC versus Control/nonIBD allow-lists; Nielsen uses IBD versus healthy allow-lists. Unknown, missing, conflicting, and cohort-incompatible labels fail closed. The shallow MLP is the manuscript comparator for post-fit analysis. Its reported run uses sample-level early stopping; the public default remains disabled because that internal validation split is not participant-grouped.

Manifest construction and post-fit participant aggregation both enforce the
participant-label invariant: one participant cannot map to more than one binary
outcome. The regression suite includes a conflicting-label fixture and fails
before any ``first()`` aggregation can occur.


