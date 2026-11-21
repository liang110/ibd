# IBD microbiome transferability benchmark

This code package supports the cross-cohort IBD post-fit evaluation protocol.
The source-side commands rebuild model predictions from permitted cohort inputs.
The registered-result replay scores the prediction registry with the saved
participant multiplicities distributed in the accompanying evidence bundle.

## Scope

The benchmark uses 126 ordered canonical taxonomic coordinates from the
HMP2–Franzosa feature intersection, six directed transfers and ten source-side
perturbations. Four reporting objects cross the scoring unit (sample or participant)
with the prediction rule (metrics averaged across fitted predictors or metric
computed once on averaged probabilities). The 500 shared target-participant
bootstrap draws quantify uncertainty conditional on the observed cohort and fits.

The 20% source reserve is used for source perturbation. Model selection uses
grouped inner cross-validation when participants contribute repeated observations.
The final fit uses all eligible source records. Target outcomes enter the evaluation.

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

The predictors use microbial abundances. Reconstructing the manuscript HMP2
panel also uses the provider's metabolomics specimen identifiers for eligibility.
The identifier intersection contained 388 collections from 105 participants.
Excluding three collections with zero abundance on the projected vocabulary
retained 385 collections. The source reconstruction guide `SOURCE_RECONSTRUCTION.md`
in the manuscript reader ZIP gives the sample-selection steps. The optional
`--include-metabolome` switch controls the additional source audit.

The source training commands below reproduce the manuscript shallow-MLP
configuration with `--allow-ungrouped-early-stopping`. The opt-in enables
scikit-learn's sample-level internal validation (`validation_fraction=0.15`,
`n_iter_no_change=30`). The public default disables this internal early stopping.
The candidate-selection folds remain participant-grouped for repeated observations.
Each training command records the flag in its run manifest.

## Input boundary

Analysis commands accept explicit paths to user-supplied cohort matrices, manifests,
and (where applicable) a Nielsen projection; the repository does not download or
embed participant-level inputs. Use the schemas in `data/input_schema.md` and the
configuration templates under `configs/`. Never place participant-level inputs or
access tokens in this repository. The source feature contract must be generated and
hashed before a locked target is projected.

## Data preparation tutorial

1. Obtain the HMP2/IBDMDB, the Franzosa study (PRISM discovery and Dutch validation) and NielsenHB_2014 source files through their
original provider accessions and permissions (PRJNA398089, PRJNA400072 and
PRJEB1220). Keep the downloaded files outside this repository and record their
SHA-256 values in a run manifest.
2. Arrange each cohort according to `data/input_schema.md`; provide a sample table,
a participant table, a binary outcome column and the cohort-specific feature matrix.
Use the explicit cohort label allow-lists described in the Input boundary section.
3. Put the permitted inputs under a local path such as
`/path/to/cohort_inputs/{hmp2,prism,nielsen}` and place the sanitized configuration
template under `configs/`. Do not change the ordered 126-coordinate contract after
`build_data_contract.py` writes `contract_manifest.json`.
4. Run the contract builder first, then the standard and reciprocal benchmark
commands below. Each command writes to a separate directory under `results/`; the
input contract, prediction table, shared participant weights and output manifest are
retained together.
5. To score the provided case study without refitting, use the commands under
   [Registered-result replay](#registered-result-replay).

## Registered-result replay

The source-side implementation for this manuscript is pinned at
[5885481c66f29356e4ca65ba801b64b5d3e2c30a](https://github.com/liang110/ibd/tree/5885481c66f29356e4ca65ba801b64b5d3e2c30a).
For scoring and figure sources, extract `IBD_PostFit_Reader_Evidence.zip` and
run the commands below from the directory containing `Source_Index.csv`.
Use Python 3.11 and the versions in `requirements-reader.txt`.

```shell
python -B verify_package.py
python -B protocol_control_20261109/check_example.py --output-dir ../reader_example
python -B analysis/declaration/score_declaration.py --declaration analysis/declaration/nielsen_example.json --output ../declared_scores.csv --contrasts-output ../declared_contrasts.csv
python -B analysis/import_case/run_import.py --bundle . --output ../reader_siamcat
python -B replay_all.py --output ../reader_full
```

The MLP example calculates twelve complete-target scores from ten runs on the
same 385 HMP2 specimens. The declaration example uses HMP2-trained logistic
regression on Nielsen to distinguish equal specimen weights, equal total weight
per person with specimen predictions retained, and one mean prediction per
person. The SIAMCAT adapter reads ten native exports. Full replay writes
`reader_validation.json`, including checks of 504 principal cell summaries
and 630 paired within-model contrasts. The ZIP's root README lists individual
inputs, expected outputs and participant multiplicities.

## Source reconstruction commands

```bash
python scripts/build_data_contract.py \
  --legacy-root /path/to/cohort_inputs \
  --output-dir results/contract

python scripts/run_standard_model_benchmark.py \
  --contract-dir results/contract \
  --legacy-root /path/to/cohort_inputs \
  --nielsen-projection /path/to/nielsen_projection.csv \
  --nielsen-manifest /path/to/nielsen_manifest.csv \
  --allow-ungrouped-early-stopping \
  --output-dir results/standard_benchmark

python scripts/run_reciprocal_nielsen_standard_benchmark.py \
  --contract-dir results/contract \
  --legacy-root /path/to/cohort_inputs \
  --nielsen-projection /path/to/nielsen_projection.csv \
  --nielsen-manifest /path/to/nielsen_manifest.csv \
  --allow-ungrouped-early-stopping \
  --output-dir results/reciprocal_nielsen

python scripts/bootstrap_grouped_predictions.py \
  --predictions results/standard_benchmark/standard_model_benchmark_predictions.csv \
  --output-dir results/cluster_bootstrap \
  --group-columns direction,model,split \
  --aggregate-columns direction
```

The bootstrap example shows how to generate shared participant weights for a
new source-side analysis. With `--aggregate-columns direction`, every model and
run in the same direction shares a draw. It writes
`shared_target_bootstrap_weights.csv`. For exact manuscript replay, use the
distributed weights and the registered replay command above.

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

## Manuscript and source data

Use the current reader ZIP's root `Source_Index.csv` to locate figures, tables
and Supplementary Data. It includes the publication item, input and output
paths, row counts, columns and SHA-256 values. Figure and table numbers used by
the source-fitting commit refer to the source-side analysis documentation.

| Manuscript item | Where to start |
|---|---|
| Main Table 2 and Figure 1 | Choose the four principal scoring rules or the participant-weighted specimen rule E. |
| Main Table 3 and Table S17 | Score the twelve worked-example metrics from the project MLP predictions. |
| Figures 2 and 3 | Check sampling denominators and the two changes along A to E to C. |
| Figure 4 | Compare conditional changes in the leading family. |
| Figure 5, Tables S20 and S23, SD24 | Inspect specimen eligibility and the twelve Nielsen source-by-family contrasts. |
| Tables S16 and S21 | Locate reader commands and check the native SIAMCAT example. |

## Tests and checks

```bash
python -m compileall -q src scripts tests
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -p no:cacheprovider -q
ruff check src scripts tests
```

## Reproducibility and disclosure

Every public run should archive the exact input checksums, package revision,
configuration, software versions, seeds, contract hash, canonicalizer hash, shared
bootstrap weight manifest and output manifest. Use `docs/runtime_manifest_template.json`
to record Python, platform, dependency versions, configuration hash, input hash,
output hash and bootstrap-weight hash rather than relying on lower-bound dependency
ranges. The source-reconstruction commands are `build_data_contract.py`,
`run_standard_model_benchmark.py` for four directions,
`run_reciprocal_nielsen_standard_benchmark.py` for two directions, and
`bootstrap_grouped_predictions.py --aggregate-columns direction`. For scoring from existing predictions and weights, see
[Registered-result replay](#registered-result-replay).
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


