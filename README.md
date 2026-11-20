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

The primary microbiome contract is independent of metabolome files. Optional
metabolome matrices can be included in the contract audit only with
`--include-metabolome`; they are not used by the IBD transfer benchmark.

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

Extract the manuscript reader distribution. From its top-level directory run
`python -B verify_package.py` to verify the file manifest. The evidence directory
inside it is named `submission_evidence_bundle_20261108`.

**S17 real-data example.** From inside the extracted evidence directory:

```shell
python -B protocol_control_20261109/check_example.py --output-dir ../reader_results/s17
```

The command checks the 3,850 PRISM→HMP2 shallow-MLP predictions against their
registry, computes the 12 full-target A–D scores with all participant multiplicities
equal to one, and compares the values to the archived outputs and independent
reference. Expected AUROCs for A–D are 0.635197098, 0.728818779, 0.653456670
and 0.749756573. The command writes `reader_check.json` and the invalid-input
check results in the selected output directory.

**Full weighted replay.** Run the following from this Git checkout in a
Python 3.11 environment. In the Bash example, set `EVIDENCE` to the extracted
evidence directory. The bundle contains `requirements-reader.txt` for the
tested environment. `scripts/replay_registered.py` and the corresponding
`validation/replay_registered.py` in the reader distribution implement the
same registered replay.

```bash
EVIDENCE=/path/to/extracted/submission_evidence_bundle_20261108
python -B scripts/replay_registered.py \
  --predictions "${EVIDENCE}/evidence/full_reproducibility_inputs/prediction_registry.csv" \
  --weights "${EVIDENCE}/evidence/full_reproducibility_inputs/shared_bootstrap_weights.csv" \
  --output-dir results/postfit
python -B "${EVIDENCE}/validation/check_replay_tables.py" --output-dir results/postfit
```

The evaluator consumes the archived participant multiplicities and writes
252,000 draw-level metric records, 504 cell summaries, 630 within-family paired
contrasts and 504 ranking entries. The optional comparison command matches the
four output tables to the archive by scientific keys at absolute tolerance
`1e-9`. The archive's top-level `README.md` also includes the independent
statistical audit entry.

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

| Manuscript item | Content |
|---|---|
| Figures 1–2 | Evaluation workflow and A-cell AUROC benchmark |
| Figure 3 | Paired changes in aggregation order and scoring unit |
| Figure 4 | Ranking impact and stability |
| Figure 5 | AUROC and Brier comparisons |
| Main Tables 1–3 | Cohort composition, A-cell benchmark, and pre-analysis choice of A–D |
| Supplementary Tables S1–S17 | Protocol definitions, sensitivity analyses, audits and the S17 point estimates |
| Supplementary Data SD01–SD22 | Machine-readable evidence, commands, tests and file bindings |

The reader distribution includes the study outputs and the authoritative figure
and table mappings. Open `Main_Table_Source_Index_Final_20261109.csv`,
`Supplementary_Source_Index_Final_20261109.csv` and
`Figure_Source_Index_Final_20261108.csv` at its top level. Source paths in
these indices are relative to its `submission_evidence_bundle_20261108/`
directory. Main Table 3 gives the four use cases. Supplementary Table S17 gives
the twelve point estimates for the PRISM→HMP2 example.

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


