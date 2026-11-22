# IBD microbiome transferability benchmark

This code package supports the cross-cohort IBD post-fit evaluation protocol.
The source-side commands rebuild model predictions from permitted cohort inputs.
The public reproducibility core contains executable fixtures, result-level source tables,
validation scripts and the source index. Permitted participant-level inputs remain outside
the repository and are obtained from the original providers.

## Scope

The benchmark uses 126 ordered canonical taxonomic coordinates from the
HMP2–Franzosa feature intersection, six directed transfers and ten source-side
perturbations. Four principal reporting objects cross the scoring unit (sample or participant)
with the prediction rule (metrics averaged across fitted predictors or metric
computed once on averaged probabilities). A fifth rule (E) keeps each specimen
prediction while giving each participant equal total weight on the complete target.
The 500 shared target-participant
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
retained 385 collections. The manuscript-specific wrapper below prepares and
checks this panel. The general contract builder and optional metabolome audit
remain available for other input populations.

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
3. Put `manual_abundance.csv`, `metabolites.csv`, `clean_metadata.csv`,
`prism_abundance_clean.csv` and `prism_metadata_clean.csv` in one private flat
directory, e.g. `/path/to/cohort_inputs/`. Provide the original IBDMDB metadata
CSV by explicit path, and keep the Nielsen projection and manifest separate.
4. Run the manuscript-specific contract command below in a **new private output
directory** before using its contract for either benchmark command. The wrapper
intersects microbial and metabolomic specimen identifiers with clean metadata
and verifies that each matched HMP2 record maps to a distinct provider collection
(`site_sub_coll`). It projects the matched records onto the shared vocabulary with
the training loader, then checks specimen IDs, participant IDs and labels against
both HMP2 target directions in the permitted registered inputs. The fixed commit's generic builder
on unfiltered inputs includes 1,638 HMP2 profiles. Its optional metabolome audit
writes a cohort manifest per modality to the same filename. Use the dedicated
matched-panel command for the manuscript counts.
5. To run the public synthetic case study without refitting, use the commands under
   [Registered-result replay](#registered-result-replay). Full registered-result
   replay uses permitted participant-level inputs described below.

## Registered-result replay

The public reproducibility core is under reproducibility/. It contains the
environment file, executable protocol-control and custom-panel examples,
fail-closed fixtures, result-level tables, figure source CSVs, regional audit
outputs and public_core_manifest.csv. Run these checks from that directory:

    python -m pip install -r requirements-reader.txt
    python examples/custom_panel/verify_example.py --output ../reader_custom
    python validation/contract_fixture/validate_fixture.py
    python validation/contract_fixture/validate_adversarial_fixture.py

The declaration, SIAMCAT adapter and regional audit scripts are included for
use with permitted cohort inputs. The protocol-control implementation is the
source-side reference for the registered control and requires its permitted
fixture and participant-level registry. The full registered-result replay
requires the participant-level prediction and multiplicity inputs described by
the input schemas and the original provider access terms.

## Source reconstruction commands

```bash
python scripts/prepare_manuscript_contract.py \
  --legacy-root /path/to/cohort_inputs \
  --provider-metadata /path/to/cohort_inputs/hmp2_metadata_2018-08-20.csv \
  --prediction-registry /path/to/reader/evidence/full_reproducibility_inputs/prediction_registry.csv \
  --output-dir /path/to/private/contract

python scripts/run_standard_model_benchmark.py \
  --contract-dir /path/to/private/contract \
  --legacy-root /path/to/cohort_inputs \
  --nielsen-projection /path/to/nielsen_projection.csv \
  --nielsen-manifest /path/to/nielsen_manifest.csv \
  --allow-ungrouped-early-stopping \
  --output-dir results/standard_benchmark

python scripts/run_reciprocal_nielsen_standard_benchmark.py \
  --contract-dir /path/to/private/contract \
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

The source-side wrapper and its input schema are available in the repository root for use with the installed source package at the pinned commit.
`cohort_manifest_hmp2.csv` contains the 388 matched records supplied to the
training loader. `hmp2_retained_manifest.csv` contains 385 positive-mass
records, and `hmp2_exclusions.csv` lists the two CD and one control record
excluded at projection. `hmp2_profile_eligibility.csv` records the status of
each inventory profile. Keep these files private. The aggregate
`hmp2_selection_check.json` records the 126-coordinate digest, original input
digests and successful comparison to the target registry.

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

The reproducibility core under reproducibility/ contains the public source CSVs
for figures, tables and Supplementary Data, together with result-level checks.
public_core_manifest.csv records each included file, its schema-level role and
SHA-256 digest. Cohort matrices, participant identifiers, specimen order and
participant multiplicities are obtained through the original providers and are
not redistributed in this repository.

| Manuscript item | Where to start |
|---|---|
| Main Table 2 and Figure 1 | Choose the four principal scoring rules or the participant-weighted specimen rule E. |
| Main Table 3 and Table S17 | Inspect the worked-example result tables and the executable fixtures. |
| Figures 2 and 3 | Check sampling denominators and the two changes along A to E to C. |
| Figure 4 | Compare conditional changes in the leading family. |
| Figure 5, Tables S20 and S23, SD24 | Inspect specimen eligibility and the Nielsen source-by-family contrasts. |
| Tables S16, S19 and S21 | Locate custom-input commands, E metric/weight checks and the native SIAMCAT example. |
| Table S24 and SD50 | Use reproducibility/analysis/pair_order/regional_results for the full-target country split and Spain-only strata. |

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
ranges. The manuscript source-reconstruction commands are `prepare_manuscript_contract.py`,
`run_standard_model_benchmark.py` for four directions,
`run_reciprocal_nielsen_standard_benchmark.py` for two directions, and
`bootstrap_grouped_predictions.py --aggregate-columns direction`. For scoring from existing predictions and weights, see
[Registered-result replay](#registered-result-replay).
The public repository contains the analysis code, fixtures, result-level tables
and source index. Manuscript documents and administrative metadata are maintained
with the submission materials.

The default license is intentionally not assumed. Replace
`LICENSE_PLACEHOLDER.txt` with the authors' approved license before public
release.


## Explicit label-contract rules

HMP2 and PRISM use explicit CD/UC versus Control/nonIBD allow-lists; Nielsen uses IBD versus healthy allow-lists. Unknown, missing, conflicting, and cohort-incompatible labels fail closed. The shallow MLP is the manuscript comparator for post-fit analysis. Its reported run uses sample-level early stopping; the public default remains disabled because that internal validation split is not participant-grouped.

Manifest construction and post-fit participant aggregation both enforce the
participant-label invariant: one participant cannot map to more than one binary
outcome. The regression suite includes a conflicting-label fixture and fails
before any ``first()`` aggregation can occur.


