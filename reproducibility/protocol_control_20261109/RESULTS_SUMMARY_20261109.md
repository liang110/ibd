# Protocol-level reporting ambiguity control — 2026-11-09

## Decision and scope

This is a computational protocol control, not a new biological experiment and not a model-training run. It uses a fixed slice of the registered prediction registry (`PRISM_to_HMP2`, `shallow_mlp`) and changes only the post-fit reporting functional. The purpose is to test whether scoring unit and aggregation order must be declared to identify the reported estimand.

The input contains 3,850 prediction rows from 10 source runs, 385 samples and 105 participants. The fixed point-estimate control validates a shared multiplicity table but does not estimate bootstrap intervals; uncertainty remains a separate post-fit operation.

## Executed objects

- **A:** sample-level metric for each run, then mean across runs.
- **B:** mean run probability at each sample, then sample-level metric.
- **C:** participant-level metric for each run after within-participant probability averaging, then mean across runs.
- **D:** mean run probability at each sample, then within-participant averaging and participant-level metric.

## Results

| Object | AUROC | AP | Brier |
|---|---:|---:|---:|
| A | 0.635197098 | 0.814791218 | 0.233980555 |
| B | 0.728818779 | 0.858090419 | 0.194410149 |
| C | 0.653456670 | 0.844721448 | 0.222365890 |
| D | 0.749756573 | 0.882928396 | 0.189904592 |

Relative to A, B changes AUROC by +0.093621681 and Brier by −0.039570405; C changes AUROC by +0.018259572; D changes AUROC by +0.114559474 and Brier by −0.044075963. These differences occur with identical labels and prediction probabilities. They therefore demonstrate non-equivalence of the reporting functionals rather than a change in model performance.

## Independent implementation and failure checks

The canonical implementation uses scikit-learn metrics and the declared schema. `reference_impl.py` is a separate implementation that does not import the project evaluator; the updated reader reference computes AUROC by tie-aware concordance, AP by threshold-wise precision–recall integration and Brier score as mean squared probability error. The manuscript-run canonical/reference maximum absolute difference across the 12 values in output/ is 2.220e-16. The reader replay in verification/ gives an independent maximum absolute difference of 1.110e-16 and a maximum absolute difference from the archived canonical values of 2.220e-16. The acceptance tolerance is 1e-12.

The pre-existing bundle regressions also pass: `validate_fixture.py` reports `fixture PASS` and `validate_adversarial_fixture.py` reports `adversarial fixture PASS 12 values`.

Seven invalid cases were tested and all failed closed: duplicate prediction key, participant-label conflict, probability outside [0,1], incomplete run coverage, participant/multiplicity mismatch, missing scoring declaration and single-class input.

## Reproduction inputs and hashes

- `input/predictions_long.csv`: normalized long-format frozen predictions; SHA-256 is recorded in `input_manifest_sha256.csv`.
- `input/shared_participant_multiplicities.csv`: deterministic all-one multiplicity fixture used to validate the interface.
- `input/scoring_specification.json`: required columns and A/B/C/D definitions.
- `output/protocol_control_outputs.csv`: 12 object–metric values.
- `output/reference_outputs.csv`: independent 12-value output.
- `output/fail_closed_report.csv`: seven validation cases.
- `output_manifest_sha256.csv`: output file hashes.

The source registry was not refit. The control is therefore executable locally and does not require HPC resources.

## Interpretation for the manuscript

The control directly supports the protocol-level claim that a prediction table alone does not identify a unique reported score. The incremental contribution is the predeclared scoring declaration, long-format schema, shared-weight interface and validator. It does not establish that any one object is universally superior, and it does not replace SIAMCAT or other model-construction workflows.
