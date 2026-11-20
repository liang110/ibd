# Reproducibility protocol

## Replay the registered results

The four-object evaluator reads two files from the accompanying reader bundle:
`evidence/full_reproducibility_inputs/prediction_registry.csv` and
`evidence/full_reproducibility_inputs/shared_bootstrap_weights.csv`. From the Git
checkout, run `scripts/replay_registered.py` with these files as `--predictions`
and `--weights`, and select a directory with `--output-dir`. It reports AUROC,
average precision and Brier score for objects A–D using each direction's saved
participant multiplicities. The [README replay guide](../README.md#registered-result-replay)
provides the commands for this analysis, the independent S17 example and the
comparison against the four archived result tables.

## Reconstruct predictions from permitted cohort inputs

1. Obtain the public cohort resources under their original access and license
   terms.
2. Record the exact file paths, release identifiers and SHA-256 values.
3. Build the HMP2–PRISM feature contract before reading Nielsen labels or target
   abundance statistics into any fitting step.
4. Prepare the Nielsen projection with `prepare_nielsenhb_2014_target.py` and
   retain its audit JSON.
5. Run `run_standard_model_benchmark.py` for the four HMP2/PRISM and Nielsen-target
   edges, then run `run_reciprocal_nielsen_standard_benchmark.py` for the two
   Nielsen-source edges and merge the prediction tables using their declared
   schemas. The source 80:20 participant split is a perturbation for model
   selection; the grouped inner CV is the source-side selection layer; the final
   refit uses all eligible source observations.
6. To apply the four-object protocol to these predictions, prepare aligned
   participant multiplicities and use the replay command above. For a separate
   grouped-bootstrap analysis, `bootstrap_grouped_predictions.py` generates its
   own draws and writes `shared_target_bootstrap_weights.csv`.
   `paired_model_bootstrap.py` is a separate marginal model-comparison diagnostic.
7. Run the leakage, participant-unit and calibration audits on the frozen
   predictions. These audits do not refit models unless the command explicitly
   states that it is a source-only diagnostic.
8. Compare row counts, column names, hashes and estimand labels against the
   result registry before producing tables or figures.
