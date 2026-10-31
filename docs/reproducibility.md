# Reproducibility protocol

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
6. Compute common-target participant-cluster bootstrap replicates with
   `bootstrap_grouped_predictions.py`. The command writes the exact shared
   participant-weight vectors used for each aggregate cell and replicate to
   `shared_target_bootstrap_weights.csv`; use `paired_model_bootstrap.py` for
   paired model contrasts when the input table has one prediction row per model
   and target sample.
7. Run the leakage, participant-unit and calibration audits on the frozen
   predictions. These audits do not refit models unless the command explicitly
   states that it is a source-only diagnostic.
8. Compare row counts, column names, hashes and estimand labels against the
   result registry before producing tables or figures.

The archived Table 11 run is an implementation audit only. It is excluded from
quantitative robustness claims unless the original inputs can be recovered and
the main pipeline is reproduced element by element.
