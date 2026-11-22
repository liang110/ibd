# Source reconstruction

Use the root repository scripts with permitted HMP2/IBDMDB, PRISM and
NielsenHB_2014 inputs obtained from PRJNA398089, PRJNA400072 and PRJEB1220.
Keep participant-level matrices, identifiers, specimen order and labels outside
the repository. Prepare the cohort files according to data/input_schema.md,
record their SHA-256 values and generate the 126-coordinate feature contract
before projection.

The manuscript source-side implementation is pinned to commit
5885481c66f29356e4ca65ba801b64b5d3e2c30a. Use
scripts/prepare_manuscript_contract.py, the standard benchmark scripts and
scripts/bootstrap_grouped_predictions.py with an explicit runtime manifest.
The manuscript comparator uses the sample-level early-stopping option
allow-ungrouped-early-stopping; the public default remains participant-grouped.
The reproducibility core supplies the scoring schemas, public fixtures and
result-level audits needed to check the post-fit interface without private
cohort files.
