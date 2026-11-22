# Source reconstruction

Use the root repository scripts with permitted HMP2/IBDMDB, PRISM and
NielsenHB_2014 inputs obtained from PRJNA398089, PRJNA400072 and PRJEB1220.
Keep participant-level matrices, identifiers, specimen order and labels outside
the repository. Prepare files according to `data/input_schema.md`, record their
SHA-256 values and generate the 126-coordinate feature contract before target
projection.

The source-side implementation is pinned to commit
5885481c66f29356e4ca65ba801b64b5d3e2c30a. Source reconstruction and prediction
registry checks are separate steps. The public core supplies the scoring schema,
synthetic custom-panel example, fail-closed fixtures, result-level summaries and
regional aggregate tables needed to check the post-fit interface without private
cohort files.
