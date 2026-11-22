# Public reproducibility core

This directory contains the scoring code, input schemas, executable fixtures,
result-level tables, figure source CSVs and validation records. Participant-level
prediction registries, participant multiplicities, specimen order and cohort
matrices are not redistributed here. Obtain permitted cohort inputs from their
original providers and use the schemas and source-reconstruction instructions.

## Quick checks

Run from this directory with Python 3.11 and requirements-reader.txt:

    python -m pip install -r requirements-reader.txt
    python examples/custom_panel/verify_example.py --output ../reader_custom
    python validation/contract_fixture/validate_fixture.py
    python validation/contract_fixture/validate_adversarial_fixture.py

The protocol-control replay additionally accepts a permitted prediction registry
through `--registry`; its default path is intentionally outside this public core.

The publication directory contains source CSVs for the figures and manuscript
tables. The pair_order regional_results directory contains the 24-row
full-target country split, the 48-row Spain-only strata result and the
36,000 draw-level regional records. The orthogonal_audit directory contains
cell summaries, paired contrasts, rankings and validation records.

Source reconstruction uses the root repository scripts and permitted HMP2/IBDMDB,
PRISM and NielsenHB_2014 inputs from PRJNA398089, PRJNA400072 and PRJEB1220.
The manuscript source implementation is pinned to commit
5885481c66f29356e4ca65ba801b64b5d3e2c30a. Record runtime, configuration, input
and output hashes in the supplied manifest template.

The public fixtures test scoring, interface validation and fail-closed behavior
without participant-level cohort files. The public core does not redistribute
raw participant data, access credentials, fitted MLP weights or manuscript
documents.
