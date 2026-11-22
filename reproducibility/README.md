# Public reproducibility core

This directory contains the minimal public interface for the post-fit evaluation
contract. It includes scoring declarations, a synthetic prediction-to-object
example, fail-closed fixtures, result-level source tables, regional summaries,
validators and SHA-256 manifests. Real participant-level prediction registries,
specimen-order mappings, multiplicities and cohort matrices are obtained under
the original providers’ access terms and remain outside this directory.

## Quick checks

Run from this directory with Python 3.11 and `requirements-reader.txt`:

    python -m pip install -r requirements-reader.txt
    python examples/custom_panel/verify_example.py --output ../reader_custom
    python validation/contract_fixture/validate_fixture.py
    python validation/contract_fixture/validate_adversarial_fixture.py

The custom-panel command takes an explicit prediction table, specimen order and
participant weights through A–E and compares the results with an independent
calculation. The protocol-control implementation is retained as the source-side
reference for the registered control. Its fixture and participant-level registry
are supplied only under the corresponding access terms; the public runnable
example is the synthetic custom panel.

`publication/` and `supplementary_source_data/` contain result-level CSVs used by
the manuscript figures and tables. `analysis/pair_order/regional_results/` contains
the 24-row full-target country summary and the 48-row Spain-only summary.
`SOURCE_INDEX.csv` maps the public result files to their manuscript roles.
`public_core_manifest.csv` records each included path, byte count and SHA-256.

The public-core release label is `postfit-contract-core-2026-11-22`; its files are listed in `public_core_manifest.csv`.

Source reconstruction uses permitted HMP2/IBDMDB, PRISM and NielsenHB_2014
inputs from PRJNA398089, PRJNA400072 and PRJEB1220. The source-side
implementation is pinned to commit
5885481c66f29356e4ca65ba801b64b5d3e2c30a.
