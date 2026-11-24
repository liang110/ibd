# IBD microbiome post-fit evaluation contract

This repository contains the source implementation and a minimal reader for the post-fit evaluation contract. The reader scores supplied prediction probabilities by specimen or participant and checks specimen order, participant labels and required fields.

## Minimal reader example

Run from the repository root:

    python -m pip install -r requirements.txt
    python reproducibility/examples/custom_panel/verify_example.py --output reader_custom
    python reproducibility/validation/contract_fixture/validate_fixture.py
    python reproducibility/validation/contract_fixture/validate_adversarial_fixture.py

The interface returns four primary reporting objects, A–D. E is an equal-participant-weighted specimen object used to decompose the A→E→C contrast. It is emitted for this audit and is not a fifth primary claim. The example writes scores and paired contrasts, while the validation fixtures reject malformed mappings, conflicting labels and other incompatible requests.

The public frozen-prediction replay uses the same prediction and participant-multiplicity fields:

    python scripts/replay_registered.py --predictions reproducibility/examples/custom_panel/predictions.csv --weights reproducibility/examples/custom_panel/weights.csv --output-dir reader_replay

The replay writes cell summaries, cell-level draws, paired contrasts and rankings. This fixture requires no restricted cohort files.

## Inputs and outputs

The reader expects direction, model, run, specimen, participant, binary outcome and positive-class probability fields. Replay additionally accepts replicate-level participant multiplicities. Result-level tables used by the manuscript are in `reproducibility/publication/` and `reproducibility/supplementary_source_data/`. Cohort matrices, participant identifiers, specimen order, prediction registries and participant multiplicities are obtained or reconstructed under the original providers’ access terms and are not redistributed here.

## Source reconstruction and prediction replay

For source reconstruction, obtain the HMP2/IBDMDB, Franzosa/PRISM and NielsenHB_2014 inputs under the original provider terms and run the permitted local workflow:

    python scripts/prepare_manuscript_contract.py --source-root <permitted-source-root> --provider-metadata <provider-metadata.csv> --prediction-registry <prediction-registry.csv> --output-dir <private-contract-dir>
    python scripts/run_standard_model_benchmark.py --contract-dir <private-contract-dir> --source-root <permitted-source-root> --nielsen-projection <nielsen-projection.csv> --nielsen-manifest <nielsen-manifest.csv> --output-dir <private-results>

This path refits the declared source models and requires permitted local cohort inputs. A frozen-prediction replay starts from an existing permitted prediction registry, retained-order table and participant multiplicity table, then applies the same declaration and validation interface without refitting a model. Both paths use the same A–D interface and emit the diagnostic E object for the A→E→C decomposition.

The published transfer analyses use a fixed 126-coordinate unlabeled transductive feature contract. This vocabulary is fixed without disease labels but is not a source-only feature discovery step. A strict source-only analysis must define and report a separate vocabulary.

`requirements.txt` pins the public reader environment. The full source benchmark uses the package dependencies in `pyproject.toml` and its optional benchmark extras.

## License

MIT.
