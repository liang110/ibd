# IBD microbiome post-fit evaluation contract

This repository contains the source implementation and a minimal reader for the post-fit evaluation contract. The reader scores supplied prediction probabilities by specimen or participant and checks specimen order, participant labels and required fields.

## Minimal reader example

Run from the repository root:

    python -m pip install -r requirements.txt
    python reproducibility/examples/custom_panel/verify_example.py --output reader_custom
    python reproducibility/validation/contract_fixture/validate_fixture.py
    python reproducibility/validation/contract_fixture/validate_adversarial_fixture.py

The interface returns four primary reporting objects, A–D. E is an equal-participant-weighted specimen object used to decompose the A→E→C contrast. It is emitted for this audit and is not a fifth primary claim. The example writes scores and paired contrasts, while the validation fixtures reject malformed mappings, conflicting labels and other incompatible requests.

## Inputs and outputs

The reader expects direction, model, run, specimen, participant, binary outcome and positive-class probability fields. Result-level tables used by the manuscript are in `reproducibility/publication/` and `reproducibility/supplementary_source_data/`. Cohort matrices, participant identifiers, specimen order, prediction registries and participant multiplicities are obtained or reconstructed under the original providers’ access terms and are not redistributed here.

## Source reconstruction and prediction replay

The source workflow accepts permitted HMP2/IBDMDB, Franzosa/PRISM and NielsenHB_2014 inputs and can recreate the prediction records used by the evaluator. A frozen-prediction replay starts from an existing permitted prediction registry and retained-order and multiplicity information, then applies the same declaration and validation interface without refitting a model. The repository documentation describes both paths and the output locations.

## License

MIT.
