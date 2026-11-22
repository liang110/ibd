# IBD microbiome post-fit evaluation contract

This repository contains the source implementation and a small reader interface for the post-fit evaluation contract. The interface scores supplied prediction probabilities by specimen or participant and checks the declared specimen order, participant labels and required fields.

## Minimal reader example

From the repository root:

    python -m pip install -r requirements.txt
    python reproducibility/examples/custom_panel/verify_example.py --output reader_custom
    python reproducibility/validation/contract_fixture/validate_fixture.py
    python reproducibility/validation/contract_fixture/validate_adversarial_fixture.py

The custom panel accepts prediction, specimen-order and participant-weight files and writes the A–E scores, paired contrasts, fixed expected-output comparisons and an independent reference check. The validation fixtures exercise malformed mappings and conflicting labels. Input-name mapping and manuscript result-file mapping are documented in `reproducibility/README.md`, including the canonical S16, S17, S19, S21 and S24 source filenames.

## Inputs and outputs

Result-level tables used by the manuscript are in `reproducibility/publication/` and `reproducibility/supplementary_source_data/`. Cohort matrices, participant identifiers, specimen order, prediction registries and bootstrap multiplicities are obtained or reconstructed under the original providers’ access terms and are not redistributed here.

## Source reconstruction

The source scripts and input schema accept permitted HMP2/IBDMDB, Franzosa and NielsenHB_2014 inputs. The same reader interface can score the resulting probabilities without refitting a model.

## License

MIT.

