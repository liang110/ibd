# Minimal reproducibility interface

The public interface contains one executable prediction-to-object example, fail-closed validation fixtures and the result-level tables used by the manuscript.

Run from the repository root:

    python -m pip install -r requirements.txt
    python reproducibility/examples/custom_panel/verify_example.py --output reader_custom
    python reproducibility/validation/contract_fixture/validate_fixture.py
    python reproducibility/validation/contract_fixture/validate_adversarial_fixture.py

The custom example uses four participants, seven ordered specimens and two prediction runs. It returns the A–E objects, paired contrasts and an independent reference comparison. The validation fixtures test incomplete mappings and conflicting participant labels.

The source workflow accepts permitted HMP2/IBDMDB, Franzosa and NielsenHB_2014 inputs. Raw cohort matrices, participant-level identifiers, prediction registries and bootstrap multiplicities remain subject to the original providers’ access terms.

Result-level outputs for the manuscript are in `reproducibility/publication/` and `reproducibility/supplementary_source_data/`.
