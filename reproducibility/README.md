# Minimal reproducibility interface

The public interface contains one executable prediction-to-object example, fail-closed validation fixtures and the result-level tables used by the manuscript.

Run from the repository root:

    python -m pip install -r requirements.txt
    python reproducibility/examples/custom_panel/verify_example.py --output reader_custom
    python reproducibility/validation/contract_fixture/validate_fixture.py
    python reproducibility/validation/contract_fixture/validate_adversarial_fixture.py

The custom example uses four participants, seven ordered specimens and two prediction runs. It writes A–E scores and paired contrasts, compares them with fixed expected-output files and checks them with a separate grouped reference calculation. The validation fixtures test incomplete mappings, conflicting participant labels and other fail-closed conditions.

## Input names

The public reader uses `direction`, `model`, `run_id`, `sample_id`, `participant_id`, `y_true` and `y_prob`. The registered replay interface uses the following equivalent names:

| Reader | Registered replay |
|---|---|
| `run_id` | `split` |
| `y_true` | `outcome_ibd` |
| `y_prob` | `probability` |
| `participant_id` | `bootstrap_cluster` |
| `replicate_id` | `replicate` |
| `multiplicity` | `weight` |

The adapter preserves direction, model, specimen and participant keys, row counts, probabilities and non-negative multiplicities. Registered replay additionally requires the permitted prediction registry, retained-order table and participant multiplicities.

## Data and result locations

The source workflow accepts permitted HMP2/IBDMDB, Franzosa and NielsenHB_2014 inputs. Raw cohort matrices, participant-level identifiers, prediction registries, specimen-order mappings and bootstrap multiplicities remain subject to the original providers’ access terms. Result-level outputs for the manuscript are in `reproducibility/publication/` and `reproducibility/supplementary_source_data/`.

The current manuscript mapping is explicit. Main tables and figures use the files in `reproducibility/publication/`. Supplementary Table S16 uses `Supplementary_Table_S16_reader_entry_points.csv`, S17 uses `Supplementary_Table_S17.csv`, S19 uses `Supplementary_Table_S19_AEC_decomposition.csv`, S21 uses `Supplementary_Table_S21.csv` and S24 uses `Supplementary_Table_S24_Nielsen_pair_attribution.csv`. The public result directory contains the canonical table sources used by the current manuscript and the native SIAMCAT result tables.


