# Minimal reproducibility interface

The public interface contains one executable prediction-to-object example, fail-closed validation fixtures and the result-level tables used by the manuscript.

Run from the repository root:

    python -m pip install -r requirements.txt
    python reproducibility/examples/custom_panel/verify_example.py --output reader_custom
    python reproducibility/validation/contract_fixture/validate_fixture.py
    python reproducibility/validation/contract_fixture/validate_adversarial_fixture.py

The custom example uses four participants, seven ordered specimens and two prediction runs. It writes the primary A–D scores, the diagnostic E score and paired contrasts, compares them with fixed expected outputs and checks them with an independent grouped calculation. The validation fixtures test incomplete mappings, conflicting participant labels and other fail-closed conditions.

The same fixture can be replayed from its prediction and participant-multiplicity files:

    python scripts/replay_registered.py --predictions reproducibility/examples/custom_panel/predictions.csv --weights reproducibility/examples/custom_panel/weights.csv --output-dir reader_replay

The replay writes cell summaries, cell-level draws, paired contrasts and rankings.

## Reproduction workflow

For source reconstruction, obtain the HMP2/IBDMDB, Franzosa/PRISM and NielsenHB_2014 inputs under the original provider terms and run:

    python scripts/prepare_manuscript_contract.py --source-root <permitted-source-root> --provider-metadata <provider-metadata.csv> --prediction-registry <prediction-registry.csv> --output-dir <private-contract-dir>
    python scripts/run_standard_model_benchmark.py --contract-dir <private-contract-dir> --source-root <permitted-source-root> --nielsen-projection <nielsen-projection.csv> --nielsen-manifest <nielsen-manifest.csv> --output-dir <private-results>

This path refits the declared source models and requires permitted local cohort inputs. For frozen-prediction replay, start from an existing permitted prediction registry, retained-order table and participant multiplicity table, validate their keys, labels and probabilities and run the replay command above. This path performs no model fitting. Both paths use the same A–D interface and emit the diagnostic E object for the A→E→C decomposition.

The manuscript transfer analyses use a fixed 126-coordinate unlabeled transductive feature contract. It is not a source-only vocabulary. A strict source-only feature analysis requires a separately declared vocabulary.

## Input names

The public reader uses `direction`, `model`, `run_id`, `sample_id`, `participant_id`, `y_true` and `y_prob`. The source-side adapter preserves direction, model, specimen and participant keys, row counts, probabilities and non-negative multiplicities. Replay inputs additionally require specimen order and participant multiplicities.

## Data and result locations

The source workflow accepts permitted HMP2/IBDMDB, Franzosa/PRISM and NielsenHB_2014 inputs. Raw cohort matrices, participant identifiers, prediction registries, specimen-order mappings and bootstrap multiplicities remain subject to the original providers’ access terms. `requirements.txt` pins the public reader environment, while the complete source benchmark uses the package dependencies in `pyproject.toml` and optional benchmark extras. Result-level outputs for the manuscript are in `reproducibility/publication/` and `reproducibility/supplementary_source_data/`.
