# S17 complete-target example

This case evaluates the ten PRISM→HMP2 shallow-MLP fits on 385 samples from 105 participants. The target participant multiplicities in `input/shared_participant_multiplicities.csv` are all one. Nonuniform weights, multiple replicates and a direction mismatch are rejected at the entry point. To score registered bootstrap weights use the weighted replay in [the reader guide](../READER_GUIDE.md).

From the extracted evidence bundle directory run

```shell
python -B protocol_control_20261109/check_example.py --output-dir ../reader_results/s17
```

The point scorer uses scikit-learn to evaluate AUROC, step-integral AP and Brier score. The independently coded reference uses the Python standard library and handles probability ties explicitly. The checker verifies exact membership against the registered prediction table, compares all 12 values within absolute tolerance `1e-12`, and checks the recorded diagnostic for each input violation. See `output/` for the archived manuscript reference outputs and `verification/` for the reader replay and validator results. The command writes the reader's own results under `../reader_results/s17/`.

For individual steps, run these commands from this protocol-control directory:

```shell
python -B run_protocol_control.py --predictions input/predictions_long.csv --weights input/shared_participant_multiplicities.csv --spec input/scoring_specification.json --output-dir ../../reader_results/s17
python -B reference_impl.py --predictions input/predictions_long.csv --output ../../reader_results/s17/reference_outputs.csv
python -B validate_fail_closed.py --predictions input/predictions_long.csv --weights input/shared_participant_multiplicities.csv --spec input/scoring_specification.json --cases invalid_cases/cases.json --report ../../reader_results/s17/fail_closed_report.csv
python -B validate_fail_closed.py --predictions input/predictions_long.csv --weights input/shared_participant_multiplicities.csv --spec input/scoring_specification.json --cases invalid_cases/interface_cases.json --report ../../reader_results/s17/interface_checks.csv
```

The full check also compares the fixture to `../evidence/full_reproducibility_inputs/prediction_registry.csv` and reads the archived expected scores. A row is identified by run and sample inside this single direction and model; sample-to-participant and endpoint mappings must remain consistent across runs.
