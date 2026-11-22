# Score native SIAMCAT predictions

From the root of the submission evidence bundle, with the reader Python requirements installed, run

```shell
python -B analysis/import_case/run_import.py --bundle . --output reader_results/siamcat
```

The ten files under `native/` contain the original two columns exported by the R workflow for the Franzosa-to-HMP2 task. Each source seed defines a five-fold LASSO ensemble. Its output is one IBD probability per specimen after averaging the five source models. A source run in the post-fit evaluator therefore represents one five-model ensemble. The exporter creates labels with IBD as the positive class. The fitted code and session are documented in `native_export.R` and `model_session.txt`.

The importer obtains the endpoint and participant mapping from the accompanying prediction registry, joins by specimen ID and keeps the registered 500 HMP2 participant multiplicities. The expected panel contains 3,850 records, 385 specimens and 105 participants. The reference calculation reads the native CSVs again and computes weighted scores with scikit-learn. It checks all 12 complete-target metrics and the same 12 at draws 0 and 499. Positive-class orientation is tested against the archived six-decimal A metrics. Duplicate records, missing panels and inconsistent labels have separate validation checks.

The four complete-target AUROCs in A–D order are 0.736230, 0.739940, 0.732376 and 0.736611. The results directory contains the three metrics with marginal percentile intervals, the five within-model paired contrasts with their intervals, draw-level scores, interface tests and the comparison with the independent implementation. The primary worked example uses the shallow MLP. These files demonstrate transfer of the post-fit interface to a second training program within the same study.

The ten model-derived probability files contain cohort specimen keys already present in the companion registry. The worked example uses these files with the manuscript evidence and the original cohort access and reuse terms. The source abundance and metabolomics matrices, restricted provider metadata and serialized source models are outside the distributed example.

The independent cohort inventory in `results/sampling_*.csv`, `results/timing_summary.csv` and `results/hmp2_selection.json` summarizes the study populations. From the bundle root the command below reproduces the group counts, specimen frequencies and case fractions from the registered data. Rechecking the upstream HMP2 intersection and provider visit times additionally requires the permitted source matrices, provider metadata and input and exclusion manifests. The optional `audit_sampling.py --help` command lists these arguments.

```shell
python -B analysis/import_case/audit_sampling.py --bundle . --output reader_results/sampling
```
