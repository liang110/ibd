# SIAMCAT holdout probabilities under three specimen-eligibility rules

Use the Python environment specified at the bundle root. From the extracted bundle, run

```text
python -B analysis/native_windows/score_native_windows.py --bundle . --output ../reader_native_windows
```

The program reads the ten native holdout probability files in `analysis/import_case/native`. It joins each specimen to its target label and participant in the Franzosa→HMP2 registry, matches the HMP2 collection order, and reuses the 500 saved participant multiplicities for that direction. Each output run is the mean of five source-fold predictions within one source split. IBD is the positive class. The data cover the same 79 IBD and 26 control participants under each specimen rule.

The three declarations select the earliest retained specimen, the first two where available, or all retained specimens. For each rule the program invokes the declaration entry point with explicit prediction, weight, order and output paths. A retains equal specimen weights, E gives each participant equal total specimen weight, and C averages the eligible probabilities within each participant. Each computes the metric within a run and then averages runs. The archived numeric outputs use full-target estimates and linear percentile intervals across participants.

| File in `results/` | Contents |
|---|---|
| `mapped_predictions.csv` | Join of 3,850 native predictions to the registered direction, outcome and participant. |
| `earliest_retained.json`, `first_two_retained.json`, `all_retained.json` | Three valid declarations. |
| `scores.csv` and `score_draws.csv` | 27 score summaries and their 13,500 bootstrap values. |
| `within_eligibility_contrasts.csv` | 27 paired contrast summaries (C−A, E−A, C−E). |
| `eligibility_contrast.csv` and `eligibility_contrast_draws.csv` | Three C(all)−C(earliest) summaries and the 1,500 within-draw differences. |
| `independent_checks.csv` | Maximum errors against an independently constructed target and reference calculation by eligibility, object and metric. Includes the complete target and all 500 shared draws. |
| `declaration_checks.csv` and `incompatible_declaration.json` | Actual CLI exit codes for three valid and two invalid requests. Both rejected requests use the same mapped predictions. |
| `validation.json` | Denominators, numerical residuals and input SHA-256 values. |

For C(all)−C(earliest), each replicate reuses one participant multiplicity vector across both rules. The reference calculation joins the raw exports separately, assigns the eligible specimen count per participant, and integrates the ROC and step precision–recall curves at the distinct probability thresholds. For diagnostics it compares each of the 27 score configurations on the full target and all 500 draws. The validation test for missing order records passes a header-only order file. The returned nonzero exit code and absence of score and contrast files are both checked in fresh directories.
