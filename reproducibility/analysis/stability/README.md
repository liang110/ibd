# Participant bootstrap and configuration sensitivities

Run these commands from the extracted evidence directory containing `validation/replay_registered.py` and `evidence/full_reproducibility_inputs/`. The example assumes this folder is installed as `analysis/stability/`. If installed elsewhere, replace that path in the commands. Use the Python environment listed in the evidence bundle.

```shell
python -B analysis/stability/run_stability.py --bundle . --out ../reader_stability
python -B analysis/stability/run_sensitivities.py --bundle . --inputs analysis/stability/inputs --out ../reader_stability
python -B analysis/stability/make_display_sources.py --bundle . --results ../reader_stability
```

The primary analysis uses the original 500 shared participant multiplicity vectors per direction. The sensitivity analysis copies them and appends 1,500 draws per direction. The draws use independent uniform sampling of target participants with replacement, accepting draws with positive weight in both endpoint classes. Directions in lexicographic order receive separate PCG64 streams with NumPy SeedSequence entropy [20261121, direction index]. The same vector serves every reporting object, source run and model in its direction. The first 252,000 metric records match the primary release exactly. In the compressed draw table, replicate −1 denotes the complete-target point, 0–499 the primary sample and 0–1999 the expanded sample. Confidence limits are pointwise 2.5th and 97.5th percentiles with linear interpolation. A difference uses matched participant draws.

| Output in `results/` | Interpretation |
|---|---|
| `cell_summaries.csv`, `paired_contrasts.csv` | All 504 reporting-object cells and 630 within-model contrasts, each summarized at both draw counts. The point is evaluated on the complete eligible target. |
| `fixed_pairs.csv`, `fixed_brier_pairs.csv` | Fixed-model AUROC and Brier contrasts. The Brier table uses the primary 500 draws. |
| `full_target_summary_ranks.csv`, `full_target_leaders.csv`, `primary_auroc_leader_frequencies.csv` | All complete-target rankings, leaders within each trainable-model set, and the primary resampling frequency with which each six-family complete-target leader was first. |
| `leader_draw_changes.csv` | Within-draw change in the AUROC leader when object A is replaced by B, C or D. Exact ties follow the registered evaluator’s fixed model order. The two candidate sets are the six trainable families and the five families remaining after removal of the shallow MLP. |
| `bootstrap_interval_status_changes.csv`, `bootstrap_stability_overview.json` | Five paired intervals switched from excluding zero with 500 draws to containing zero with 2,000. The five illustrative fixed-family AUROC comparisons retained their zero-inclusion status. |
| `feature_S7.csv`, `feature_cell_comparisons.csv` | Five source-side feature specifications scored against their common sensitivity-runner reference (cap 126 and prevalence 0.10). A separate field gives the combined configuration difference relative to the manuscript predictions. Each average gives equal weight to the ten source runs and then the 42 direction–model cells. |
| `retained_mass_S8.csv`, `retained_mass_cells.csv`, `retained_mass_denominators.csv` | Complete-target cell A scores after restricting specimens. Each full-target participant vector is projected to eligible specimens and retained for that threshold if both classes remain; the point is evaluated before resampling. The S8 mean weights each of the 42 direction–model cells equally. |
| `mlp_cells.csv`, `mlp_contrasts.csv`, `mlp_run_sd.csv` | Complete-target MLP scores and matched 500-draw intervals for the manuscript and disabled-internal-stopping configurations, their paired differences and the sample standard deviation (ddof=1) of ten run-specific AUROCs. |

The `inputs/` directory contains 2,100 full-target feature-branch run summaries, 20,020 paired MLP probabilities for the disabled-stopping configuration, 1,001 specimen retained fractions and fitted MLP stopping flags. Its source audit records file digests and a direct recalculation of all 2,100 × 3 feature-branch metrics from source predictions. `prepare_sensitivity_inputs.py` regenerates those compact files when the original project archive is locally available. The three reader commands above need only the distributed data.
