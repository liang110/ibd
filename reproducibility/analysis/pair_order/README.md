# Case–control audit of the Nielsen participant AUROC

Run from the extracted reader archive with the Python environment listed in its root README.

```shell
python -B analysis/pair_order/audit.py --bundle . --output ../reader_pair_order
```

The command reads Nielsen predictions, participant multiplicities and retained specimen order. For each source run it compares one mean probability per participant under the all-retained and earliest-retained rules (object C). Concordance earns 1 for a higher case score, 0.5 for a tie and 0 for a lower case score. The result averages over the ten source runs.

| Output in `results/` | Contents |
|---|---|
| `contributions.csv` | Twelve complete-target rows. Each sampling-history component uses all 19,352 case–control pairs as its denominator. |
| `draws.csv.gz` | Six thousand model–source/draw records. The four components sum to the AUROC difference from the tied-score sorting implementation. |
| `stratum_diagnostics.csv` | Forty-eight model–source/stratum records. Each gives group size, share of all pairs, signed total contribution, signed within-group mean, its finite-draw count and unadjusted participant-resampling percentile interval, and ranks on both scales. The increased/unchanged/decreased counts refer to pair-by-run comparisons. |
| `stratum_draws.csv.gz` | Twenty-four thousand records with the draw-specific total and group pair mass, group share, total contribution and mean signed credit change. A case–control pair has the product of the two participant multiplicities as its weight. A group with zero mass has zero contribution and an undefined mean. |
| `checks.json`, `stratum_checks.json` | Whole-panel agreement, hash checks on the original scoring inputs, counts of pair opportunities and draw-level checks. |

The first letter of the printed stratum code describes the case, the second the control. S denotes one retained specimen and R two. The group share is its pair mass divided by total pair mass. Multiplying this share by the group's signed mean change gives its contribution to the AUROC difference. Positive and negative within-group changes can cancel. With 318 participants, the complete target contains 82 cases times 236 controls. Repeated comparisons involving the same people share observations.

The first computation combines run-level pair counts and checks the sum against the sorted AUROC calculation. The second routine in `stratum_diagnostics.py` rebuilds the two prediction matrices and participant multiplicities directly from the input rows, checks every additive component and calculates the normalized quantities. Each resampling interval uses the same 500 accepted target-participant draws for that direction. The publication index links these files to Table S24 and SD50.

### Where do the comparisons come from?

Run the regional audit from the archive root. Each command writes to a new directory.

```shell
python -B analysis/pair_order/regional.py --bundle . --output ../reader_regional
```

The two summaries in `regional_results/` use different target populations. The complete Nielsen panel contributes 19,352 case–control comparisons. In `rs_by_control_country.csv` (24 rows), the 14,784 repeat-case/single-control comparisons include 11,682 with a Danish control and 3,102 with a Spanish control. Each contribution uses **all 19,352** comparisons in the plug-in denominator. The Spain-only target in `spain_sampling_strata.csv` (48 rows) contains the same 82 cases and 59 Spanish controls. Its four contributions use **4,838** comparisons in the plug-in denominator. A group’s mean signed credit change uses the comparisons in that group.

The code repeats each calculation on the 500 registered participant draws for its own direction and target. Percentile intervals for the signed contribution and the signed mean use the corresponding draw-specific denominators. `regional_draws.csv.gz` contains 12,000 complete-target and 24,000 Spain-only component records with their weighted masses. `regional_checks.json` records checks against the original full-panel repeat-case/single-control contributions, independently sorted AUROCs, the Spain-only SD24 point estimates and interval endpoints, and input SHA-256 values.
