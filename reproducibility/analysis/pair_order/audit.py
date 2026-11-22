"""Account for the change in Nielsen participant AUROC by case/control sampling history."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from stratum_diagnostics import run as check_strata


STRATA = (
    ("single_single", False, False, 3584),
    ("repeat_single", True, False, 14784),
    ("single_repeat", False, True, 192),
    ("repeat_repeat", True, True, 792),
)


def half_credit_units(prob, cases, controls):
    difference = prob[cases, None] - prob[None, controls]
    return 1 + np.sign(difference).astype(int)


def run(bundle, out):
    files = {
        "registry": "evidence/full_reproducibility_inputs/prediction_registry.csv",
        "draws": "evidence/full_reproducibility_inputs/shared_bootstrap_weights.csv",
        "order": "analysis/unit/input/retained_sampling_order.csv",
        "scorer": "validation/replay_registered.py",
        "figure": "publication/Figure5c_Nielsen.csv",
    }
    load = lambda key: pd.read_csv(bundle / files[key], float_precision="round_trip")
    predictions = load("registry")
    draws = load("draws")
    order = load("order").query("target == 'NielsenHB_2014'").set_index("sample_id").retained_order
    displayed = load("figure").set_index(["direction", "model"])
    spec = importlib.util.spec_from_file_location("scorer", bundle / files["scorer"])
    score = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(score)
    summary, details = [], []
    max_error = 0.0
    for (direction, model), fig in displayed.iterrows():
        selected = predictions.loc[(predictions.direction == direction) & (predictions.model == model)]
        multiplicities = draws[draws.direction == direction]
        samples, people, groups, _, pw, replicates = score.direction_inputs(selected, multiplicities)
        ix = samples.sample_id
        probabilities = selected.pivot(index="split", columns="sample_id", values="probability").loc[:, ix]
        ranks = order.loc[ix].to_numpy()
        sizes = np.array([len(group) for group in groups])
        if not all(sorted(ranks[g]) == list(range(1, len(g)+1)) for g in groups):
            raise ValueError("Missing or ambiguous within-person order")
        first = np.array([g[ranks[g].argmin()] for g in groups])
        label = people.outcome_ibd.to_numpy()
        cases, controls = np.flatnonzero(label == 1), np.flatnonzero(label == 0)
        if probabilities.shape != (10, 396) or probabilities.isna().any().any() or len(replicates) != 500:
            raise ValueError("Incomplete registered panel")
        delta = np.zeros((len(cases), len(controls)), dtype=int)
        before, after = [], []
        for row in probabilities.to_numpy():
            earliest = row[first]
            all_specimens = np.array([score.mean_probability(row[g]) for g in groups])
            delta += half_credit_units(all_specimens, cases, controls) - half_credit_units(earliest, cases, controls)
            before.append(score.score(label, earliest, pw))
            after.append(score.score(label, all_specimens, pw))
        # Independent whole-panel comparison is calculated by sorting scores into tied groups.
        registered = score.average_scores(after, label, pw)[0] - score.average_scores(before, label, pw)[0]
        denominator = pw[:, cases].sum(axis=1) * pw[:, controls].sum(axis=1)
        terms = {}
        for name, rep_case, rep_ctrl, expected_n in STRATA:
            mask_case = (sizes[cases] > 1) == rep_case
            mask_ctrl = (sizes[controls] > 1) == rep_ctrl
            block = delta[np.ix_(mask_case, mask_ctrl)]
            if block.size != expected_n or (name == "single_single" and np.any(block)):
                raise ValueError("Unexpected pair opportunities")
            wc = pw[:, cases[mask_case]]
            wn = pw[:, controls[mask_ctrl]]
            numerator = ((wc @ block) * wn).sum(axis=1)
            terms[name] = numerator / (2 * len(probabilities) * denominator)
        total = sum(terms.values())
        residual = float(np.abs(registered - total).max())
        max_error = max(residual, max_error)
        np.testing.assert_allclose(total, registered, atol=1e-12, rtol=0)
        np.testing.assert_allclose(
            [total[0], *np.quantile(total[1:], [.025, .975])],
            [fig.all_minus_earliest, fig.paired_ci_low, fig.paired_ci_high], atol=1e-12, rtol=0,
        )
        cols = dict(direction=direction, model=model)
        summary.append(dict(**cols, **{k: v[0] for k, v in terms.items()}, total=total[0],
                            max_abs_residual=residual, accepted_draws=len(replicates)))
        details.append(pd.DataFrame(dict(**cols, replicate=replicates,
                          **{k: v[1:] for k, v in terms.items()}, total=total[1:], sorted_auc_difference=registered[1:])))
    assert len(summary) == 12
    out.mkdir(parents=True, exist_ok=False)
    pd.DataFrame(summary).to_csv(out / "contributions.csv", index=False)
    pd.concat(details).to_csv(out / "draws.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    audit = dict(status="PASS", opportunity_counts={k: n for k, _, _, n in STRATA}, total_pairs=19352,
                 cells=12, accepted_draws_per_cell=500, full_target_and_draws_checked=6012,
                 max_abs_residual=max_error, sha256={f: hashlib.sha256((bundle / f).read_bytes()).hexdigest() for f in files.values()})
    (out / "checks.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf8")
    check_strata(bundle, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.bundle.resolve(), args.output.resolve())
