"""Split the Nielsen C AUROC change by control origin and by sampling history in Spain."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "validation"))
from replay_registered import average_scores, direction_inputs, mean_probability, score


REGISTRY = "evidence/full_reproducibility_inputs/prediction_registry.csv"
WEIGHTS = {"full": "evidence/full_reproducibility_inputs/shared_bootstrap_weights.csv",
           "ESP": "evidence/nielsen_esp_audit/esp_only_shared_weights.csv"}
ORDER = "analysis/unit/input/retained_sampling_order.csv"
S24 = "analysis/pair_order/results/contributions.csv"
S24_DRAWS = "analysis/pair_order/results/draws.csv.gz"
SD24 = "analysis/unit/results/window_paired.csv"


def read(path):
    return pd.read_csv(path, float_precision="round_trip")


def interval(draws):
    finite = draws[np.isfinite(draws)]
    return np.quantile(finite, [0.025, 0.975]) if len(finite) else (np.nan, np.nan)


def analyze_panel(panel, weights, order, scope, direction, model):
    samples, people, groups, _, w, rep = direction_inputs(panel, weights)
    person_origin = panel[["bootstrap_cluster", "origin"]].drop_duplicates().set_index("bootstrap_cluster")
    if not person_origin.index.is_unique:
        raise ValueError("A participant has conflicting origins")
    origins = person_origin.loc[people.bootstrap_cluster, "origin"].to_numpy()
    repeated = np.array([len(g) > 1 for g in groups])
    ranks = order.loc[samples.sample_id].to_numpy()
    if any(sorted(ranks[g]) != list(range(1, len(g) + 1)) for g in groups):
        raise ValueError("Missing or ambiguous specimen order")
    labels = people.outcome_ibd.to_numpy(int)
    if set(labels) != {0, 1} or set(origins[labels == 1]) != {"ESP"}:
        raise ValueError("Inconsistent case labels or origins")
    matrix = panel.pivot(index="split", columns="sample_id", values="probability").loc[:, samples.sample_id]
    expected_n = (318, 396) if scope == "full" else (141, 219)
    if (len(people), len(samples)) != expected_n or matrix.shape != (10, len(samples)) or len(rep) != 500:
        raise ValueError("Unexpected size of the registered Nielsen target")
    probability = matrix.to_numpy()
    if not np.isfinite(probability).all() or np.any((probability < 0) | (probability > 1)):
        raise ValueError("A registered probability is missing or outside [0,1]")
    first_idx = [g[np.argmin(ranks[g])] for g in groups]
    case, ctrl = np.flatnonzero(labels == 1), np.flatnonzero(labels == 0)
    changes = np.zeros((len(case), len(ctrl)), dtype=np.int64)
    after, before = [], []
    for row in probability:
        early = row[first_idx]
        all_retained = np.array([mean_probability(row[g]) for g in groups])
        # In half-credit units a tied comparison has zero signed order.
        changes += np.sign(all_retained[case,None] - all_retained[None,ctrl]).astype(int)
        changes -= np.sign(early[case,None] - early[None,ctrl]).astype(int)
        before.append(score(labels, early, w))
        after.append(score(labels, all_retained, w))
    whole = average_scores(after, labels, w)[0] - average_scores(before, labels, w)[0]
    denominator = w[:, case].sum(axis=1) * w[:,ctrl].sum(axis=1)
    if scope == "full":
        strata = [(country, repeated[case], ~repeated[ctrl] & (origins[ctrl] == country), n)
                  for country, n in [("DNK", 11682), ("ESP", 3102)]]
    else:
        strata = [("single_single", ~repeated[case], ~repeated[ctrl], 752),
                  ("repeat_single", repeated[case], ~repeated[ctrl], 3102),
                  ("single_repeat", ~repeated[case], repeated[ctrl], 192),
                  ("repeat_repeat", repeated[case], repeated[ctrl], 792)]
    points, records, components = [], [], []
    for name, include_case, include_ctrl, expected_count in strata:
        wc, wn = w[:, case[include_case]], w[:, ctrl[include_ctrl]]
        block = changes[np.ix_(include_case, include_ctrl)]
        mass = wc.sum(axis=1) * wn.sum(axis=1)
        if mass[0] != expected_count or (name == "single_single" and np.any(block)):
            raise ValueError("Unexpected pair count or change among singly sampled participants")
        numerator = ((wc @ block) * wn).sum(axis=1) / (2 * len(matrix))
        contribution = numerator / denominator
        mean = np.divide(numerator, mass, out=np.full_like(numerator, np.nan), where=mass > 0)
        share = mass / denominator
        np.testing.assert_allclose(np.nan_to_num(mean) * share, contribution, atol=1e-12, rtol=0)
        cl, ch = interval(contribution[1:])
        ml, mh = interval(mean[1:])
        key = dict(direction=direction, model=model, scope=scope, stratum=name)
        point = dict(**key, unique_pairs=int(mass[0]), target_pairs=int(denominator[0]),
                     opportunity_share=share[0], total_contribution=contribution[0],
                     contribution_ci_lo=cl, contribution_ci_hi=ch, net_mean=mean[0],
                     net_ci_lo=ml, net_ci_hi=mh, finite_mean_draws=int(np.isfinite(mean[1:]).sum()),
                     target_C_difference=whole[0])
        if scope == "full":
            point["share_of_repeat_case_single_control_pairs"] = mass[0] / 14784
        else:
            point["target_C_ci_lo"], point["target_C_ci_hi"] = interval(whole[1:])
        points.append(point)
        records.append(pd.DataFrame(dict(**key, replicate=rep, target_pair_mass=denominator[1:],
                                         group_pair_mass=mass[1:], opportunity_share=share[1:],
                                         total_contribution=contribution[1:], net_mean=mean[1:])))
        components.append(contribution)
    return pd.DataFrame(points), pd.concat(records, ignore_index=True), whole, rep, np.sum(components, axis=0)


def check_error(observed, expected):
    np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-12)
    return float(np.max(np.abs(np.asarray(observed) - np.asarray(expected))))


def run(bundle, output):
    paths = [REGISTRY, ORDER, *WEIGHTS.values(), S24, S24_DRAWS, SD24]
    registry = read(bundle / REGISTRY)
    order = read(bundle / ORDER).query("target == 'NielsenHB_2014'").set_index("sample_id").retained_order
    weights = {s: read(bundle / name) for s, name in WEIGHTS.items()}
    s24 = read(bundle / S24).set_index(["direction", "model"])
    parent_draws = read(bundle / S24_DRAWS).set_index(["direction","model","replicate"]).sort_index()
    sd24 = read(bundle / SD24).query(
        "scope == 'ESP' and object == 'C' and metric == 'roc_auc' and left_window == 'all' and right_window == 'first'"
    ).set_index(["direction", "model"])
    if not (len(s24) == 12 and s24.index.is_unique and sd24.index.is_unique and set(s24.index) <= set(sd24.index)):
        raise ValueError("Missing model-source reference")
    summaries = {s: [] for s in WEIGHTS}
    all_draws, errors = [], {}
    for direction, model in s24.index:
        cell = registry.loc[registry.direction.eq(direction) & registry.model.eq(model)]
        checks = {}
        for scope in WEIGHTS:
            selected = cell if scope == "full" else cell.loc[cell.origin.eq(scope)]
            w = weights[scope]
            points, draws, sorted_delta, replicates, total = analyze_panel(
                selected, w.loc[w.direction.eq(direction)], order, scope, direction, model)
            if scope == "full":
                parent = s24.loc[(direction,model)]
                parent_boot = parent_draws.loc[(direction,model)].loc[replicates]
                checks["full_RS_reconstruction"] = check_error(total, np.r_[parent.repeat_single, parent_boot.repeat_single])
                checks["full_sorted_AUROC"] = check_error(sorted_delta, np.r_[parent.total, parent_boot.sorted_auc_difference])
                inferred = parent.total - (59 / 236) * sd24.loc[(direction,model),"point"]
                checks["full_country_identity"] = check_error(points.loc[points.stratum.eq("DNK"),"total_contribution"], inferred)
            else:
                checks["Spain_four_group_sum_vs_sorted_AUROC"] = check_error(total, sorted_delta)
                checks["Spain_point_and_interval_vs_SD24"] = check_error(
                    np.r_[total[0], interval(total[1:])], sd24.loc[(direction,model),["point","ci_low","ci_high"]].to_numpy(float))
            summaries[scope].append(points)
            all_draws.append(draws)
        for key, error in checks.items():
            errors[key] = max(error, errors.get(key, 0.0))
    full, spain = (pd.concat(summaries[s], ignore_index=True) for s in ["full","ESP"])
    for field, rank in [("total_contribution","rank_abs_total"), ("net_mean","rank_abs_mean")]:
        spain[rank] = spain.groupby(["direction","model"])[field].transform(
            lambda x: x.abs().rank(ascending=False,method="min")).astype(int)
    draws = pd.concat(all_draws, ignore_index=True)
    if (len(full),len(spain),len(draws)) != (24,48,36000):
        raise ValueError("Regional output dimensions differ from the registered panel")
    output.mkdir(parents=True, exist_ok=False)
    full.to_csv(output / "rs_by_control_country.csv", index=False, float_format="%.17g")
    spain.to_csv(output / "spain_sampling_strata.csv", index=False, float_format="%.17g")
    draws.to_csv(output / "regional_draws.csv.gz", index=False, float_format="%.17g",compression={"method":"gzip","mtime":0})
    report = dict(status="PASS", absolute_tolerance=1e-12, country_rows=len(full), spain_rows=len(spain),
                  country_component_draws=12000, spain_component_draws=24000,
                  validated_draw_sums_per_scope=6000, minimum_finite_group_draws=int(draws.groupby(
                      ["scope","direction","model","stratum"]).net_mean.count().min()),
                  maximum_absolute_errors=errors, source_sha256={p:hashlib.sha256((bundle/p).read_bytes()).hexdigest() for p in paths})
    (output / "regional_checks.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.bundle.resolve(), args.output.resolve())
