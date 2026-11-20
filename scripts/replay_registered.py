"""Score registered binary predictions using saved, direction-specific participant draws."""
import argparse
import csv
from contextlib import ExitStack
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd


METRICS = ("roc_auc", "average_precision", "brier")
MODEL_ORDER = ("prevalence_baseline", "logistic_l2", "elastic_net", "linear_svm",
               "random_forest", "xgboost", "shallow_mlp")
OBJECTS = {"A": ("sample", "split_metric_mean"), "B": ("sample", "ensemble_probability"),
           "C": ("participant", "split_metric_mean"), "D": ("participant", "ensemble_probability")}
CONTRASTS = {"aggregation_order_sample": ("B", "A"), "aggregation_order_participant": ("D", "C"),
             "participant_unit_split_metric": ("C", "A"), "participant_unit_ensemble": ("D", "B"),
             "joint_D_minus_A": ("D", "A")}
KEY = ["direction", "model", "metric"]
FIELDS = {
    "orthogonal_cell_draws": KEY + ["cell", "replicate", "value", "n_target_samples", "n_target_participants"],
    "orthogonal_cell_summary": KEY + ["evaluation_unit", "source_aggregation", "cell", "plugin_point_estimate",
                                      "bootstrap_mean", "ci_low", "ci_high", "n_target_samples",
                                      "n_target_participants", "n_source_perturbations", "n_bootstrap"],
    "orthogonal_paired_deltas": KEY + ["contrast", "left_cell", "right_cell", "delta_bootstrap_mean",
                                       "ci_low", "ci_high", "ci_excludes_zero"],
    "orthogonal_rankings": ["direction", "metric", "cell", "rank", "model", "bootstrap_mean", "directionality"],
}


def mean_probability(values):
    """Correctly round an arithmetic mean of binary64 probabilities, preserving exact ties."""
    first = float(values[0])
    if all(v == first for v in values):
        return first
    return float(sum((Fraction.from_float(float(v)) for v in values), Fraction()) / len(values))


def score(y, probability, weights):
    """All rows of weights have positive mass in both classes; ties are grouped before scoring."""
    order = np.argsort(-probability, kind="stable")
    sorted_p = probability[order]
    starts = np.r_[0, np.flatnonzero(sorted_p[:-1] != sorted_p[1:]) + 1]
    ordered_weights = weights[:, order]
    positives = np.add.reduceat(ordered_weights * y[order], starts, axis=1)
    negatives = np.add.reduceat(ordered_weights * (1 - y[order]), starts, axis=1)
    tp, fp = np.cumsum(positives, axis=1), np.cumsum(negatives, axis=1)
    npos, nneg = tp[:, -1], fp[:, -1]
    if np.any(npos == 0) or np.any(nneg == 0):
        raise ValueError("A declared shared draw is missing one of the two classes")
    # A positive beats negatives with lower scores and ties half of its group's negatives.
    auc = (positives * (nneg[:, None] - fp + 0.5 * negatives)).sum(axis=1) / (npos * nneg)
    precision = np.divide(tp, tp + fp, out=np.zeros_like(tp, dtype=float), where=(tp + fp > 0))
    ap = (positives * precision).sum(axis=1) / npos
    bs = (weights * np.square(probability - y)).sum(axis=1) / weights.sum(axis=1)
    return np.stack([auc, ap, bs])


def direction_inputs(predictions, draws):
    identities = predictions[["sample_id", "bootstrap_cluster", "outcome_ibd"]].drop_duplicates()
    if identities.sample_id.duplicated().any() or identities.groupby("bootstrap_cluster").outcome_ibd.nunique().max() != 1:
        raise ValueError("A sample or participant has conflicting ownership or labels")
    samples = identities.sort_values("sample_id")
    people = samples[["bootstrap_cluster", "outcome_ibd"]].drop_duplicates().sort_values("bootstrap_cluster")
    participant_ids = pd.Index(people.bootstrap_cluster)
    slot = participant_ids.get_indexer(samples.bootstrap_cluster)
    per_person = [np.flatnonzero(slot == k) for k in range(len(people))]
    matrix = draws.pivot(index="replicate", columns="bootstrap_cluster", values="weight")
    if set(matrix.columns) != set(participant_ids) or matrix.isna().any().any():
        raise ValueError("Each draw must explicitly cover the full participant panel")
    matrix = matrix.reindex(columns=participant_ids).sort_index()
    weights = matrix.to_numpy(dtype=float)
    if (not np.isfinite(weights).all() or np.any(weights < 0) or np.any(weights != np.floor(weights))
            or np.any(weights.sum(axis=1) != len(people))):
        raise ValueError("Each registered participant draw must contain exactly U nonnegative integer draws")
    # Row zero is the full-target plug-in, not an extra bootstrap observation.
    rows = np.vstack([np.ones(len(people), dtype=int), weights.astype(np.int64)])
    return samples, people, per_person, rows[:, slot], rows, matrix.index


def four_scores(prob, y, labels, groups, stool_weights, person_weights):
    sample_results, person_results = [], []
    for run in prob:
        sample_results.append(score(y, run, stool_weights))
        patient = np.array([mean_probability(run[group]) for group in groups])
        person_results.append(score(labels, patient, person_weights))
    ensemble = np.array([mean_probability(col) for col in prob.T])
    person_ensemble = np.array([mean_probability(ensemble[group]) for group in groups])
    return {"A": average_scores(sample_results, y, stool_weights), "B": score(y, ensemble, stool_weights),
            "C": average_scores(person_results, labels, person_weights), "D": score(labels, person_ensemble, person_weights)}


def average_scores(results, labels, weights):
    """Use half-pair numerators for AUROC; preserve exact equality for identical metric inputs."""
    measured = np.stack(results)
    averaged = np.where((measured == measured[0]).all(axis=0), measured[0], measured.mean(axis=0))
    denom = 2 * (weights * labels).sum(axis=1) * (weights * (1 - labels)).sum(axis=1)
    # All runs see the same labelled sample panel and multiplicities, so the
    # concordance denominators coincide. Numerators are integers after doubling.
    numerators = np.rint(measured[:, 0] * denom)
    averaged[0] = numerators.sum(axis=0) / (len(results) * denom)
    return averaged


def describe(values):
    low, high = np.quantile(values, [0.025, 0.975], method="linear")
    return float(np.mean(values)), float(low), float(high)


def write_direction(direction, pred, draws, writers):
    samples, people, groups, sw, pw, replicates = direction_inputs(pred, draws)
    ids = sorted(pred.split.unique())
    y = samples.outcome_ibd.to_numpy(dtype=int)
    labels = people.outcome_ibd.to_numpy(dtype=int)
    ranks = {}
    for model, frame in pred.groupby("model", sort=True):
        matrix = frame.pivot(index="split", columns="sample_id", values="probability")
        matrix = matrix.reindex(index=ids, columns=samples.sample_id)
        if matrix.isna().any().any():
            raise ValueError(f"Incomplete run/sample panel for {direction}, {model}")
        result = four_scores(matrix.to_numpy(), y, labels, groups, sw, pw)
        for mi, metric in enumerate(METRICS):
            for cell, (unit, aggregation) in OBJECTS.items():
                plugin, *boot = result[cell][mi]
                mean, low, high = describe(boot)
                writers["orthogonal_cell_summary"].writerow([direction, model, metric, unit, aggregation,
                    cell, plugin, mean, low, high, len(samples), len(people), len(ids), len(replicates)])
                writers["orthogonal_cell_draws"].writerows(
                    (direction, model, metric, cell, rep, value, len(samples), len(people))
                    for rep, value in zip(replicates, boot))
                ranks.setdefault((metric, cell), []).append((mean, model))
            for name, (left, right) in CONTRASTS.items():
                differences = result[left][mi, 1:] - result[right][mi, 1:]
                mean, low, high = describe(differences)
                writers["orthogonal_paired_deltas"].writerow([direction, model, metric, name, left, right,
                    mean, low, high, low > 0 or high < 0])
    for (metric, cell), pairs in ranks.items():
        reverse = metric != "brier"
        priority = {name: i for i, name in enumerate(MODEL_ORDER)}
        for rank, (mean, model) in enumerate(sorted(pairs, key=lambda p: (
                -p[0] if reverse else p[0], priority.get(p[1], len(priority)), p[1])), 1):
            writers["orthogonal_rankings"].writerow([direction, metric, cell, rank, model, mean,
                                                     "higher" if reverse else "lower"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--weights", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    pcols = ["direction", "split", "model", "sample_id", "outcome_ibd", "probability", "bootstrap_cluster"]
    wcols = ["direction", "replicate", "bootstrap_cluster", "weight"]
    types = {"direction": str, "model": str, "sample_id": str, "bootstrap_cluster": str}
    pred = pd.read_csv(args.predictions, usecols=pcols, dtype=types, float_precision="round_trip", keep_default_na=False)
    draws = pd.read_csv(args.weights, usecols=wcols, dtype={k:v for k,v in types.items() if k in wcols}, keep_default_na=False)
    for data, keys in [(pred, ["direction", "model", "split", "sample_id"]), (draws, wcols[:-1])]:
        if data.empty or data.isna().any().any() or (data.astype(str) == "").any().any() or data.duplicated(keys).any():
            raise ValueError("Inputs have empty or duplicate keys or missing fields")
    if (not pred.outcome_ibd.isin([0, 1]).all() or not pred.probability.between(0, 1).all()
            or set(pred.direction) != set(draws.direction)):
        raise ValueError("Inputs have invalid probabilities, labels or direction sets")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.output_dir) as tmp:
        with ExitStack() as stack:
            writers = {}
            for name, fields in FIELDS.items():
                file = stack.enter_context((Path(tmp) / (name + ".csv")).open("w", encoding="utf-8", newline=""))
                writers[name] = csv.writer(file, lineterminator="\n")
                writers[name].writerow(fields)
            for direction, grp in pred.groupby("direction", sort=True):
                write_direction(direction, grp, draws[draws.direction == direction], writers)
        for name in FIELDS:
            (Path(tmp) / (name + ".csv")).replace(args.output_dir / (name + ".csv"))


if __name__ == "__main__":
    main()
