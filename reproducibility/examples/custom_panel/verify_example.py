"""Run the public A–E example against fixed outputs and an independent calculation."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

EXAMPLE = Path(__file__).resolve().parent
BUNDLE = EXAMPLE.parents[1]
RULES = {"earliest": (1, 4), "first_two": (2, 6), "all": (None, 7)}
METRICS = {"roc_auc": roc_auc_score, "average_precision": average_precision_score, "brier": brier_score_loss}
CONTRASTS = (("B", "A"), ("D", "C"), ("C", "A"), ("D", "B"), ("D", "A"), ("E", "A"), ("C", "E"))
FIELDS = ["point", "bootstrap_mean", "ci_low", "ci_high"]
ATOL = 1e-12


def metric_block(predictions, labels, weights):
    return {name: float(metric(labels, predictions, sample_weight=weights))
            for name, metric in METRICS.items()}


def formal_scores(panel, counts):
    """Formal scoring path used to write the result-level output."""
    samples = panel.drop_duplicates("sample_id").sort_values("sample_id").reset_index(drop=True)
    sample_ids = samples.sample_id.tolist()
    runs = sorted(panel.run_id.unique())
    matrix = panel.pivot(index="run_id", columns="sample_id", values="y_prob").loc[runs, sample_ids].to_numpy()
    people = list(dict.fromkeys(samples.participant_id.tolist()))
    person_index = {person: i for i, person in enumerate(people)}
    specimen_people = samples.participant_id.tolist()
    labels = samples.y_true.to_numpy(dtype=int)
    person_labels = np.asarray([samples.loc[samples.participant_id.eq(person), "y_true"].iloc[0] for person in people], dtype=int)
    person_counts = np.asarray([specimen_people.count(person) for person in people], dtype=float)
    multiplicities = counts.reindex(people).to_numpy(dtype=float)
    if np.isnan(multiplicities).any():
        raise ValueError("Missing participant multiplicity")
    person_weights = np.asarray([multiplicities[person_index[person]] for person in specimen_people])
    balanced_weights = person_weights / np.asarray([person_counts[person_index[person]] for person in specimen_people])
    person_means = np.column_stack([
        matrix[:, [i for i, person in enumerate(specimen_people) if person == target]].mean(axis=1)
        for target in people
    ])
    scores = {}
    scores["A"] = {metric: float(np.mean([metric_block(row, labels, person_weights)[metric] for row in matrix])) for metric in METRICS}
    scores["B"] = metric_block(matrix.mean(axis=0), labels, person_weights)
    scores["C"] = {metric: float(np.mean([metric_block(row, person_labels, multiplicities)[metric] for row in person_means])) for metric in METRICS}
    scores["D"] = metric_block(person_means.mean(axis=0), person_labels, multiplicities)
    scores["E"] = {metric: float(np.mean([metric_block(row, labels, balanced_weights)[metric] for row in matrix])) for metric in METRICS}
    person_variance = np.column_stack([
        matrix[:, [i for i, person in enumerate(specimen_people) if person == target]].var(axis=1)
        for target in people
    ]).mean(axis=0)
    residual = abs(scores["E"]["brier"] - scores["C"]["brier"] - float(np.average(person_variance, weights=multiplicities)))
    return scores, residual


def independent_scores(panel, counts):
    """Independent reference path using grouped pandas operations."""
    per_run = panel.pivot(index="run_id", columns="sample_id", values="y_prob")
    specimens = panel.drop_duplicates("sample_id").set_index("sample_id").loc[per_run.columns]
    person = specimens.participant_id
    person_counts = person.value_counts()
    person_labels = specimens.groupby("participant_id").y_true.first()
    counts = counts.reindex(person_labels.index)
    weights = person.map(counts)
    balanced = weights / person.map(person_counts)
    np.testing.assert_allclose(balanced.groupby(person).sum().reindex(counts.index), counts, rtol=0, atol=ATOL)
    person_means = per_run.T.groupby(person).mean().T
    person_variance = per_run.T.groupby(person).var(ddof=0).mean(axis=1)
    within_variance = np.average(person_variance, weights=counts)
    scores = {}
    for code, block, labels, weight in [
        ("A", per_run, specimens.y_true, weights),
        ("B", [per_run.mean()], specimens.y_true, weights),
        ("C", person_means, person_labels, counts),
        ("D", [person_means.mean()], person_labels, counts),
        ("E", per_run, specimens.y_true, balanced),
    ]:
        runs = block.to_numpy() if isinstance(block, pd.DataFrame) else block
        for name, metric in METRICS.items():
            scores[code, name] = float(np.mean([metric(labels, prediction, sample_weight=weight) for prediction in runs]))
    residual = abs(scores["E", "brier"] - scores["C", "brier"] - within_variance)
    return scores, residual


def summaries(values):
    draws = values[1:]
    return [values[0], np.mean(draws), *np.quantile(draws, [0.025, 0.975], method="linear")]


def validate_metadata(frame, n_specimens):
    for field, count in {"n_runs": 2, "n_specimens": n_specimens, "n_participants": 4, "accepted_draws": 3}.items():
        if not frame[field].eq(count).all():
            raise ValueError(f"Incorrect {field}")


def compare_fixed(actual_path, expected_path, keys, n_specimens):
    actual = pd.read_csv(actual_path, float_precision="round_trip")
    expected = pd.read_csv(expected_path, float_precision="round_trip")
    validate_metadata(actual, n_specimens)
    validate_metadata(expected, n_specimens)
    actual = actual.set_index(keys).sort_index()
    expected = expected.set_index(keys).sort_index()
    if not actual.index.is_unique or not expected.index.is_unique:
        raise ValueError(f"Duplicate result keys for {actual_path.name}")
    if not actual.index.equals(expected.index):
        raise ValueError(f"Expected-output keys differ for {actual_path.name}")
    errors = []
    for field in FIELDS:
        estimate = actual[field].to_numpy(dtype=float)
        target = expected[field].to_numpy(dtype=float)
        np.testing.assert_allclose(estimate, target, rtol=0, atol=ATOL)
        errors.append(float(np.max(np.abs(estimate - target))))
    return max(errors, default=0.0)


def compare_independent(actual_path, keys, values, n_specimens):
    actual = pd.read_csv(actual_path, float_precision="round_trip").set_index(keys).sort_index()
    validate_metadata(actual.reset_index(), n_specimens)
    if not actual.index.is_unique:
        raise ValueError(f"Duplicate result keys for {actual_path.name}")
    errors = []
    for key, reference in values.items():
        estimate = actual.loc[key, FIELDS].to_numpy(dtype=float)
        target = np.asarray(summaries(reference), dtype=float)
        np.testing.assert_allclose(estimate, target, rtol=0, atol=ATOL)
        errors.append(float(np.max(np.abs(estimate - target))))
    return max(errors, default=0.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New directory outside the extracted bundle")
    out = parser.parse_args().output.resolve()
    if out.is_relative_to(BUNDLE):
        raise ValueError("Use a new directory outside the extracted bundle")
    out.mkdir(parents=True, exist_ok=False)
    records = pd.read_csv(EXAMPLE / "predictions.csv", dtype={"participant_id": str}).merge(
        pd.read_csv(EXAMPLE / "order.csv"), on="sample_id", validate="many_to_one")
    weights = pd.read_csv(EXAMPLE / "weights.csv", dtype={"participant_id": str}).pivot(
        index="replicate_id", columns="participant_id", values="multiplicity")
    multiplicities = [pd.Series(1, index=weights.columns), *(row for _, row in weights.iterrows())]
    results = []
    for stem, (cap, n_specimens) in RULES.items():
        panel = records if cap is None else records.loc[records.retained_order <= cap]
        if panel.sample_id.nunique() != n_specimens:
            raise ValueError(f"Unexpected size for {stem}")
        declaration = json.loads((EXAMPLE / f"{stem}.json").read_text(encoding="utf-8"))
        formal = [formal_scores(panel, counts) for counts in multiplicities]
        independent = [independent_scores(panel, counts) for counts in multiplicities]
        formal_values = [{(obj, metric): scores[obj][metric] for obj in "ABCDE" for metric in METRICS}
                         for scores, _ in formal]
        independent_values = [scores for scores, _ in independent]
        score_rows = {key: [v[key] for v in formal_values] for key in formal_values[0]}
        pair_rows = {(a, b, name): [v[a, name] - v[b, name] for v in formal_values]
                     for a, b in CONTRASTS for name in METRICS}
        independent_score_rows = {(obj, metric): [v[obj, metric] for v in independent_values]
                                 for obj in "ABCDE" for metric in METRICS}
        independent_pair_rows = {(a, b, name): [v[a, name] - v[b, name] for v in independent_values]
                                 for a, b in CONTRASTS for name in METRICS}
        identity_error = float(max(residual for _, residual in formal))
        np.testing.assert_allclose(identity_error, 0, rtol=0, atol=ATOL)
        common = dict(direction=declaration["direction"], model=declaration["model"],
                      eligibility=declaration["eligibility"], n_specimens=n_specimens,
                      n_participants=4, n_runs=2, accepted_draws=3)
        rules = {rule["object"]: rule for rule in declaration["rules"]}
        score_output = []
        for (obj, metric), draws in score_rows.items():
            point, mean, low, high = summaries(np.asarray(draws, dtype=float))
            score_output.append(dict(**common, object=obj, evaluation_record=rules[obj]["evaluation_record"],
                                     target_weighting=rules[obj]["target_weighting"],
                                     within_person=rules[obj]["within_person"],
                                     across_runs=rules[obj]["across_runs"], metric=metric,
                                     point=point, bootstrap_mean=mean, ci_low=low, ci_high=high))
        cells = out / f"{stem}_scores.csv"
        pairs = out / f"{stem}_contrasts.csv"
        pd.DataFrame(score_output).to_csv(cells, index=False)
        contrast_output = []
        for (left, right, metric), draws in pair_rows.items():
            point, mean, low, high = summaries(np.asarray(draws, dtype=float))
            contrast_output.append(dict(**common, left_object=left, right_object=right,
                                       metric=metric, point=point, bootstrap_mean=mean,
                                       ci_low=low, ci_high=high))
        pd.DataFrame(contrast_output).to_csv(pairs, index=False)
        expected_dir = EXAMPLE / "expected"
        expected_score_error = compare_fixed(cells, expected_dir / f"{stem}_scores.csv", ["object", "metric"], n_specimens)
        expected_contrast_error = compare_fixed(pairs, expected_dir / f"{stem}_contrasts.csv",
                                                ["left_object", "right_object", "metric"], n_specimens)
        independent_score_error = compare_independent(cells, ["object", "metric"], independent_score_rows, n_specimens)
        independent_contrast_error = compare_independent(pairs, ["left_object", "right_object", "metric"], independent_pair_rows, n_specimens)
        results.append(dict(declaration=stem, specimens=n_specimens, scores=len(score_rows), contrasts=len(pair_rows),
                            max_expected_score_error=expected_score_error,
                            max_expected_contrast_error=expected_contrast_error,
                            max_independent_score_error=independent_score_error,
                            max_independent_contrast_error=independent_contrast_error,
                            max_brier_identity_error=identity_error))
    report = dict(status="PASS", inputs="explicit synthetic prediction, weight, order and declaration files",
                  expected_outputs="fixed files in expected/", independent_reference="grouped pandas calculation",
                  n_participants=4, n_runs=2, n_draws=3, absolute_tolerance=ATOL, checks=results)
    (out / "validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("PASS: fixed expected outputs and an independent reference calculation agree for 45 object scores and 63 paired contrasts.")


if __name__ == "__main__":
    main()
