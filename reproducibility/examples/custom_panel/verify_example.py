"""Run the supplied declarations and independently score their small synthetic panel."""
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


def reference_scores(panel, counts):
    """Reconstruct target weights and participant means from the selected raw rows."""
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
            scores[code, name] = np.mean([metric(labels, prediction, sample_weight=weight) for prediction in runs])
    residual = abs(scores["E", "brier"] - scores["C", "brier"] - within_variance)
    return scores, residual


def summaries(values):
    draws = values[1:]
    return [values[0], np.mean(draws), *np.quantile(draws, [0.025, 0.975], method="linear")]


def compare(file, keys, values, n_specimens):
    actual = pd.read_csv(file, float_precision="round_trip").set_index(keys)
    if not actual.index.is_unique or set(actual.index) != set(values):
        raise ValueError(f"Incomplete result keys in {file.name}")
    for field, count in {"n_runs": 2, "n_specimens": n_specimens, "n_participants": 4, "accepted_draws": 3}.items():
        if not actual[field].eq(count).all():
            raise ValueError(f"Incorrect {field} in {file.name}")
    errors = []
    for key, reference in values.items():
        estimate = actual.loc[key, FIELDS].to_numpy(dtype=float)
        expected = summaries(reference)
        np.testing.assert_allclose(estimate, expected, rtol=0, atol=ATOL)
        errors.append(np.max(np.abs(estimate - expected)))
    return float(max(errors))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New directory outside the extracted bundle")
    out = parser.parse_args().output.resolve()
    if out.is_relative_to(BUNDLE):
        raise ValueError("Use a new directory outside the extracted bundle")
    out.mkdir(parents=True, exist_ok=False)
    records = pd.read_csv(EXAMPLE/"predictions.csv", dtype={"participant_id": str}).merge(
        pd.read_csv(EXAMPLE/"order.csv"), on="sample_id", validate="many_to_one")
    weights = pd.read_csv(EXAMPLE/"weights.csv", dtype={"participant_id": str}).pivot(
        index="replicate_id", columns="participant_id", values="multiplicity")
    multiplicities = [pd.Series(1, index=weights.columns), *(row for _, row in weights.iterrows())]
    results = []
    for stem, (cap, n_specimens) in RULES.items():
        panel = records if cap is None else records.loc[records.retained_order <= cap]
        if panel.sample_id.nunique() != n_specimens:
            raise ValueError(f"Unexpected size for {stem}")
        cells = out/f"{stem}_scores.csv"
        pairs = out/f"{stem}_contrasts.csv"
        declaration = json.loads((EXAMPLE/f"{stem}.json").read_text(encoding="utf-8"))
        references = [reference_scores(panel, counts) for counts in multiplicities]
        values = [scores for scores, _ in references]
        score_rows = {key: [v[key] for v in values] for key in values[0]}
        pair_rows = {(a, b, name): [v[a, name] - v[b, name] for v in values]
                     for a, b in CONTRASTS for name in METRICS}
        identity_error = float(max(residual for _, residual in references))
        np.testing.assert_allclose(identity_error, 0, rtol=0, atol=ATOL)

        # Write the same result-level schema used by the public evaluator.
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
        pd.DataFrame(score_output).to_csv(cells, index=False)
        contrast_output = []
        for (left, right, metric), draws in pair_rows.items():
            point, mean, low, high = summaries(np.asarray(draws, dtype=float))
            contrast_output.append(dict(**common, left_object=left, right_object=right,
                                       metric=metric, point=point, bootstrap_mean=mean,
                                       ci_low=low, ci_high=high))
        pd.DataFrame(contrast_output).to_csv(pairs, index=False)
        results.append(dict(declaration=stem, specimens=n_specimens, scores=len(score_rows), contrasts=len(pair_rows),
                            max_score_error=compare(cells, ["object", "metric"], score_rows, n_specimens),
                            max_contrast_error=compare(pairs, ["left_object", "right_object", "metric"], pair_rows, n_specimens),
                            max_brier_identity_error=identity_error))
    report = dict(status="PASS", inputs="explicit synthetic prediction, weight, order and declaration files",
                  n_participants=4, n_runs=2, n_draws=3, absolute_tolerance=ATOL, checks=results)
    (out/"validation.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print("PASS: 45 object scores and 63 paired contrasts checked against the independent calculation.")


if __name__ == "__main__":
    main()
