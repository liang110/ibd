"""Independent, standard-library implementation of the S17 complete-target point estimates."""
import argparse
import csv
from collections import defaultdict
from itertools import groupby
from math import fsum, isfinite
from pathlib import Path

def mean(values):
    # A constant score is an exact tie, independent of the cluster size.
    if min(values) == max(values):
        return values[0]
    return fsum(values) / len(values)

def metrics(entries):
    ordered = sorted(entries, key=lambda row: row[2], reverse=True)
    positives = sum(y for _, y, _ in ordered)
    negatives = len(ordered) - positives
    if not positives or not negatives:
        raise ValueError("Both endpoint classes must occur")
    seen_pos = seen_neg = 0
    concordance, ap_terms = [], []
    for _, tied in groupby(ordered, key=lambda row: row[2]):
        ys = [y for _, y, _ in tied]
        npos, nneg = sum(ys), len(ys) - sum(ys)
        seen_pos += npos
        seen_neg += nneg
        concordance.append(npos * (negatives - seen_neg + nneg / 2))
        ap_terms.append(npos * seen_pos / (seen_pos + seen_neg))
    return {"AUROC": fsum(concordance) / (positives * negatives), "AP": fsum(ap_terms) / positives,
            "Brier": mean([(p - y) ** 2 for _, y, p in ordered])}

def people(panel):
    scores = defaultdict(list)
    for person, y, p in panel:
        scores[person, y].append(p)
    return [(person, y, mean(vals)) for (person, y), vals in scores.items()]

def calculate(path):
    samples, person_labels, runs = {}, {}, defaultdict(dict)
    with path.open(encoding="utf-8-sig", newline="") as f:
        selection = set()
        for row in csv.DictReader(f):
            s, person = row["sample_id"], row["participant_id"]
            p, y = float(row["y_prob"]), int(row["y_true"])
            if not isfinite(p) or not 0 <= p <= 1 or y not in (0, 1):
                raise ValueError("Invalid probability or label")
            selection.add((row["direction"], row["model"]))
            if samples.setdefault(s, (person, y)) != (person, y) or person_labels.setdefault(person, y) != y:
                raise ValueError("Inconsistent sample or participant mapping")
            panel = runs[row["run_id"]]
            if s in panel:
                raise ValueError("Duplicate sample in run")
            panel[s] = (person, y, p)
    if len(selection) != 1 or any(panel.keys() != samples.keys() for panel in runs.values()):
        raise ValueError("Expected a complete, single-model and single-direction panel")
    sample_runs = [metrics(panel.values()) for panel in runs.values()]
    person_runs = [metrics(people(panel.values())) for panel in runs.values()]
    ensemble = [(person, y, mean([panel[s][2] for panel in runs.values()])) for s, (person, y) in samples.items()]
    def combine(scores):
        return {m: mean([run[m] for run in scores]) for m in scores[0]}
    return {"A": combine(sample_runs), "B": metrics(ensemble), "C": combine(person_runs),
            "D": metrics(people(ensemble))}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    values = calculate(args.predictions)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["object", "metric", "value"])
        for obj, scores in values.items():
            for metric, value in scores.items():
                writer.writerow([obj, metric, value])
    print("S17 independent reference PASS: 12 values")

if __name__ == "__main__":
    main()
