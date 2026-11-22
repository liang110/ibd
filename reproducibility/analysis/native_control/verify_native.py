"""Check recorded SIAMCAT readouts and the identities of all shared source inputs."""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def verify(root, output):
    sources = rows(output / "native_inputs.csv")
    splits = [f"analysis/import_case/native/split_{i:02d}.csv" for i in range(1, 11)]
    registry = "evidence/full_reproducibility_inputs/prediction_registry.csv"
    reference = "analysis/import_case/results/scores.csv"
    if [row["path"] for row in sources] != splits + [registry, reference]:
        raise ValueError("native input index mismatch")
    for row in sources:
        with (root / row["path"]).open("rb") as handle:
            digest = hashlib.sha256()
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        if digest.hexdigest() != row["sha256"]:
            raise ValueError(f"input SHA-256 mismatch: {row['path']}")
    panel = set()
    with (root / registry).open(encoding="utf-8-sig", newline="") as handle:
        for record in csv.DictReader(handle):
            if record["direction"] == "PRISM_to_HMP2":
                panel.add((record["sample_id"], record["bootstrap_cluster"], int(record["outcome_ibd"])))
    ids, people, labels = zip(*sorted(panel))
    if (len(ids), len(set(ids)), len(set(people)), sum(labels)) != (385, 385, 105, 281):
        raise ValueError("native target mismatch")
    if len(set(zip(people, labels))) != len(set(people)):
        raise ValueError("inconsistent labels within participants")
    runs = []
    for rel in splits:
        records = rows(root / rel)
        by_id = {r["sample_id"]: float(r["probability"]) for r in records}
        if len(records) != len(ids) or set(by_id) != set(ids):
            raise ValueError(f"unaligned sample IDs: {rel}")
        runs.append([by_id[sample] for sample in ids])
    runs = np.asarray(runs)
    if not np.isfinite(runs).all() or runs.min() < 0 or runs.max() > 1:
        raise ValueError("invalid native probability")
    individual = [roc_auc_score(labels, p) for p in runs]
    repeat_rows = rows(output / "native_run_auroc.csv")
    if [r["run"] for r in repeat_rows] != [Path(rel).stem for rel in splits]:
        raise ValueError("native AUROC run order mismatch")
    np.testing.assert_allclose(individual, [float(r["auroc"]) for r in repeat_rows], atol=1e-12, rtol=0)
    native = rows(output / "native_auroc.csv")
    if [r["object"] for r in native] != ["A", "B"] or [r["statistic"] for r in native] != ["mean(auroc.all)", "auroc"]:
        raise ValueError("native readout keys mismatch")
    protocol = {r["object"]: float(r["complete_target"]) for r in rows(root / reference) if r["metric"] == "roc_auc"}
    for r, independently_scored in zip(native, (np.mean(individual), roc_auc_score(labels, runs.mean(axis=0)))):
        if (r["n_specimens"], r["n_people"], r["n_runs"], r["positive_class"], float(r["tolerance"])) != (
            "385", "105", "10", "IBD", 1e-12
        ):
            raise ValueError("native readout panel or tolerance mismatch")
        native_value = float(r["native_auroc"])
        expected = protocol[r["object"]]
        np.testing.assert_allclose(native_value, independently_scored, atol=1e-12, rtol=0)
        np.testing.assert_allclose(native_value, expected, atol=1e-12, rtol=0)
        np.testing.assert_allclose(float(r["protocol_auroc"]), expected, atol=1e-15, rtol=0)
        np.testing.assert_allclose(float(r["absolute_difference"]), abs(native_value - expected), atol=1e-15, rtol=0)
    print("PASS: 12 input hashes, 10 individual AUROCs, 2 protocol-matched SIAMCAT readouts.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("native_results", type=Path)
    args = parser.parse_args()
    verify(args.bundle, args.native_results)
