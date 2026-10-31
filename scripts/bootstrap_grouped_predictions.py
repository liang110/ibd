#!/usr/bin/env python
"""Shared-target participant-cluster bootstrap for frozen predictions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


METRICS = ("roc_auc", "average_precision", "brier")


def _metric_values(frame: pd.DataFrame, weights: np.ndarray | None = None) -> dict[str, float]:
    return {
        "roc_auc": float(
            roc_auc_score(frame.outcome_ibd, frame.probability, sample_weight=weights)
        ),
        "average_precision": float(
            average_precision_score(frame.outcome_ibd, frame.probability, sample_weight=weights)
        ),
        "brier": float(
            brier_score_loss(frame.outcome_ibd, frame.probability, sample_weight=weights)
        ),
    }


def _shared_weight_map(
    frame: pd.DataFrame, rng: np.random.Generator
) -> dict[str, int]:
    """Draw one participant bootstrap vector for an aggregate target cell."""
    clusters = np.sort(frame.bootstrap_cluster.astype(str).unique())
    outcomes = frame.assign(bootstrap_cluster=frame.bootstrap_cluster.astype(str)).groupby(
        "bootstrap_cluster", sort=True
    ).outcome_ibd.first()
    for _ in range(100):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        counts = pd.Series(sampled).value_counts().to_dict()
        retained = outcomes.index.intersection(counts)
        if outcomes.loc[retained].nunique() == 2:
            return {str(cluster): int(counts.get(cluster, 0)) for cluster in clusters}
    raise RuntimeError("Unable to retain both outcomes in shared participant bootstrap")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--group-columns", required=True)
    parser.add_argument("--aggregate-columns", required=True)
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260818)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    group_columns = args.group_columns.split(",")
    aggregate_columns = args.aggregate_columns.split(",")
    predictions = pd.read_csv(args.predictions)
    required = {"outcome_ibd", "probability", "bootstrap_cluster", *group_columns}
    missing = required.difference(predictions.columns)
    if missing:
        raise ValueError(f"Prediction table missing columns: {sorted(missing)}")
    rng = np.random.default_rng(args.seed)
    summaries, draws = [], []
    shared_weights: dict[tuple, list[dict[str, int]]] = {}
    grouped_frames = {
        (keys,) if len(group_columns) == 1 else keys: frame.copy()
        for keys, frame in predictions.groupby(group_columns, sort=True)
    }
    aggregate_to_groups: dict[tuple, list[tuple]] = {}
    for group_key in grouped_frames:
        aggregate_key = tuple(group_key[group_columns.index(column)] for column in aggregate_columns)
        aggregate_to_groups.setdefault(aggregate_key, []).append(group_key)
    for aggregate_key, group_keys in aggregate_to_groups.items():
        reference = grouped_frames[group_keys[0]]
        shared_weights[aggregate_key] = [
            _shared_weight_map(reference, rng) for _ in range(args.n_bootstrap)
        ]
    for keys, frame in predictions.groupby(group_columns, sort=True):
        keys = (keys,) if len(group_columns) == 1 else keys
        key_values = dict(zip(group_columns, keys))
        aggregate_key = tuple(key_values[column] for column in aggregate_columns)
        observed = _metric_values(frame)
        boot = []
        for replicate in range(args.n_bootstrap):
            weight_map = shared_weights[aggregate_key][replicate]
            weights = frame.bootstrap_cluster.astype(str).map(weight_map).fillna(0.0).to_numpy(float)
            values = _metric_values(frame, weights)
            boot.append(values)
            draws.append({**key_values, "replicate": replicate, **values})
        boot_frame = pd.DataFrame(boot)
        summaries.append(
            {
                **key_values,
                "n_samples": len(frame),
                "n_clusters": int(frame.bootstrap_cluster.nunique()),
                **{f"{metric}_observed": observed[metric] for metric in METRICS},
                **{
                    f"{metric}_ci_low": float(boot_frame[metric].quantile(0.025))
                    for metric in METRICS
                },
                **{
                    f"{metric}_ci_high": float(boot_frame[metric].quantile(0.975))
                    for metric in METRICS
                },
            }
        )
    summary = pd.DataFrame(summaries)
    draw_frame = pd.DataFrame(draws)
    summary.to_csv(out / "cluster_bootstrap_by_split.csv", index=False)
    draw_frame.to_csv(out / "cluster_bootstrap_draws.csv", index=False)
    weight_rows = []
    for aggregate_key, vectors in shared_weights.items():
        aggregate_values = dict(zip(aggregate_columns, aggregate_key))
        for replicate, vector in enumerate(vectors):
            for cluster, multiplicity in vector.items():
                weight_rows.append(
                    {**aggregate_values, "replicate": replicate, "bootstrap_cluster": cluster,
                     "multiplicity": multiplicity}
                )
    pd.DataFrame(weight_rows).to_csv(out / "shared_target_bootstrap_weights.csv", index=False)
    aggregate_draws = draw_frame.groupby([*aggregate_columns, "replicate"], as_index=False)[
        list(METRICS)
    ].mean()
    aggregate_rows = []
    for keys, frame in aggregate_draws.groupby(aggregate_columns, sort=True):
        keys = (keys,) if len(aggregate_columns) == 1 else keys
        row = dict(zip(aggregate_columns, keys))
        for metric in METRICS:
            row[f"{metric}_bootstrap_mean"] = float(frame[metric].mean())
            row[f"{metric}_ci_low"] = float(frame[metric].quantile(0.025))
            row[f"{metric}_ci_high"] = float(frame[metric].quantile(0.975))
        aggregate_rows.append(row)
    aggregate = pd.DataFrame(aggregate_rows)
    aggregate.to_csv(out / "aggregate_cluster_bootstrap.csv", index=False)
    report = {
        "n_bootstrap": args.n_bootstrap,
        "seed": args.seed,
        "group_columns": group_columns,
        "aggregate_columns": aggregate_columns,
        "n_groups": int(len(summary)),
        "n_aggregate_cells": int(len(aggregate)),
        "weight_protocol": "one_shared_participant_weight_vector_per_aggregate_cell_and_replicate",
        "weight_manifest": "shared_target_bootstrap_weights.csv",
    }
    (out / "bootstrap_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
