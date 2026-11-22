"""Score the declared S17 full-target example. Bootstrap weights use a separate entry."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

METRICS = ["AUROC", "AP", "Brier"]
PRED_COLUMNS = {"direction", "model", "run_id", "sample_id", "participant_id", "y_true", "y_prob"}
WEIGHT_COLUMNS = {"direction", "replicate_id", "participant_id", "multiplicity"}
OBJECTS = {
    "A": ("sample", "none", "mean_over_runs"),
    "B": ("sample", "mean_over_runs", "single_metric"),
    "C": ("participant", "none", "mean_over_runs"),
    "D": ("participant", "mean_over_runs", "single_metric"),
}

def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def validate(pred, weights, spec):
    if spec.get("schema_version") != "protocol-control-1.0":
        raise ValueError("unsupported_scoring_schema")
    if "objects" not in spec:
        raise ValueError("missing_scoring_declaration")
    fields = ("scoring_unit", "probability_aggregation", "metric_aggregation")
    expected_objects = {obj: dict(zip(fields, rules)) for obj, rules in OBJECTS.items()}
    if spec["objects"] != expected_objects or spec.get("metrics") != METRICS:
        raise ValueError("unsupported_scoring_declaration")
    if not str(spec.get("uncertainty_target", "")).startswith("point_estimate_only;"):
        raise ValueError("expected_point_estimate_only")
    required = set(spec.get("required_prediction_columns", []))
    weight_cols = set(spec.get("required_multiplicity_columns", []))
    if not PRED_COLUMNS.issubset(required) or not WEIGHT_COLUMNS.issubset(weight_cols):
        raise ValueError("incomplete_schema_declaration")
    if not required.issubset(pred):
        raise ValueError("missing_prediction_columns")
    if not weight_cols.issubset(weights):
        raise ValueError("missing_multiplicity_columns")
    if pred.empty or pred[list(required)].isna().any().any():
        raise ValueError("missing_prediction_values")
    if weights.empty or weights[list(weight_cols)].isna().any().any():
        raise ValueError("missing_multiplicity_values")
    selection = spec.get("input_selection", {})
    if not pred.direction.eq(selection.get("direction")).all() or not pred.model.eq(selection.get("model")).all():
        raise ValueError("prediction_selection_mismatch")
    if not weights.direction.eq(selection["direction"]).all():
        raise ValueError("weight_direction_mismatch")
    if pred.duplicated(["direction", "model", "run_id", "sample_id"]).any():
        raise ValueError("duplicate_prediction_key")
    if not pred.y_true.isin([0, 1]).all():
        raise ValueError("nonbinary_target")
    if pred.y_true.nunique() != 2:
        raise ValueError("single_class_draw")
    if not np.isfinite(pred.y_prob).all() or not pred.y_prob.between(0, 1).all():
        raise ValueError("probability_out_of_range")
    if pred.groupby("participant_id").y_true.nunique().gt(1).any():
        raise ValueError("participant_label_conflict")
    per_sample = pred.groupby("sample_id")
    if per_sample.participant_id.nunique().ne(1).any():
        raise ValueError("sample_participant_conflict")
    runs = pred.run_id.nunique()
    if runs != selection.get("n_runs") or per_sample.run_id.nunique().ne(runs).any():
        raise ValueError("incomplete_run_coverage")
    if len(per_sample) != selection.get("n_samples") or pred.participant_id.nunique() != selection.get("n_participants"):
        raise ValueError("selection_size_mismatch")
    if weights.duplicated(["direction", "replicate_id", "participant_id"]).any():
        raise ValueError("duplicate_multiplicity_key")
    if weights.replicate_id.nunique() != 1:
        raise ValueError("expected_one_point_fixture")
    if set(weights.participant_id) != set(pred.participant_id):
        raise ValueError("multiplicity_participant_mismatch")
    if not weights.multiplicity.eq(1).all():
        raise ValueError("expected_all_one_point_fixture")

def metric_dict(frame):
    labels, scores = frame.y_true, frame.y_prob
    return dict(zip(METRICS, (roc_auc_score(labels, scores), average_precision_score(labels, scores),
                              brier_score_loss(labels, scores))))

def aggregate(frame, keys):
    group = frame.groupby(keys, as_index=False, sort=True).y_prob
    summary = group.agg(["min", "max", "mean"])
    # Preserve an exact tie when a participant or sample contributes identical scores.
    summary["y_prob"] = summary["mean"].where(summary["min"].ne(summary["max"]), summary["min"])
    return summary[keys + ["y_prob"]]

def four_objects(pred):
    per_run, per_person_run = [], []
    for _, group in pred.groupby("run_id", sort=True):
        per_run.append(metric_dict(group))
        per_person_run.append(metric_dict(aggregate(group, ["participant_id", "y_true"])))
    ensemble = aggregate(pred, ["sample_id", "participant_id", "y_true"])
    return {"A": pd.DataFrame(per_run).mean().to_dict(), "B": metric_dict(ensemble),
            "C": pd.DataFrame(per_person_run).mean().to_dict(),
            "D": metric_dict(aggregate(ensemble, ["participant_id", "y_true"]))}

def read_predictions(path):
    return pd.read_csv(path, dtype={k: "string" for k in ("direction", "model", "run_id", "participant_id", "sample_id")},
                       float_precision="round_trip")

def read_weights(path):
    return pd.read_csv(path, dtype={k: "string" for k in ("direction", "participant_id", "replicate_id")})

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("predictions", "weights", "spec", "output-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    pred = read_predictions(args.predictions)
    weights = read_weights(args.weights)
    spec = json.loads(args.spec.read_text(encoding="utf-8-sig"))
    validate(pred, weights, spec)
    values = four_objects(pred)
    meta = {"direction": str(pred.direction.iloc[0]), "model": str(pred.model.iloc[0]),
            "n_runs": int(pred.run_id.nunique()), "n_samples": int(pred.sample_id.nunique()),
            "n_participants": int(pred.participant_id.nunique()), "input_sha256": sha256_file(args.predictions),
            "weights_sha256": sha256_file(args.weights), "spec_sha256": sha256_file(args.spec),
            "uncertainty_target": spec["uncertainty_target"]}
    rows = [{"object": obj, "metric": metric, "value": values[obj][metric], **spec["objects"][obj], **meta}
            for obj in OBJECTS for metric in METRICS]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "protocol_control_outputs.csv", index=False, lineterminator="\n")
    (args.output_dir / "protocol_control_summary.json").write_text(
        json.dumps({"metadata": meta, "objects": values}, indent=2) + "\n", encoding="utf-8")
    print("S17 full target PASS: 12 values")

if __name__ == "__main__":
    main()
