"""Run the named input mutations through the production S17 validation routine."""
import argparse
import copy
import json
from pathlib import Path

import pandas as pd
from run_protocol_control import read_predictions, read_weights, validate

def mutate(name, pred, weights, spec):
    i = pred.index[0]
    wi = weights.index[0]
    if name == "duplicate_key":
        pred = pd.concat([pred, pred.iloc[[0]]], ignore_index=True)
    elif name == "label_conflict":
        pred.loc[i, "y_true"] = 1 - int(pred.loc[i, "y_true"])
    elif name == "probability_out_of_range":
        pred.loc[i, "y_prob"] = 1.2
    elif name == "incomplete_runs":
        pred = pred[pred.run_id != pred.loc[i, "run_id"]]
    elif name == "weight_mismatch":
        weights = weights.iloc[:-1]
    elif name == "missing_declaration":
        del spec["objects"]
    elif name == "single_class":
        pred["y_true"] = 0
    elif name == "nonuniform_weights":
        weights.loc[wi, "multiplicity"] = 0
        weights.loc[weights.index[1], "multiplicity"] = 2
    elif name == "duplicate_weight":
        weights = pd.concat([weights, weights.iloc[[0]]], ignore_index=True)
    elif name == "direction_mismatch":
        weights.loc[wi, "direction"] = "unexpected_direction"
    elif name == "multiple_draws":
        weights.loc[wi, "replicate_id"] = "another_rep"
    elif name == "model_mismatch":
        pred.loc[i, "model"] = "unexpected_model"
    elif name == "sample_owner_mismatch":
        alternative = pred[(pred.y_true == pred.loc[i, "y_true"]) & (pred.participant_id != pred.loc[i, "participant_id"])]
        pred.loc[i, "participant_id"] = alternative.participant_id.iloc[0]
    elif name == "fractional_label":
        pred["y_true"] = pred.y_true.astype(float)
        pred.loc[i, "y_true"] = 0.5
    elif name == "missing_sample_run":
        pred = pred.drop(index=i)
    elif name == "objects_swapped":
        spec["objects"]["B"] = spec["objects"]["A"]
    elif name == "changed_count":
        spec["input_selection"]["n_samples"] += 1
    else:
        raise ValueError("Unknown test mutation: " + name)
    return pred, weights, spec

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("predictions", "weights", "spec", "cases", "report"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    pred = read_predictions(args.predictions)
    weights = read_weights(args.weights)
    spec = json.loads(args.spec.read_text(encoding="utf-8-sig"))
    validate(pred, weights, spec)
    configs = json.loads(args.cases.read_text(encoding="utf-8-sig"))
    results = []
    for case, cfg in configs.items():
        trial = mutate(cfg["mutation"], pred.copy(), weights.copy(), copy.deepcopy(spec))
        reason = "unexpected_pass"
        try:
            validate(*trial)
        except ValueError as error:
            reason = str(error)
        rejected = reason != "unexpected_pass"
        results.append({"case": case, "expected": cfg["expected"], "observed": "FAIL" if rejected else "PASS",
                        "reason": reason, "pass": cfg["expected"] == "FAIL" and reason == cfg.get("reason", case)})
    report = pd.DataFrame(results)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.report, index=False, lineterminator="\n")
    if not report["pass"].all():
        raise SystemExit("Validation disagreement. See " + str(args.report))
    print("Rejected all", len(report), "invalid cases with their expected diagnostic")

if __name__ == "__main__":
    main()
