"""Reproduce the registered S17 slice, its independent reference and all input checks."""
import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

CASE = Path(__file__).resolve().parent
DEFAULT_REGISTRY = CASE.parent / "evidence/full_reproducibility_inputs/prediction_registry.csv"
TOLERANCE = 1e-12

def records(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f)

def check_registration(registry, fixture, selection):
    expected, found = {}, {}
    for r in records(fixture):
        key = (r["run_id"], r["sample_id"])
        if key in expected:
            raise ValueError("Duplicate prediction key in fixture")
        expected[key] = (r["participant_id"], int(r["y_true"]), float(r["y_prob"]), r["origin"])
    for r in records(registry):
        if r["direction"] == selection["direction"] and r["model"] == selection["model"]:
            key = (r["split"], r["sample_id"])
            if key in found:
                raise ValueError("Duplicate prediction key in registry")
            found[key] = (r["bootstrap_cluster"], int(r["outcome_ibd"]), float(r["probability"]), r["origin"])
    if found != expected or len(expected) != selection["n_runs"] * selection["n_samples"]:
        raise ValueError("S17 rows differ from the selected registry slice")
    return len(found)

def scores(path):
    entries = list(records(path))
    values = {(r["object"], r["metric"]): float(r["value"]) for r in entries}
    expected = {(obj, metric) for obj in "ABCD" for metric in ("AUROC", "AP", "Brier")}
    if len(entries) != 12 or set(values) != expected or not all(map(math.isfinite, values.values())):
        raise ValueError("Incomplete four-object results")
    return values

def run(name, *args):
    subprocess.run([sys.executable, "-B", str(CASE / name), *map(str, args)], check=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output == CASE or CASE in output.parents:
        raise SystemExit("Choose a replay output directory outside the archived protocol-control folder")
    fixture = CASE / "input/predictions_long.csv"
    spec = CASE / "input/scoring_specification.json"
    selection = json.loads(spec.read_text(encoding="utf-8-sig"))["input_selection"]
    n = check_registration(args.registry, fixture, selection)
    common = ("--predictions", fixture, "--weights", CASE / "input/shared_participant_multiplicities.csv", "--spec", spec)
    run("run_protocol_control.py", *common, "--output-dir", output)
    run("reference_impl.py", "--predictions", fixture, "--output", output / "reference_outputs.csv")
    tests = {}
    for case_file, report in (("cases.json", "fail_closed_report.csv"), ("interface_cases.json", "interface_checks.csv")):
        run("validate_fail_closed.py", *common, "--cases", CASE / "invalid_cases" / case_file, "--report", output / report)
        tests[case_file] = len(list(records(output / report)))
    calculated = scores(output / "protocol_control_outputs.csv")
    independent = scores(output / "reference_outputs.csv")
    archived = scores(CASE / "output/protocol_control_outputs.csv")
    differences = [{"object": k[0], "metric": k[1], "value_canonical": a, "value_reference": independent[k],
                    "abs_diff": abs(a - independent[k])} for k, a in calculated.items()]
    reference_error = max(row["abs_diff"] for row in differences)
    archive_error = max(abs(calculated[k] - archived[k]) for k in calculated)
    if max(reference_error, archive_error) > TOLERANCE:
        raise SystemExit("Numeric replay exceeds tolerance")
    with (output / "independent_comparison.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(differences[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(differences)
    report = {"registry_slice_rows": n, "values_checked": len(calculated),
              "max_abs_difference_to_reference": reference_error, "max_abs_difference_to_archived": archive_error,
              "tolerance": TOLERANCE, "rejected_cases": tests, "result": "PASS"}
    (output / "reader_check.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))

if __name__ == "__main__":
    main()
