from pathlib import Path
import importlib.util
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "replay_registered.py"
VALIDATION_SCRIPT = Path(__file__).resolve().parents[1] / "reproducibility" / "validation" / "replay_registered.py"


def load_replay_module():
    spec = importlib.util.spec_from_file_location("replay_registered", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def replay_inputs(tmp_path):
    # Three identical positive-sample scores must still tie a single negative
    # score after participant aggregation, even for a value such as 0.1.
    records = [("case", 1, 3), ("control_a", 0, 1), ("control_b", 0, 1)]
    predictions = [
        ("toy_transfer", split, "constant", f"{person}_{visit}", outcome, 0.1, person)
        for split in (1, 2)
        for person, outcome, count in records
        for visit in range(count)
    ]
    frame = pd.DataFrame(predictions, columns=[
        "direction", "split", "model", "sample_id", "outcome_ibd", "probability", "bootstrap_cluster",
    ])
    frame.to_csv(tmp_path / "predictions.csv", index=False)
    weights = pd.DataFrame([
        ("toy_transfer", draw, person, weight)
        for draw, counts in enumerate([(1, 1, 1), (2, 1, 0)])
        for (person, _, _), weight in zip(records, counts)
    ], columns=["direction", "replicate", "bootstrap_cluster", "weight"])
    weights.to_csv(tmp_path / "weights.csv", index=False)
    return tmp_path


def run_replay(directory):
    return subprocess.run([
        sys.executable, str(SCRIPT), "--predictions", str(directory / "predictions.csv"),
        "--weights", str(directory / "weights.csv"), "--output-dir", str(directory / "out"),
    ], capture_output=True, text=True)


def test_registered_cli_scores_ties_with_stored_multiplicities(replay_inputs):
    result = run_replay(replay_inputs)
    assert result.returncode == 0, result.stderr
    draw_rows = pd.read_csv(replay_inputs / "out" / "orthogonal_cell_draws.csv")
    scored = draw_rows.set_index(["cell", "metric", "replicate"]).value
    assert len(scored) == 30
    for draw, (sample_rate, person_rate) in enumerate([(3 / 5, 1 / 3), (6 / 7, 2 / 3)]):
        for cell, rate in [("A", sample_rate), ("B", sample_rate), ("C", person_rate),
                           ("D", person_rate), ("E", person_rate)]:
            assert scored[cell, "roc_auc", draw] == 0.5
            assert scored[cell, "average_precision", draw] == pytest.approx(rate)
            assert scored[cell, "brier", draw] == pytest.approx(rate * 0.9**2 + (1 - rate) * 0.1**2)


def test_registered_cli_rejects_incomplete_participant_draw(replay_inputs):
    path = replay_inputs / "weights.csv"
    rows = pd.read_csv(path)
    # Omit even the zero-weight participant: each draw must cover the full panel.
    rows.loc[~((rows.replicate == 1) & (rows.bootstrap_cluster == "control_b"))].to_csv(path, index=False)
    result = run_replay(replay_inputs)
    assert result.returncode != 0
    assert "full participant panel" in result.stderr
    assert not list((replay_inputs / "out").glob("*.csv"))


def test_equal_participant_weights_follow_interleaved_sample_indices():
    replay = load_replay_module()
    rows = [
        ("toy", "1", "model", "s1", 1, 0.90, "P1"),
        ("toy", "1", "model", "s2", 0, 0.80, "P2"),
        ("toy", "1", "model", "s3", 1, 0.60, "P1"),
        ("toy", "1", "model", "s4", 0, 0.20, "P3"),
        ("toy", "1", "model", "s5", 0, 0.10, "P3"),
    ]
    predictions = pd.DataFrame(rows, columns=[
        "direction", "split", "model", "sample_id", "outcome_ibd", "probability",
        "bootstrap_cluster",
    ])
    draws = pd.DataFrame(
        [("toy", 0, person, 1.0) for person in ("P1", "P2", "P3")],
        columns=["direction", "replicate", "bootstrap_cluster", "weight"],
    )
    samples, people, groups, stool_weights, person_weights, _ = replay.direction_inputs(
        predictions, draws
    )
    y = samples.outcome_ibd.to_numpy(dtype=int)
    labels = people.outcome_ibd.to_numpy(dtype=int)
    probability = predictions.set_index("sample_id").loc[samples.sample_id, "probability"].to_numpy()
    run_matrix = probability[None, :]
    counts = np.empty(len(y), dtype=float)
    for group in groups:
        counts[group] = len(group)
    expected_weights = stool_weights / counts[None, :]
    expected = replay.score(y, probability, expected_weights)[0]
    actual = replay.four_scores(
        run_matrix, y, labels, groups, stool_weights, person_weights
    )["E"][0]
    assert actual == pytest.approx(expected)
    assert all(np.array_equal(np.sum(expected_weights[:, group], axis=1), person_weights[:, i])
               for i, group in enumerate(groups))


def test_fractional_weight_auroc_is_not_integerized():
    replay = load_replay_module()
    labels = np.array([1, 1, 0, 0], dtype=int)
    probability = np.array([0.90, 0.60, 0.80, 0.10], dtype=float)
    weights = np.full((1, 4), 0.5, dtype=float)
    run_result = replay.score(labels, probability, weights)
    averaged = replay.average_scores([run_result], labels, weights)
    assert averaged[0, 0] == pytest.approx(0.75)


def test_validation_entry_accepts_current_public_schema(tmp_path):
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run([
        sys.executable, str(VALIDATION_SCRIPT),
        "--predictions", str(root / "reproducibility" / "examples" / "custom_panel" / "predictions.csv"),
        "--weights", str(root / "reproducibility" / "examples" / "custom_panel" / "weights.csv"),
        "--output-dir", str(tmp_path / "validation_replay"),
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "validation_replay" / "orthogonal_cell_summary.csv").exists()
