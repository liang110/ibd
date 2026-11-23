from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "replay_registered.py"


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
