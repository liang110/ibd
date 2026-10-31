from __future__ import annotations

import pandas as pd

from scripts.bootstrap_grouped_predictions import main


def test_bootstrap_writes_shared_target_weights(tmp_path, monkeypatch):
    predictions = pd.DataFrame(
        [
            {"direction": "A_to_B", "model": "m1", "split": 0, "sample_id": "s1", "outcome_ibd": 0, "probability": 0.1, "bootstrap_cluster": "p1"},
            {"direction": "A_to_B", "model": "m1", "split": 0, "sample_id": "s2", "outcome_ibd": 1, "probability": 0.9, "bootstrap_cluster": "p2"},
            {"direction": "A_to_B", "model": "m2", "split": 0, "sample_id": "s1", "outcome_ibd": 0, "probability": 0.2, "bootstrap_cluster": "p1"},
            {"direction": "A_to_B", "model": "m2", "split": 0, "sample_id": "s2", "outcome_ibd": 1, "probability": 0.8, "bootstrap_cluster": "p2"},
            {"direction": "A_to_B", "model": "m1", "split": 1, "sample_id": "s1", "outcome_ibd": 0, "probability": 0.3, "bootstrap_cluster": "p1"},
            {"direction": "A_to_B", "model": "m1", "split": 1, "sample_id": "s2", "outcome_ibd": 1, "probability": 0.7, "bootstrap_cluster": "p2"},
        ]
    )
    prediction_path = tmp_path / "predictions.csv"
    prediction_path.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(prediction_path, index=False)
    output_dir = tmp_path / "out"
    monkeypatch.setattr(
        "sys.argv",
        [
            "bootstrap_grouped_predictions.py",
            "--predictions",
            str(prediction_path),
            "--output-dir",
            str(output_dir),
            "--group-columns",
            "direction,model,split",
            "--aggregate-columns",
            "direction",
            "--n-bootstrap",
            "5",
            "--seed",
            "7",
        ],
    )
    main()
    weights = pd.read_csv(output_dir / "shared_target_bootstrap_weights.csv")
    assert set(weights.columns) == {"direction", "replicate", "bootstrap_cluster", "multiplicity"}
    vectors = weights.pivot_table(
        index=["direction", "replicate"], columns="bootstrap_cluster", values="multiplicity"
    )
    assert vectors.shape == (5, 2)
    assert (vectors.sum(axis=1) == 2).all()
