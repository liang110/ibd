from __future__ import annotations

import pandas as pd
import pytest

from ibd_microbiome_framework.data.contracts import (
    MatrixSpec,
    build_manifest,
    validate_participant_labels,
    write_contracts,
)


def _write_inputs(root, prefix: str, feature_column: str):
    matrix = pd.DataFrame({"S1": [1.0, 2.0], "S2": [2.0, 1.0]}, index=["taxon_a", "taxon_b"])
    matrix.index.name = feature_column
    matrix_path = root / f"{prefix}_abundance.csv"
    matrix.to_csv(matrix_path)
    metadata = pd.DataFrame(
        {"sample_id": ["S1", "S2"], "group": ["Control", "IBD"], "Participant ID": ["P1", "P2"]}
    )
    metadata_path = root / f"{prefix}_metadata.csv"
    metadata.to_csv(metadata_path, index=False)
    return matrix_path, metadata_path


def test_contract_summary_supports_microbiome_only_and_optional_metabolome(tmp_path):
    hmp_matrix, hmp_metadata = _write_inputs(tmp_path, "hmp", "Feature")
    prism_matrix, prism_metadata = _write_inputs(tmp_path, "prism", "Species")
    specs = [
        MatrixSpec("HMP2", "microbiome", hmp_matrix, "Feature", hmp_metadata),
        MatrixSpec("PRISM", "microbiome", prism_matrix, "Species", prism_metadata),
    ]
    summary = write_contracts(specs, tmp_path / "contract", tmp_path)
    assert summary["n_source_matrices"] == 2
    assert summary["prism_manifest_origins"] == {"unknown": 2}


def _spec(tmp_path, cohort, group="Control"):
    metadata = tmp_path / f"{cohort.lower()}_metadata.csv"
    pd.DataFrame(
        {"sample_id": ["s1"], "group": [group], "Participant ID": ["p1"]}
    ).to_csv(metadata, index=False)
    return MatrixSpec(
        cohort=cohort,
        modality="microbiome",
        path=tmp_path / "matrix.csv",
        feature_column="feature",
        metadata_path=metadata,
    )


@pytest.mark.parametrize(
    "cohort,group,expected",
    [("HMP2", "Control", 0), ("HMP2", "CD", 1), ("PRISM", "UC", 1),
     ("NielsenHB_2014", "healthy", 0), ("NielsenHB_2014", "IBD", 1)],
)
def test_explicit_label_allow_list(tmp_path, cohort, group, expected):
    manifest = build_manifest(_spec(tmp_path, cohort, group), ["s1"])
    assert manifest.loc[0, "outcome_ibd"] == expected


@pytest.mark.parametrize("bad_label", ["", "unknown", "case"])
def test_unknown_or_missing_labels_fail_closed(tmp_path, bad_label):
    with pytest.raises(ValueError, match="label mapping is not closed"):
        build_manifest(_spec(tmp_path, "HMP2", bad_label), ["s1"])


def test_missing_labels_fail_closed(tmp_path):
    spec = _spec(tmp_path, "HMP2")
    frame = pd.read_csv(spec.metadata_path)
    frame.loc[0, "group"] = pd.NA
    frame.to_csv(spec.metadata_path, index=False)
    with pytest.raises(ValueError, match="label mapping is not closed"):
        build_manifest(spec, ["s1"])


def test_manifest_rejects_conflicting_participant_labels(tmp_path):
    matrix = tmp_path / "matrix.csv"
    matrix.write_text(",s1,s2\nf1,1,2\n", encoding="utf-8")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "sample_id,group,Participant ID\n"
        "s1,Control,p1\n"
        "s2,IBD,p1\n",
        encoding="utf-8",
    )
    spec = MatrixSpec("HMP2", "microbiome", matrix, "feature", metadata)
    with pytest.raises(ValueError, match="Conflicting participant-level labels"):
        build_manifest(spec, ["s1", "s2"])


def test_validate_participant_labels_rejects_conflicts():
    frame = pd.DataFrame(
        {"participant_id": ["p1", "p1"], "outcome_ibd": [0, 1]}
    )
    with pytest.raises(ValueError, match="Conflicting participant-level labels"):
        validate_participant_labels(frame, context="fixture")
