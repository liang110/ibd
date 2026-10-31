from __future__ import annotations

import pandas as pd

from ibd_microbiome_framework.data.contracts import MatrixSpec, write_contracts


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
