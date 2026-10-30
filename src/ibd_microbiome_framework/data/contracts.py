"""Build auditable cohort and feature contracts from frozen delimited matrices."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Sequence

import numpy as np
import pandas as pd

from ibd_microbiome_framework.harmonization.features import (
    canonical_metabolite,
    canonical_taxon,
    is_forbidden_feature,
)


@dataclass(frozen=True)
class MatrixSpec:
    cohort: str
    modality: str
    path: Path
    feature_column: str
    metadata_path: Path


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def matrix_axes(path: Path) -> tuple[List[str], List[str]]:
    """Read header and first-column IDs without loading the numeric matrix."""
    samples: List[str] = []
    features: List[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        samples = [cell.strip() for cell in header[1:] if cell.strip()]
        for row in reader:
            if row and row[0].strip():
                features.append(row[0].strip())
    return samples, features


def abundance_library_sizes(path: Path) -> pd.Series:
    """Compute sample sums for a modest abundance matrix using numeric coercion."""
    frame = pd.read_csv(path, index_col=0, low_memory=False)
    numeric = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return numeric.sum(axis=0)


def infer_prism_origin(sample_id: str) -> str:
    if sample_id.startswith("Validation|"):
        return "validation"
    if sample_id.startswith("PRISM|"):
        return "discovery"
    return "unknown"


def build_manifest(spec: MatrixSpec, matrix_samples: Sequence[str]) -> pd.DataFrame:
    metadata = pd.read_csv(spec.metadata_path)
    if "sample_id" not in metadata.columns:
        raise ValueError(f"{spec.metadata_path} lacks sample_id")
    manifest = metadata.copy()
    manifest = manifest[manifest["sample_id"].isin(matrix_samples)].copy()
    manifest["cohort"] = spec.cohort
    if spec.cohort == "PRISM":
        manifest["origin"] = manifest["sample_id"].map(infer_prism_origin)
        manifest["participant_id"] = manifest["sample_id"]
        manifest["site"] = "unknown"
        manifest["timepoint"] = "cross_sectional"
    else:
        manifest["origin"] = "HMP2"
        manifest["participant_id"] = manifest.get("Participant ID", manifest["sample_id"])
        manifest["site"] = manifest.get("site_name", "unknown")
        manifest["timepoint"] = manifest.get("week_num", np.nan)
    manifest["outcome_ibd"] = (~manifest["group"].isin(["Control", "nonIBD"])).astype(int)
    manifest["diagnosis_subtype"] = manifest.get("diagnosis", manifest["group"])
    if manifest["sample_id"].duplicated().any():
        raise ValueError(f"Duplicate sample IDs in {spec.cohort} manifest")
    return manifest.sort_values("sample_id").reset_index(drop=True)


def build_feature_registry(
    specs: Sequence[MatrixSpec],
    feature_lists: Mapping[str, Sequence[str]],
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for spec in specs:
        canonicalizer: Callable[[object], object]
        canonicalizer = canonical_taxon if spec.modality == "microbiome" else canonical_metabolite
        for raw_feature in feature_lists[spec.cohort + ":" + spec.modality]:
            canonical = canonicalizer(raw_feature)
            rows.append(
                {
                    "cohort": spec.cohort,
                    "modality": spec.modality,
                    "source_feature_id": raw_feature,
                    "canonical_feature_id": canonical,
                    "is_forbidden_clinical": is_forbidden_feature(raw_feature),
                    "mapping_status": "exact_normalized" if canonical else "excluded_or_unresolved",
                }
            )
    registry = pd.DataFrame(rows)
    valid = registry.dropna(subset=["canonical_feature_id"])
    presence = (
        valid.groupby(["modality", "canonical_feature_id"])["cohort"]
        .agg(lambda values: ";".join(sorted(set(values))))
        .rename("present_in")
        .reset_index()
    )
    presence["shared_across_hmp2_prism"] = presence["present_in"].eq("HMP2;PRISM")
    return registry.merge(presence, on=["modality", "canonical_feature_id"], how="left")


def write_contracts(
    specs: Sequence[MatrixSpec], output_dir: Path, source_root: Path
) -> Dict[str, object]:
    """Write manifests and a feature registry; never alter source matrices."""
    output_dir.mkdir(parents=True, exist_ok=True)
    axes: Dict[str, tuple[List[str], List[str]]] = {}
    feature_lists: Dict[str, Sequence[str]] = {}
    source_inventory: List[Dict[str, object]] = []
    for spec in specs:
        samples, features = matrix_axes(spec.path)
        key = spec.cohort + ":" + spec.modality
        axes[key] = (samples, features)
        feature_lists[key] = features
        source_inventory.append(
            {
                "cohort": spec.cohort,
                "modality": spec.modality,
                "path": str(spec.path),
                "metadata_path": str(spec.metadata_path),
                "sha256": sha256_file(spec.path),
                "n_samples_header": len(samples),
                "n_features_raw": len(features),
            }
        )
        manifest = build_manifest(spec, samples)
        manifest.to_csv(output_dir / f"cohort_manifest_{spec.cohort.lower()}.csv", index=False)

    registry = build_feature_registry(specs, feature_lists)
    registry.to_csv(output_dir / "feature_registry.csv", index=False)
    pd.DataFrame(source_inventory).to_csv(output_dir / "source_inventory.csv", index=False)

    library_rows: List[Dict[str, object]] = []
    for spec in specs:
        if spec.modality != "microbiome":
            continue
        sums = abundance_library_sizes(spec.path)
        for sample_id, value in sums.items():
            library_rows.append(
                {
                    "cohort": spec.cohort,
                    "sample_id": sample_id,
                    "library_sum": float(value),
                }
            )
    library = pd.DataFrame(library_rows)
    library.to_csv(output_dir / "microbiome_library_sums.csv", index=False)

    shared = registry.loc[
        registry["shared_across_hmp2_prism"].fillna(False) & ~registry["is_forbidden_clinical"],
        ["modality", "canonical_feature_id"],
    ].drop_duplicates()
    shared.to_csv(output_dir / "shared_feature_candidates.csv", index=False)

    summary = {
        "source_root": str(source_root),
        "n_source_matrices": len(specs),
        "shared_microbiome_candidates": int((shared["modality"] == "microbiome").sum()),
        "shared_metabolite_exact_name_candidates": int((shared["modality"] == "metabolome").sum()),
        "prism_manifest_origins": build_manifest(specs[2], axes["PRISM:microbiome"][0])["origin"]
        .value_counts()
        .to_dict(),
        "forbidden_clinical_rows": int(registry["is_forbidden_clinical"].sum()),
    }
    (output_dir / "contract_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary
