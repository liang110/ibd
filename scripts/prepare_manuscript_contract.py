"""Build the manuscript's matched HMP2 source contract from permitted local files."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from ibd_microbiome_framework.data.contracts import MatrixSpec, sha256_file, write_contracts
from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list


# Fingerprints of the matched sample/person/endpoint mapping and ordered source features.
MATCHED_SHA = "53bb1cc1756ac1087fc7487c1449e1459eee018f0782ba6b88bfc4d6a0077415"
FEATURE_SHA = "f1a8cda85f11fb9094743b49924b6ba93515c1f0cc09e947f76d60711b0f23b4"
KEY = ["sample_id", "participant_id", "outcome_ibd"]
TARGET_DIRECTIONS = ("PRISM_to_HMP2", "NielsenHB_2014_to_HMP2")


def check(condition, message):
    if not condition:
        raise ValueError(message)


def specimen_columns(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        samples = next(csv.reader(handle))[1:]
    check(all(samples) and len(set(samples)) == len(samples), f"Invalid specimen IDs: {path.name}")
    return set(samples)


def compare_registered_targets(registry, retained):
    columns = ["direction", "sample_id", "bootstrap_cluster", "outcome_ibd"]
    rows = pd.read_csv(registry, usecols=columns, dtype=str).rename(
        columns={"bootstrap_cluster": "participant_id"}
    )
    reference = retained[KEY].astype(str).sort_values("sample_id").reset_index(drop=True)
    for direction in TARGET_DIRECTIONS:
        actual = rows.loc[rows.direction.eq(direction), KEY].drop_duplicates()
        actual = actual.sort_values("sample_id").reset_index(drop=True)
        check(actual.equals(reference), f"Specimen, participant or outcome mismatch in {direction}")


def build(root, provider_metadata, registry, work, destination):
    abundance = root / "manual_abundance.csv"
    metabolites = root / "metabolites.csv"
    microbial_ids = specimen_columns(abundance)
    metabolite_ids = specimen_columns(metabolites)
    provider_cols = ["External ID", "Participant ID", "diagnosis", "site_sub_coll", "data_type"]
    raw = pd.read_csv(provider_metadata, usecols=provider_cols, dtype=str)
    provider = raw.loc[raw.data_type.str.lower().eq("metagenomics") & raw["External ID"].notna()]
    provider = provider.set_index("External ID", verify_integrity=True)
    check(set(provider.index) == microbial_ids, "Provider and abundance profile IDs disagree")
    check(provider.notna().all().all(), "Missing provider participant, diagnosis or collection")
    check(set(provider.diagnosis) == {"CD", "UC", "nonIBD"}, "Unexpected provider diagnosis")
    by_collection = provider.groupby("site_sub_coll")[["Participant ID", "diagnosis"]].nunique()
    check(by_collection.eq(1).all().all(), "Collection has conflicting participants or diagnoses")
    check(provider.groupby("Participant ID").diagnosis.nunique().eq(1).all(), "Diagnosis conflict")
    check(
        (len(provider), len(by_collection), provider["Participant ID"].nunique(), len(metabolite_ids))
        == (1638, 1595, 130, 546),
        "Input inventory differs from the manuscript release",
    )

    clean = pd.read_csv(root / "clean_metadata.csv", dtype=str)
    check(clean.sample_id.notna().all() and clean.sample_id.is_unique, "Invalid clean sample IDs")
    matched = clean.loc[clean.sample_id.isin(microbial_ids & metabolite_ids)].copy()
    selected = provider.loc[matched.sample_id]
    check((len(selected), selected["Participant ID"].nunique()) == (388, 105), "Unexpected panel")
    check(selected.site_sub_coll.is_unique, "Multiple matched profiles map to the same collection")
    for field in ("Participant ID", "diagnosis"):
        check((selected[field].to_numpy() == matched[field].to_numpy()).all(), f"Conflicting {field}")

    filtered = work / "hmp2_matched_metadata.csv"
    matched.to_csv(filtered, index=False)
    specs = [
        MatrixSpec("HMP2", "microbiome", abundance, "Feature", filtered),
        MatrixSpec("PRISM", "microbiome", root / "prism_abundance_clean.csv",
                   "Species", root / "prism_metadata_clean.csv"),
    ]
    write_contracts(specs, work, root)
    inventory_file = work / "source_inventory.csv"
    inventory = pd.read_csv(inventory_file)
    inventory.loc[inventory.cohort.eq("HMP2"), "metadata_path"] = str(destination / filtered.name)
    inventory.to_csv(inventory_file, index=False)

    manifest = pd.read_csv(work / "cohort_manifest_hmp2.csv")
    encoded = manifest[KEY].sort_values("sample_id").to_csv(index=False, lineterminator="\n")
    matched_sha = hashlib.sha256(encoded.encode("utf-8")).hexdigest()
    check(matched_sha == MATCHED_SHA, "Matched mapping differs from the original source manifest")
    features = load_shared_feature_list(work)
    feature_sha = hashlib.sha256("\n".join(features).encode("utf-8")).hexdigest()
    check(len(features) == 126 and feature_sha == FEATURE_SHA, "Ordered feature contract differs")
    projected = load_canonical_microbiome(abundance, manifest, features)
    retained = manifest.loc[manifest.sample_id.isin(projected.index)].copy()
    excluded = manifest.loc[~manifest.sample_id.isin(projected.index)].copy()
    people = retained.drop_duplicates("participant_id")
    check((len(retained), len(people), retained.outcome_ibd.sum(), people.outcome_ibd.sum())
          == (385, 105, 281, 79), "Unexpected retained population")
    check(excluded.diagnosis.value_counts().to_dict() == {"CD": 2, "nonIBD": 1},
          "Unexpected zero-mass exclusions")
    compare_registered_targets(registry, retained)

    eligibility = provider[["site_sub_coll"]].copy()
    eligibility["status"] = "metabolomics_id_absent"
    eligibility.loc[eligibility.index.isin(metabolite_ids), "status"] = "clean_metadata_absent"
    eligibility.loc[excluded.sample_id, "status"] = "zero_shared_mass"
    eligibility.loc[retained.sample_id, "status"] = "retained"
    eligibility.to_csv(work / "hmp2_profile_eligibility.csv", index_label="sample_id")
    retained.to_csv(work / "hmp2_retained_manifest.csv", index=False)
    excluded["reason"] = "zero_shared_mass"
    excluded.to_csv(work / "hmp2_exclusions.csv", index=False)
    report = {
        "metagenomic_profiles": len(provider), "collections": len(by_collection),
        "inventory_participants": int(provider["Participant ID"].nunique()),
        "metabolomics_ids": len(metabolite_ids),
        "matched_profiles_and_unique_collections": len(matched),
        "matched_participants": int(selected["Participant ID"].nunique()),
        "retained_profiles_and_unique_collections": len(retained), "retained_participants": len(people),
        "retained_case_control_specimens": [int(retained.outcome_ibd.eq(i).sum()) for i in (1, 0)],
        "retained_case_control_participants": [int(people.outcome_ibd.eq(i).sum()) for i in (1, 0)],
        "excluded_zero_shared_mass": excluded.diagnosis.value_counts().to_dict(),
        "matched_mapping_sha256": matched_sha, "ordered_126_features_sha256": feature_sha,
        "prediction_registry_sha256": sha256_file(registry),
        "both_hmp2_target_mappings_equal": True,
        "eligibility_source_sha256": {
            path.name: sha256_file(path)
            for path in (abundance, metabolites, provider_metadata, root / "clean_metadata.csv")
        },
    }
    (work / "hmp2_selection_check.json").write_text(json.dumps(report, indent=2) + "\n",
                                                     encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path, help="Flat permitted source directory")
    parser.add_argument("--provider-metadata", required=True, type=Path, help="Original IBDMDB metadata CSV")
    parser.add_argument("--prediction-registry", required=True, type=Path, help="Reader ZIP registry CSV")
    parser.add_argument("--output-dir", required=True, type=Path, help="New private source contract directory")
    args = parser.parse_args()
    destination = args.output_dir.resolve()
    if destination.exists():
        parser.error("--output-dir already exists; choose a new directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=destination.parent) as temporary:
        staging = Path(temporary) / "contract"
        staging.mkdir()
        build(args.source_root.resolve(), args.provider_metadata.resolve(),
              args.prediction_registry.resolve(), staging, destination)
        staging.rename(destination)
    print(f"Verified HMP2 matched and retained populations. Contract: {destination}")


if __name__ == "__main__":
    main()
