# Minimal package inventory

## Included

- The complete reusable `src/ibd_microbiome_framework` package and its unit,
  contract, split, synthetic-control and Nielsen preparation tests.
- Primary standard model benchmark and reciprocal Nielsen runners.
- Source-only contract sensitivity, learning-curve, repeated-subject leakage,
  participant-unit, calibration and paired-bootstrap audits.
- Nielsen public-resource export and locked-target preparation utilities.
- SIAMCAT input export and native protocol-reference runner, clearly separated
  from the primary CLR benchmark and requiring its own external environment.
- P0-R closure, participant-unit, calibration, paired-bootstrap and subgroup
  feasibility audits used to verify the manuscript's declared estimands.
- Input schemas, sanitized configuration templates and reproducibility notes.

## Intentionally excluded

- Participant-level abundance matrices, manifests, predictions and frozen result
  CSV files.
- Manuscript DOCX/PDF files, figures, tables and author/funding metadata.
- Local workstation paths, temporary environments and generated caches.
- SIAMCAT executable/R environment. Its runner is included only as an optional
  protocol-level reference and requires a separately declared native environment.
- Exploratory domain-interaction, subtype, PRISM-origin and historical release
  runners that are not required for the frozen primary benchmark.

The formal six-direction interaction script is also excluded deliberately. The
manuscript reports descriptive observed cohort-edge/model-conditional ordering;
the excluded script could be misread as a universal interaction test.

The excluded runners remain in the local scaffold for provenance but are not
part of the smallest public code path. They should only be added in a separate
release if their estimands, dependencies and publication role are documented.
