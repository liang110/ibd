# Input schema

The repository contains no participant-level data. The primary runners expect
the following user-supplied files:

| File | Required fields or orientation |
|---|---|
| Microbiome abundance matrix | taxa/features in rows; sample IDs in columns; first column is the raw feature ID |
| Cohort manifest | `sample_id`, `group`; `participant_id` is required for grouped splitting |
| Shared-feature contract | `modality`, `canonical_feature_id`; microbiome IDs must be unique |
| Nielsen projection | samples in rows and the frozen canonical features in columns |
| Nielsen manifest | `sample_id`, `outcome_ibd`, `participant_id`; optional `origin` and `site` are retained for audits |

Abundance values must be finite and non-negative. Clinical variables,
participant identifiers and outcome-derived columns are rejected as features.
The canonicalizer is deliberately conservative: it normalizes taxonomic
strings and does not use an external synonym database.

All externally supplied files should be accompanied by SHA-256 values in the
run manifest. The target cohort is projected to the already-frozen contract;
the target does not define prevalence, pseudocount, scaling, feature ranking or
model parameters.
