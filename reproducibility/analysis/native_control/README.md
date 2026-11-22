# Native SIAMCAT specimen-AUROC positive control

The ten archived probabilities per HMP2 specimen came from ten source seeds. Each export already averages the predictions of five Franzosa source folds. Their identifiers and probabilities remain in the original ten split files in ../import_case/native. The scoring adapter uses the registry only for the sample key, participant identity and target label. The registry direction name PRISM_to_HMP2 denotes Franzosa to HMP2. Across the ten files there are 3,850 predictions for 385 specimens from 105 people. The scoring labels are 281 IBD and 104 control specimens.

Run the following from the unpacked bundle root with R 4.5 and SIAMCAT 2.14.0 available on the R library path. The second argument is a new output directory.

```sh
Rscript --vanilla analysis/native_control/score_native.R . rerun/native_control
python analysis/native_control/verify_native.py . rerun/native_control
```

The script constructs a valid instance of the exported SIAMCAT S4 class using create.label() and the two scoring slots, label and pred_matrix. The class constructor and validObject() check the object. SIAMCAT::evaluate.predictions() scores it. All ten prediction columns match named IBD-positive labels by sample ID. mean(eval_data$auroc.all) is the average specimen AUROC over ten runs (A). eval_data$auroc is the AUROC after averaging the probabilities by specimen (B). This comparison uses complete-target AUROC. SIAMCAT's precision–recall curve area and the protocol's non-interpolated average precision have separate definitions.

results/native_run_auroc.csv records the ten run values. results/native_auroc.csv compares A and B with the full-precision registered values to an absolute tolerance of 1e-12. results/native_inputs.csv records the SHA-256 of every split, the registry and the protocol score reference. The Python check independently joins the original records and recomputes the ten run AUROCs and two corresponding summaries using scikit-learn. To check the supplied readouts without running R, use the second command with analysis/native_control/results as the output location argument.

## R environment

The actual scoring run used R 4.5.3 for Windows x86-64, Bioconductor release 3.22, SIAMCAT 2.14.0, and pROC 1.19.1. Package metadata records SIAMCAT commit cf500c96fcb1ae58030a027a468028ac385e0690 on RELEASE_3_22. The SIAMCAT Windows binary SHA-256 is 1b6970e1d00634f60074a0b9f94913055674b5dcd1ac6609f9d3350c53b10a97, and the signed Windows R 4.5.3 installer SHA-256 is 768ae31bb0b6056def5b1a9789a7dc49306bd037d69b0a99cdd90183aa0c1a31. The 95 downloaded Windows dependency archives, versions, binary source URLs and checksums are listed in native_windows_packages.csv. Recommended packages included with R are recorded in the actual scoring session in results/native_session.txt. The original model export environment has its own ../import_case/model_session.txt.

For Bioconductor-compatible R 4.5 systems, use a separate user library and run:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(version = "3.22", ask = FALSE)
BiocManager::install(c("SIAMCAT", "digest"), ask = FALSE, update = FALSE)
stopifnot(as.character(packageVersion("SIAMCAT")) == "2.14.0")
```

On Windows, the listed R 4.5 binaries identify the exact dependency builds used here. Package repositories may later move individual older dependency builds to archives. Installation from Bioconductor release 3.22 on another operating system can give different dependency build identifiers; verify the scoring CSVs with the specified tolerance and preserve the runtime-specific session record.

Package entry: https://bioconductor.posit.co/packages/3.22/bioc/html/SIAMCAT.html  
Scorer source: https://github.com/bioc/SIAMCAT/blob/cf500c96fcb1ae58030a027a468028ac385e0690/R/evaluate_predictions.r
