# Run from the unpacked evidence-bundle root; no feature matrix or trained object is needed.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript --vanilla score_native.R BUNDLE OUTPUT")
bundle <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
suppressPackageStartupMessages(library(SIAMCAT))
stopifnot(as.character(packageVersion("SIAMCAT")) == "2.14.0")
options(digits = 17)
# pROC computes auxiliary curve intervals internally; only the complete-panel AUROC is compared.
set.seed(20261121L)
regpath <- "evidence/full_reproducibility_inputs/prediction_registry.csv"
scorepath <- "analysis/import_case/results/scores.csv"
splits <- sprintf("analysis/import_case/native/split_%02d.csv", seq_len(10L))
read_input <- function(rel) read.csv(file.path(bundle, rel), check.names = FALSE)
registry <- read_input(regpath)
panel <- unique(registry[registry$direction == "PRISM_to_HMP2",
                         c("sample_id", "outcome_ibd", "bootstrap_cluster")])
panel <- panel[order(panel$sample_id, method = "radix"), ]
stopifnot(nrow(panel) == 385L, !anyNA(panel), !anyDuplicated(panel$sample_id),
          length(unique(panel$bootstrap_cluster)) == 105L,
          all(panel$outcome_ibd %in% 0:1), sum(panel$outcome_ibd) == 281L,
          all(tapply(panel$outcome_ibd, panel$bootstrap_cluster, function(x) length(unique(x))) == 1L))
runs <- lapply(splits, function(rel) {
    records <- read_input(rel)
    stopifnot(identical(names(records), c("sample_id", "probability")),
              nrow(records) == nrow(panel), !anyDuplicated(records$sample_id),
              setequal(records$sample_id, panel$sample_id),
              all(is.finite(records$probability)), all(records$probability >= 0 & records$probability <= 1))
    records$probability[match(panel$sample_id, records$sample_id)]
})
pred <- do.call(cbind, runs)
dimnames(pred) <- list(panel$sample_id, sprintf("split_%02d", seq_len(ncol(pred))))
outcome <- setNames(ifelse(panel$outcome_ibd == 1L, "IBD", "Control"), panel$sample_id)
lab <- SIAMCAT::create.label(label = outcome, case = "IBD", control = "Control", verbose = 0)
stopifnot(identical(names(lab$label), rownames(pred)), all((lab$label > 0) == (outcome == "IBD")))
# Instantiate the exported S4 class with just the scoring slots. The package validates the
# object and its public evaluator consumes these two slots directly.
sc <- methods::new("siamcat", label = lab, pred_matrix = pred)
stopifnot(methods::validObject(sc))
evaluation <- SIAMCAT::eval_data(SIAMCAT::evaluate.predictions(sc, verbose = 0))
native <- c(mean(evaluation$auroc.all), as.numeric(evaluation$auroc))
stopifnot(length(evaluation$auroc.all) == ncol(pred), all(is.finite(native)))
refs <- read_input(scorepath)
refs <- refs[refs$metric == "roc_auc" & refs$object %in% c("A", "B"), ]
stopifnot(nrow(refs) == 2L, !anyDuplicated(refs$object))
registered <- refs$complete_target[match(c("A", "B"), refs$object)]
difference <- abs(registered - native)
stopifnot(all(difference <= 1e-12))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) write.csv(x, file.path(output, name), row.names = FALSE, quote = TRUE)
precise <- function(x) sprintf("%.17g", x)
save_table(data.frame(object = c("A", "B"), statistic = c("mean(auroc.all)", "auroc"),
                      n_specimens = nrow(pred), n_people = length(unique(panel$bootstrap_cluster)),
                      n_runs = ncol(pred), positive_class = "IBD",
                      native_auroc = precise(native), protocol_auroc = precise(registered),
                      absolute_difference = precise(difference), tolerance = "1e-12"),
           "native_auroc.csv")
save_table(data.frame(run = colnames(pred), auroc = precise(as.numeric(evaluation$auroc.all))),
           "native_run_auroc.csv")
sources <- c(splits, regpath, scorepath)
save_table(data.frame(path = sources, sha256 = vapply(sources, function(x)
    digest::digest(readBin(file.path(bundle, x), "raw", n = file.size(file.path(bundle, x))), algo = "sha256", serialize = FALSE), character(1), USE.NAMES = FALSE)),
    "native_inputs.csv")
writeLines(capture.output(sessionInfo()), file.path(output, "native_session.txt"))
cat("SIAMCAT", as.character(packageVersion("SIAMCAT")), "AUROC A =", precise(native[1]),
    "B =", precise(native[2]), "max difference =", precise(max(difference)), "\n")
