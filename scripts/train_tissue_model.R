#!/usr/bin/env Rscript
# Train the cross-context sharing model of the Variant interpretation tab from
# TOPMed cis-QTL summary statistics (tensorQTL layout). Reads only; writes one
# JSON model file and its working tables.
#
#   Rscript scripts/train_tissue_model.R --type eqtl|sqtl --topmed DIR --out MODEL.json --work DIR
#     [--tabix PATH] [--strong 6000] [--random 5000] [--seed 20261112] [--cores 6] [--reuse]
# --reuse refits from the training table already in --work instead of reading the QTL files again.
#
# What it does (see docs/INTERPRETATION.md, "Cross-context model"):
#   1. Strong pairs: the lead variant of phenotypes with q <= 0.05, sampled
#      evenly from every context's permutation file.
#   2. Random pairs: for phenotypes drawn regardless of significance, a variant
#      20-200 kb from the lead variant, i.e. chosen without looking at its p.
#   3. Every pair's slope and standard error in all contexts, by tabix.
#   4. Null correlation between contexts from the random pairs; learned sharing
#      patterns from the strong pairs; mixture weights by EM on both. Four
#      fifths of the phenotypes are used. Fitted on the slope and on the z scale.
#   5. Check on the held-out fifth: each context's estimate predicted from the
#      other contexts, against two simple baselines, and interval coverage. The
#      scale that predicts held-out strong pairs better is kept.
suppressPackageStartupMessages({ library(parallel); library(jsonlite) })
argv <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  i <- match(paste0("--", name), argv)
  if (is.na(i)) { if (is.null(default)) stop("missing --", name, call. = FALSE); return(default) }
  argv[[i + 1L]]
}
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
backend <- new.env(); sys.source(file.path(dirname(dirname(normalizePath(self))), "src", "interpretation_backend.R"), envir = backend)
type <- match.arg(opt("type"), c("eqtl", "sqtl"))
topmed <- normalizePath(opt("topmed")); out_file <- opt("out"); work <- opt("work")
tabix <- opt("tabix", backend$interp_local_sources()$tabix)
n_strong <- as.integer(opt("strong", "6000")); n_random <- as.integer(opt("random", "5000"))
seed <- as.integer(opt("seed", "20261112")); cores <- as.integer(opt("cores", "6"))
dir.create(work, recursive = TRUE, showWarnings = FALSE); dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(tabix), dir.exists(topmed))
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

contexts <- backend$interp_topmed_contexts()            # file stem -> label
sub_dir <- if (type == "eqtl") "cis-eqtl" else "cis-sqtl2"
qtl_label <- if (type == "eqtl") "eQTL" else "sQTL"
nominal <- vapply(names(contexts), function(stem) file.path(topmed, sub_dir, "nominal", paste0(stem, ".maf001.tsv.gz")), character(1))
hits_file <- vapply(names(contexts), function(stem) {
  f <- Sys.glob(file.path(topmed, sub_dir, "top-hits", paste0(stem, "*.maf001.permutations.txt.gz")))
  if (length(f) != 1L) stop("expected one maf001 permutation file for ", stem); f
}, character(1))
stopifnot(all(file.exists(nominal)), all(file.exists(paste0(nominal, ".tbi"))))

table_file <- file.path(work, paste0("training_", type, ".rds"))
reuse <- "--reuse" %in% argv && file.exists(table_file)
if (reuse) {
  saved <- readRDS(table_file); pairs <- saved$pairs; Z <- saved$Z; B <- saved$B; SE <- saved$SE
  say("reusing", nrow(pairs), "pairs from", table_file)
} else {
# ---- 1-2. choose pairs -------------------------------------------------------
set.seed(seed)
hits <- do.call(rbind, lapply(names(contexts), function(stem) {
  x <- utils::read.delim(gzfile(hits_file[[stem]]), stringsAsFactors = FALSE)
  data.frame(context = stem, key = backend$interp_phenotype_key(rep(qtl_label, nrow(x)), x$phenotype_id), variant = x$variant_id, qval = x$qval, stringsAsFactors = FALSE)
}))
hits <- hits[grepl("^chr[0-9XY]+_[0-9]+_[ACGT]+_[ACGT]+$", hits$variant), ]
hits$chrom <- sub("_.*$", "", hits$variant); hits$pos <- as.integer(sub("^[^_]+_([0-9]+)_.*$", "\\1", hits$variant))
significant <- hits[is.finite(hits$qval) & hits$qval <= 0.05, ]
strong <- do.call(rbind, lapply(split(significant, significant$context), function(x) x[sample(nrow(x), min(nrow(x), ceiling(n_strong / length(contexts)))), ]))
strong <- strong[!duplicated(paste(strong$key, strong$variant)), c("key", "variant", "chrom", "pos", "context")]
say(qtl_label, "strong pairs:", nrow(strong), "from", nrow(significant), "significant lead pairs")

seeds <- hits[sample(nrow(hits), min(nrow(hits), ceiling(n_random * 1.6))), ]
seeds$start <- pmax(1L, seeds$pos + sample(c(-1L, 1L), nrow(seeds), TRUE) * sample(20000:200000, nrow(seeds), TRUE))
run_tabix <- function(file, bed, keep_awk, out) {
  # Rows of the regions in `bed`, reduced by an awk program to the wanted columns.
  status <- system2("bash", c("-c", shQuote(sprintf("set -o pipefail; %s -R %s %s 2>/dev/null | LC_ALL=C awk -F'\\t' %s > %s",
    shQuote(tabix), shQuote(bed), shQuote(file), shQuote(keep_awk), shQuote(out)))))
  if (status != 0L) stop("tabix failed on ", file)
  out
}
key_awk <- if (type == "eqtl") "k = $4" else "n = split($4, a, \":\"); k = a[1] \":\" a[2] \":\" a[3] \":\" a[n]"
random_parts <- mclapply(names(contexts), function(stem) {
  x <- seeds[seeds$context == stem, ]; if (!nrow(x)) return(NULL)
  bed <- file.path(work, paste0("random_", type, "_", stem, ".bed")); wanted <- file.path(work, paste0("random_", type, "_", stem, ".wanted"))
  x <- x[order(x$chrom, x$start), ]
  utils::write.table(data.frame(x$chrom, x$start - 1L, x$start + 500L), bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(unique(x$key), wanted)
  out <- run_tabix(nominal[[stem]], bed, sprintf("BEGIN { while ((getline line < \"%s\") > 0) want[line] = 1 } { %s; if ((k in want) && !(k in seen)) { seen[k] = 1; print k \"\\t\" $5 } }",
    wanted, key_awk), file.path(work, paste0("random_", type, "_", stem, ".pairs")))
  if (!file.info(out)$size) return(NULL)
  y <- utils::read.delim(out, header = FALSE, col.names = c("key", "variant"), stringsAsFactors = FALSE)
  y$context <- stem; y
}, mc.cores = cores)
failed <- vapply(random_parts, inherits, logical(1), "try-error")
if (any(failed)) stop("random-pair lookup failed: ", paste(unlist(random_parts[failed]), collapse = "; "))
random <- do.call(rbind, random_parts)
random <- random[grepl("^chr[0-9XY]+_[0-9]+_[ACGT]+_[ACGT]+$", random$variant) & !duplicated(paste(random$key, random$variant)), ]
random <- random[!paste(random$key, random$variant) %in% paste(strong$key, strong$variant), ]
if (nrow(random) > n_random) random <- random[sample(nrow(random), n_random), ]
random$chrom <- sub("_.*$", "", random$variant); random$pos <- as.integer(sub("^[^_]+_([0-9]+)_.*$", "\\1", random$variant))
say(qtl_label, "random pairs:", nrow(random))

# ---- 3. every pair in every context -------------------------------------------
pairs <- rbind(cbind(strong[, c("key", "variant", "chrom", "pos")], set = "strong"), cbind(random[, c("key", "variant", "chrom", "pos")], set = "random"))
pairs$id <- paste(pairs$variant, pairs$key, sep = "|")
bed <- file.path(work, paste0("pairs_", type, ".bed")); wanted <- file.path(work, paste0("pairs_", type, ".wanted"))
pos <- unique(pairs[, c("chrom", "pos")]); pos <- pos[order(pos$chrom, pos$pos), ]
utils::write.table(data.frame(pos$chrom, pos$pos - 1L, pos$pos), bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(pairs$id, wanted)
say("reading", nrow(pos), "positions in", length(contexts), "contexts")
stats <- mclapply(names(contexts), function(stem) {
  out <- run_tabix(nominal[[stem]], bed, sprintf("BEGIN { while ((getline line < \"%s\") > 0) want[line] = 1 } { %s; id = $5 \"|\" k; if (id in want) print id \"\\t\" $11 \"\\t\" $12 }",
    wanted, key_awk), file.path(work, paste0("stats_", type, "_", stem, ".tsv")))
  y <- if (file.info(out)$size) utils::read.delim(out, header = FALSE, col.names = c("id", "slope", "se"), stringsAsFactors = FALSE) else
    data.frame(id = character(), slope = numeric(), se = numeric())
  # An intron listed under two cluster numbers in one run: keep its first row.
  y[!duplicated(y$id) & is.finite(y$slope) & is.finite(y$se) & y$se > 0, ]
}, mc.cores = cores)
failed <- vapply(stats, inherits, logical(1), "try-error")
if (any(failed)) stop("pair lookup failed: ", paste(unlist(stats[failed]), collapse = "; "))
names(stats) <- names(contexts)
B <- SE <- matrix(NA_real_, nrow(pairs), length(contexts), dimnames = list(NULL, unname(contexts)))
for (i in seq_along(contexts)) { k <- match(pairs$id, stats[[i]]$id); B[, i] <- stats[[i]]$slope[k]; SE[, i] <- stats[[i]]$se[k] }
Z <- B / SE
tested <- rowSums(is.finite(Z))
say("pairs found in >= 2 contexts:", sum(tested >= 2L), "of", nrow(pairs))
keep <- tested >= 2L; pairs <- pairs[keep, ]; Z <- Z[keep, , drop = FALSE]; B <- B[keep, , drop = FALSE]; SE <- SE[keep, , drop = FALSE]
saveRDS(list(pairs = pairs, Z = Z, B = B, SE = SE), table_file)
}

# ---- 4. fit --------------------------------------------------------------------
# Hold out by phenotype, so a held-out pair's gene or intron was never fitted on.
phenotypes <- sort(unique(pairs$key)); set.seed(seed + 1L)
fold <- pairs$key %in% sample(phenotypes, round(length(phenotypes) / 5))
is_strong <- pairs$set == "strong"
ctx <- unname(contexts)
V <- backend$interp_mash_null_correlation(Z[!is_strong & !fold, , drop = FALSE]); dimnames(V) <- list(ctx, ctx)
chunks <- function(rows) split(rows, cut(seq_along(rows), min(cores, length(rows)), labels = FALSE))
fit <- function(scale) {
  X <- if (scale == "slope") B else Z
  U <- c(backend$interp_mash_canonical(ctx), backend$interp_mash_learned(X[is_strong & !fold, , drop = FALSE], V, ctx,
    SE = if (scale == "slope") SE[is_strong & !fold, , drop = FALSE] else NULL))
  model <- backend$interp_mash_model(ctx, V, U, scale = scale)
  # Weights are fitted on the random and the strong pairs together: the model is applied to
  # variants that are significant for the phenotype somewhere, which lie between the two.
  rows <- which(!fold)
  say(scale, "scale: EM on", length(rows), "pairs,", length(model$pi), "components")
  loglik <- do.call(rbind, mclapply(chunks(rows), function(r) backend$interp_mash_loglik(X[r, , drop = FALSE], model,
    if (scale == "slope") SE[r, , drop = FALSE] else NULL), mc.cores = cores))
  model$pi <- backend$interp_mash_fit_pi(loglik)
  model
}

# ---- 5. check on held-out pairs --------------------------------------------------
# Each tested context in turn is predicted from the others and compared in z units.
evaluate <- function(model, rows, label) {
  slope <- identical(model$scale, "slope")
  res <- do.call(rbind, mclapply(rows, function(i) {
    z <- Z[i, ]; o <- which(is.finite(z)); if (length(o) < 2L) return(NULL)
    do.call(rbind, lapply(o, function(leave) {
      p <- if (slope) backend$interp_mash_predict_left_out(B[i, ], model, leave, SE[i, ]) / SE[i, leave] else backend$interp_mash_predict_left_out(z, model, leave)
      others <- setdiff(o, leave); w <- 1 / SE[i, others]^2
      data.frame(context = ctx[[leave]], observed = z[[leave]], model = p[["mean"]], sd = p[["sd"]], mean_other_z = mean(z[others]),
        pooled_slope = sum(w * B[i, others]) / sum(w) / SE[i, leave], stringsAsFactors = FALSE)
    }))
  }, mc.cores = cores))
  summary_of <- function(x, name) data.frame(scale = model$scale, set = label, context = name, n = nrow(x),
    r_model = stats::cor(x$observed, x$model), r_mean_other_z = stats::cor(x$observed, x$mean_other_z), r_pooled_slope = stats::cor(x$observed, x$pooled_slope),
    rmse_model = sqrt(mean((x$observed - x$model)^2)), rmse_mean_other_z = sqrt(mean((x$observed - x$mean_other_z)^2)),
    rmse_pooled_slope = sqrt(mean((x$observed - x$pooled_slope)^2)), rmse_zero = sqrt(mean(x$observed^2)),
    coverage_95 = mean(abs(x$observed - x$model) <= 1.96 * x$sd),
    sign_agreement = { big <- abs(x$observed) > 3; if (any(big)) mean(sign(x$model[big]) == sign(x$observed[big])) else NA_real_ },
    stringsAsFactors = FALSE)
  rbind(summary_of(res, "all"), do.call(rbind, lapply(split(res, res$context), function(x) summary_of(x, x$context[[1L]]))))
}
models <- list(); evaluations <- list()
for (scale in c("slope", "z")) {
  models[[scale]] <- fit(scale)
  say(scale, "scale: checking on", sum(is_strong & fold), "held-out strong and", sum(!is_strong & fold), "held-out random pairs")
  evaluations[[scale]] <- rbind(evaluate(models[[scale]], which(is_strong & fold), "strong, held out"),
    evaluate(models[[scale]], which(!is_strong & fold), "random, held out"))
}
evaluation <- do.call(rbind, evaluations); rownames(evaluation) <- NULL
print(format(evaluation[evaluation$context == "all", ], digits = 3), row.names = FALSE)
# Keep the scale with the smaller error on held-out strong pairs.
error <- vapply(evaluations, function(e) e$rmse_model[e$set == "strong, held out" & e$context == "all"], numeric(1))
chosen <- names(error)[[which.min(error)]]
model <- models[[chosen]]
say("kept scale:", chosen, "(held-out strong RMSE:", paste(names(error), round(error, 3), collapse = ", "), ")")
comps <- backend$interp_mash_components(model)
weight_by_class <- tapply(model$pi, backend$interp_mash_class(comps$pattern), sum)

model$meta <- list(qtl_type = qtl_label, trained_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), seed = seed, source = topmed,
  files = basename(nominal), strong_pairs_fitted = sum(is_strong & !fold), random_pairs_fitted = sum(!is_strong & !fold),
  strong_pairs_held_out = sum(is_strong & fold), random_pairs_held_out = sum(!is_strong & fold),
  pairs_tested_per_context = as.list(stats::setNames(colSums(is.finite(Z)), ctx)),
  prior_weight_by_class = as.list(round(weight_by_class, 5)), evaluation = evaluation[evaluation$scale == chosen, ],
  evaluation_other_scale = evaluation[evaluation$scale != chosen & evaluation$context == "all", ],
  method = paste0("Mixture of multivariate normals on the ", chosen, " scale (x_hat ~ N(x, S V S); x ~ sum pi N(0, w^2 U)); fixed and learned sharing patterns; ",
    "weights by EM on random and strong pairs with null penalty 10; the other scale was fitted too and predicted held-out strong pairs less well."))
backend$interp_tissue_model_write(model, out_file)
say("model written:", out_file)
print(format(model$meta$evaluation, digits = 3), row.names = FALSE)
print(round(sort(weight_by_class, decreasing = TRUE), 4))
print(round(V, 3))
