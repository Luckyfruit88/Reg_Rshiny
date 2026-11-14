#!/usr/bin/env Rscript
# Interpret a list of variants from the local sources only, outside the app: one CSV row per variant with the
# level of each question, the predicted splice outcome and the mechanism lines of the Variant interpretation
# tab. Reads the GENCODE index, the reference FASTA and the local reference copies; no BAM, no genotype, no
# saved result and no network.
#
#   Rscript scripts/interpret_variants.R --variants FILE --annotation-db annotation.sqlite --reference genome.fa --out TABLE.csv
#     [--release 48] [--json-dir DIR]
#
# FILE holds one variant per line as chrom:pos:REF:ALT (GRCh38, 1-based, as in the app); empty lines and lines
# starting with # are skipped. --annotation-db is a GENCODE index prepared by the app
# (.reg_shiny/<user>/resources/gencode<release>-*/annotation.sqlite) and --reference its indexed FASTA; both
# default to REGSHINY_ANNOTATION_DB and REGSHINY_REFERENCE_FASTA. Rows are appended as they are computed, so an
# interrupted run keeps what it has; variants already in TABLE.csv are skipped. --json-dir also writes each
# variant's full tables as JSON. A variant takes about 10 to 30 seconds: run long lists as a cluster job.
argv <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  i <- match(paste0("--", name), argv)
  if (is.na(i) || i == length(argv)) { if (is.null(default)) stop("missing --", name, call. = FALSE); return(default) }
  argv[[i + 1L]]
}
if (!length(argv) || any(argv %in% c("-h", "--help"))) {
  cat("Rscript scripts/interpret_variants.R --variants FILE --annotation-db annotation.sqlite --reference genome.fa --out TABLE.csv [--release 48] [--json-dir DIR]\n")
  quit(status = if (length(argv)) 0L else 2L)
}
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
root <- dirname(dirname(normalizePath(self)))
variants_file <- opt("variants"); out_file <- opt("out")
db <- opt("annotation-db", Sys.getenv("REGSHINY_ANNOTATION_DB", "")); fasta <- opt("reference", Sys.getenv("REGSHINY_REFERENCE_FASTA", ""))
release <- opt("release", ""); json_dir <- opt("json-dir", "")
if (!nzchar(db) || !nzchar(fasta)) stop("give --annotation-db and --reference (or set REGSHINY_ANNOTATION_DB and REGSHINY_REFERENCE_FASTA)", call. = FALSE)
ids <- trimws(readLines(variants_file, warn = FALSE)); ids <- unique(ids[nzchar(ids) & !startsWith(ids, "#")])
bad <- ids[!grepl("^[^:[:space:]]+:[0-9]+:[ACGTNacgtn]+:[ACGTNacgtn]+$", ids)]
if (length(bad)) stop("not chrom:pos:REF:ALT: ", paste(utils::head(bad, 5L), collapse = ", "), call. = FALSE)
# The app's backends resolve their helper scripts from the repository root.
old <- setwd(root)
e <- new.env(); sys.source(file.path("src", "interpretation_backend.R"), envir = e)
a <- new.env(); sys.source(file.path("src", "annotation_backend.R"), envir = a)
s <- new.env(); sys.source(file.path("src", "splice_evidence_backend.R"), envir = s)
setwd(old)
annotation <- a$load_annotation_resources(db, if (nzchar(release)) release else NULL)
reference <- s$load_splice_reference(fasta)
fetch <- function(chrom, start1, end1) s$splice_fetch_reference(reference, chrom, start1, end1)$sequence
gene_names <- function(gene_ids) a$annotation_gene_names(annotation, gene_ids)
done <- if (file.exists(out_file)) tryCatch(utils::read.csv(out_file, stringsAsFactors = FALSE)$variant, error = function(err) character()) else character()
dir.create(dirname(normalizePath(out_file, mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
if (nzchar(json_dir)) dir.create(json_dir, recursive = TRUE, showWarnings = FALSE)
todo <- setdiff(ids, done)
cat(length(ids), "variant(s),", length(ids) - length(todo), "already in", out_file, "\n")
plain <- function(x) if (inherits(x, "error")) list(error = conditionMessage(x)) else x
for (id in todo) {
  f <- strsplit(id, ":", fixed = TRUE)[[1L]]; started <- Sys.time()
  row <- tryCatch({
    v <- e$interp_variant(f[[1L]], as.integer(f[[2L]]), f[[3L]], f[[4L]])
    report <- e$interp_local_report(v, function(chrom, start1, end1) a$query_annotation(annotation, chrom, start1, end1, "GRCh38"), fetch, gene_names = gene_names)
    if (nzchar(json_dir)) {
      full <- report; full$local <- lapply(full$local, plain); full$splice_protein <- plain(full$splice_protein); full$colocalisation <- plain(full$colocalisation)
      writeLines(jsonlite::toJSON(full, dataframe = "rows", auto_unbox = TRUE, na = "null", null = "null", digits = NA, pretty = TRUE),
        file.path(json_dir, paste0("interpretation_", gsub(":", "_", id, fixed = TRUE), ".json")))
    }
    e$interp_report_row(report)
  }, error = function(err) { r <- e$interp_report_row(list(variant = list(id = id), summary = data.frame(question = character(), level = character(), finding = character()),
      gene_summary = data.frame(), local = list(), linked_genes = data.frame(gene = character(), gene_id = character()), splice_protein = NULL, colocalisation = NULL,
      mechanism = character(), failed = paste("failed:", conditionMessage(err)))); r })
  utils::write.table(row, out_file, sep = ",", row.names = FALSE, col.names = !file.exists(out_file), append = file.exists(out_file), qmethod = "double", na = "")
  cat(format(Sys.time(), "%H:%M:%S"), id, round(as.numeric(Sys.time() - started, units = "secs")), "s", if (nzchar(row$sources_failed)) paste("|", row$sources_failed) else "", "\n")
}
cat("written:", out_file, "\n")
