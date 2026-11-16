# Golden outputs captured from e6d2ed6 using the existing synthetic scenarios.
# The file holds every scenario in call order as run with bcftools, tabix and bgzip (75 on SCC);
# an environment that skips the tool-backed scenarios produces a subset, which is why order is
# compared only when the counts agree.
# Compare legacy columns, and verify that metadata never reconstructs text.
expressions <- parse("tests/test_interpretation.R")
eval(expressions[[1L]], .GlobalEnv)
original_summary <- interp_evidence_summary
golden <- jsonlite::fromJSON("tests/fixtures/evidence_summary_golden.json", simplifyVector = FALSE)
canonical <- function(x) as.character(jsonlite::toJSON(x, dataframe = "rows", auto_unbox = TRUE, na = "null"))
expected <- vapply(golden, canonical, character(1)); contract_seen <- character()
interp_evidence_summary <- function(...) {
  out <- original_summary(...)
  signature <- canonical(out)
  stopifnot(signature %in% expected, identical(names(out), c("question", "level", "finding", "evidence", "sources")))
  contract_seen <<- c(contract_seen, signature)
  args <- list(...); names(args) <- names(match.call(definition = original_summary, call = sys.call()))[-1L]
  for (step in c("public", "spliceai", "ld")) {
    changed_args <- args; changed_args[step] <- list(NULL)
    changed <- do.call(original_summary, changed_args)
    outside <- !vapply(out$question, function(q) step %in% interp_summary_scope()[[q]], logical(1))
    stopifnot(identical(out[outside, ], changed[outside, ]))
  }
  detail <- original_summary(..., structured = TRUE)
  stopifnot(identical(out, detail$summary))
  for (i in seq_len(nrow(out))) {
    items <- detail$items[[i]]
    split <- strsplit(out$evidence[[i]], " | ", fixed = TRUE)[[1L]]
    stopifnot(identical(unname(items$text), split[nzchar(split)]), !any(items$source == "unattributed"),
      all(items$source %in% interp_summary_sources()$key), all(items$role %in% c("evidence", "notice", "context")))
  }
  out
}
for (expression in expressions[-1L]) eval(expression, .GlobalEnv)
if (length(contract_seen) == length(expected)) stopifnot(identical(contract_seen, expected))
cat("SUMMARY_CONTRACT_OK:", length(contract_seen), "legacy outputs unchanged; structured lines match byte-for-byte\n")
