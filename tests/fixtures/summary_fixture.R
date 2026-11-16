# Entirely synthetic UI stress fixture: no cohort rows or private identifiers.
summary_stress_fixture <- function(states = interp_request_states()) {
  questions <- names(interp_summary_scope()); rows <- list()
  text <- function(prefix, n) substr(paste(prefix, paste(rep("synthetic evidence describes a reproducible public observation", 12), collapse = " ")), 1, n)
  for (i in seq_along(questions)) {
    q <- questions[[i]]; steps <- interp_summary_scope()[[q]]
    lines <- list(interp_evidence_tag(text(paste0("Case ", i, ".1:"), 120), "gencode"),
      interp_evidence_tag(text(paste0("Case ", i, ".2:"), 330), "topmed"))
    for (step in steps) lines[[length(lines) + 1L]] <- interp_evidence_tag(
      text(paste0("Case ", i, ".", step, if (states[[step]] == "after") ": measured evidence." else ": operation not run."), 90),
      switch(step, public = "open_targets", spliceai = "spliceai", ld = "cohort_ld"), step,
      if (states[[step]] == "after") "evidence" else "notice")
    lines[[length(lines) + 1L]] <- interp_evidence_tag(text(paste0("Case ", i, ".context:"), 90), "local_pathways", role = "context")
    while (length(lines) < 8L) lines[[length(lines) + 1L]] <- interp_evidence_tag(text(paste0("Case ", i, ".", length(lines) + 1L, ":"), 90), "fhs")
    rows[[i]] <- interp_summary_row(q, if (i >= 6L) "Context" else if (any(states[steps] == "after")) "Strong" else "Moderate",
      text(paste0("Synthetic finding ", i, ":"), c(60, 120, 180, 240, 300, 360, 480)[[i]]),
      do.call(interp_evidence_combine, lines), c("GENCODE + reference (computed here)", "TOPMed"))
  }
  items <- lapply(rows, attr, "evidence_items"); rows <- lapply(rows, function(r) {attr(r, "evidence_items") <- NULL; r})
  list(summary = do.call(rbind, rows), items = items,
    identifiers = list(chrom = "22", pos = 19963748L, ref = "G", alt = "A", rs_id = "rs4680", gene_id = "ENSG00000093010", gene = "COMT"))
}
