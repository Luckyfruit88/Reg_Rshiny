# Entirely synthetic UI stress fixture: no cohort rows or private identifiers.
summary_stress_fixture <- function(states = interp_request_states()) {
  questions <- names(interp_summary_scope()); rows <- list()
  text <- function(prefix, n) substr(paste(prefix, paste(rep("synthetic evidence describes a reproducible public observation", 12), collapse = " ")), 1, n)
  for (i in seq_along(questions)) {
    q <- questions[[i]]; steps <- interp_summary_scope()[[q]]
    primary <- c("cadd", "pqtl", "topmed", "gtex_local", "clinvar_local", "reactome_local", "gtex_local")
    secondary <- c("gencode", "cohort_rna", "fhs", "topmed", "gwas_catalog", "msigdb", "topmed")
    lines <- list(interp_evidence_tag(text(paste0("Case ", i, ".1:"), 120), primary[[i]], gene = "COMT", clinvar_id = "17591"),
      interp_evidence_tag(text(paste0("Case ", i, ".2:"), 330), secondary[[i]]))
    for (step in steps) lines[[length(lines) + 1L]] <- interp_evidence_tag(
      text(paste0("Case ", i, ".", step, if (states[[step]] == "after") ": measured evidence." else ": operation not run."), 90),
      switch(step, public = "open_targets", spliceai = "spliceai", ld = "cohort_ld"), step,
      if (states[[step]] == "after") "evidence" else "notice")
    lines[[length(lines) + 1L]] <- interp_evidence_tag(text(paste0("Case ", i, ".context:"), 90), "local_pathways", role = "context")
    while (length(lines) < 8L) lines[[length(lines) + 1L]] <- interp_evidence_tag(text(paste0("Case ", i, ".", length(lines) + 1L, ":"), 90), "fhs")
    # Pathway deliberately has context only, matching the reported production case.
    if (q == "Pathway") lines <- lapply(lines, function(x) {
      meta <- attr(x, "evidence_items"); meta$role[meta$role == "evidence"] <- "context"
      attr(x, "evidence_items") <- meta; x
    })
    rows[[i]] <- interp_summary_row(q, if (i >= 6L) "Context" else if (any(states[steps] == "after")) "Strong" else "Moderate",
      text(paste0("Synthetic finding ", i, ":"), c(60, 120, 180, 240, 300, 360, 480)[[i]]),
      do.call(interp_evidence_combine, lines), c("GENCODE + reference (computed here)", "TOPMed"))
  }
  items <- lapply(rows, attr, "evidence_items"); rows <- lapply(rows, function(r) {attr(r, "evidence_items") <- NULL; r})
  list(summary = do.call(rbind, rows), items = items,
    identifiers = list(chrom = "22", pos = 19963748L, ref = "G", alt = "A", rs_id = "rs4680", gene_id = "ENSG00000093010", gene = "COMT", clinvar_id = "17591"))
}

# Two gene-specific Reactome targets and one unlinked MSigDB source. Tissue
# repeats MSigDB to check that node merging preserves all question connections.
summary_source_dedup_fixture <- function() {
  entries <- interp_evidence_combine(
    interp_evidence_tag("COMT Reactome context", "reactome_local", role = "context", gene = "COMT"),
    interp_evidence_tag("TXNRD2 Reactome context", "reactome_local", role = "context", gene = "TXNRD2"),
    interp_evidence_tag("COMT MSigDB context", "msigdb", role = "context", gene = "COMT"),
    interp_evidence_tag("TXNRD2 MSigDB context", "msigdb", role = "context", gene = "TXNRD2"))
  items <- lapply(names(interp_summary_scope()), function(q) {
    if (q == "Pathway") interp_evidence_items(entries) else if (q == "Tissue") interp_evidence_items(entries)[3:4, ] else interp_empty_evidence_items()
  })
  rows <- lapply(seq_along(items), function(i) data.frame(question = names(interp_summary_scope())[[i]], level = "Context",
    finding = "Synthetic source identity check.", evidence = paste(items[[i]]$text, collapse = " | "),
    sources = if (i == 6L) "Reactome; MSigDB" else if (i == 7L) "MSigDB" else "", stringsAsFactors = FALSE))
  list(summary = do.call(rbind, rows), items = items, identifiers = list())
}
