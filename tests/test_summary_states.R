source("src/interpretation_backend.R")
source("src/interpretation_module.R")
source("tests/fixtures/summary_fixture.R")
expected_scope <- list(`Protein sequence` = c("public", "spliceai"), `Protein abundance (pQTL)` = c("public", "spliceai", "ld"),
  Splicing = c("public", "spliceai"), `Gene expression (eQTL)` = "public", `Disease and traits` = c("public", "ld"), Pathway = "public", Tissue = "public")
stopifnot(identical(interp_summary_scope(), expected_scope))
scenarios <- list(before = interp_request_states(), public = interp_request_states(list(queried_at = "synthetic time")),
  after = interp_request_states(list(), list(), list()), failed = interp_request_states(simpleError("synthetic failure")))
for (name in names(scenarios)) {
  states <- scenarios[[name]]; f <- summary_stress_fixture(states)
  for (i in seq_len(7L)) {
    q <- f$summary$question[[i]]
    html <- as.character(interpretation_summary_cards(f$summary[i, ], items = f$items[i], states = states, identifiers = f$identifiers))
    for (step in interp_request_steps()$key) {
      present <- grepl(paste0('data-step="', step, '" data-state="', states[[step]], '"'), html, fixed = TRUE)
      stopifnot(identical(present, step %in% expected_scope[[q]]))
      if (present) stopifnot(grepl(interp_step_label(step, states[[step]]), html, fixed = TRUE))
    }
    for (line in f$items[[i]]$text) stopifnot(lengths(regmatches(html, gregexpr(htmltools::htmlEscape(line), html, fixed = TRUE))) == 1L)
    # Notices stay folded; context is a preview only when the row has no evidence.
    preview <- sub('.*<div class="regshiny-evidence-preview">', '', html)
    preview <- strsplit(preview, '</div>', fixed = TRUE)[[1L]][[1L]]
    stopifnot(lengths(regmatches(preview, gregexpr('<li ', preview, fixed = TRUE))) <= 2L,
      !grepl('data-role="notice"', preview, fixed = TRUE))
    if (any(f$items[[i]]$role == "evidence")) stopifnot(!grepl('data-role="context"', preview, fixed = TRUE)) else
      stopifnot(grepl('data-role="context"', preview, fixed = TRUE))
  }
  graph <- as.character(interpretation_summary_map(f$summary, f$items, states, f$identifiers))
  stopifnot(lengths(regmatches(graph, gregexpr('class="regshiny-map-question ', graph, fixed = TRUE))) == 7L)
  if (name == "before") stopifnot(grepl('regshiny-map-edge pending', graph, fixed = TRUE))
  if (name == "after") stopifnot(!grepl('regshiny-map-edge pending', graph, fixed = TRUE), grepl('data-step="public"', graph, fixed = TRUE))
  csv <- interp_summary_csv(f$summary, states)
  stopifnot(identical(csv[, names(f$summary)], f$summary))
  for (i in 1:7) stopifnot(csv$on_request_steps[[i]] == paste(paste0(expected_scope[[f$summary$question[[i]]]], ": ", states[expected_scope[[f$summary$question[[i]]]]]), collapse = "; "))
}
f <- summary_stress_fixture(scenarios$after)
for (key in interp_source_links()$key) {
  u <- interp_source_url(key, f$identifiers);stopifnot(!is.null(u))
  item <- data.frame(source = key, gene_id = "ENSG00000093010", gene = "COMT", clinvar_id = "17591")
  html <- as.character(interpretation_evidence_link("a verified source", item, f$identifiers))
  stopifnot(grepl(htmltools::htmlEscape(u, attribute = TRUE), html, fixed = TRUE), grepl('target="_blank"', html, fixed = TRUE), grepl('rel="noopener noreferrer"', html, fixed = TRUE))
}
for (key in c("topmed", "fhs", "cohort_ld", "gencode")) stopifnot(is.null(interp_source_url(key, f$identifiers)),
  !grepl('href=', as.character(interpretation_evidence_link("local", data.frame(source = key, gene_id = "", gene = ""), f$identifiers)), fixed = TRUE))
stopifnot(is.null(interp_source_url("open_targets", list(chrom = "22", pos = "19963748?sample=private", ref = "G", alt = "A"))),
  identical(interp_source_url("open_targets", c(f$identifiers, list(sample = "private", cohort = "private"))), interp_source_url("open_targets", f$identifiers)))
b <- summary_stress_fixture()$summary
html <- as.character(interpretation_summary_cards(f$summary, items = f$items, states = scenarios$after, baseline = b))
stopifnot(grepl("Moderate before the on-request steps", html, fixed = TRUE),
  !grepl("before the on-request steps", as.character(interpretation_summary_cards(f$summary, items = f$items, states = scenarios$after)), fixed = TRUE))
cat("SUMMARY_STATES_OK: four states, scopes, DOM preservation, links, Map, CSV, exact baseline display\n")
