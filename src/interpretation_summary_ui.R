# Table and Map views of the same structured evidence; never perform a request.
interpretation_step_tags <- function(states, steps = interp_request_steps()$key) {
  shiny::tagList(lapply(steps, function(k) shiny::span(class = paste("regshiny-step", paste0("regshiny-step-", states[[k]])),
    `data-step` = k, `data-state` = states[[k]], interp_step_label(k, states[[k]]))))
}
interpretation_summary_identifiers <- function(variant = NULL, rs_id = "") {
  if (is.null(variant)) return(list())
  list(chrom = variant$chrom, pos = variant$pos1, ref = variant$ref, alt = variant$alt, rs_id = rs_id)
}
interpretation_source_info <- function(key) {
  registry <- interp_summary_sources(); i <- match(key, registry$key)
  if (is.na(i)) list(label = key, kind = "local") else as.list(registry[i, , drop = FALSE])
}
interpretation_item_identifiers <- function(item, identifiers) {
  for (key in c("gene_id", "gene", "clinvar_id")) {
    value <- item[[key]]
    if (length(value) == 1L && !is.na(value) && nzchar(value)) identifiers[[key]] <- value
  }
  identifiers
}
interpretation_source_note <- function(source) {
  registry <- interp_source_links(); i <- match(source, registry$key)
  if (!is.na(i) && nzchar(registry$note[[i]])) registry$note[[i]] else NULL
}
interpretation_evidence_link <- function(text, item, identifiers, ...) {
  # A not-run notice is never a link, even if its future source has a URL.
  if (isTRUE(item$role == "notice")) return(shiny::span(text, ...))
  ids <- interpretation_item_identifiers(item, identifiers)
  url <- interp_source_url(item$source, ids)
  if (is.null(url)) shiny::span(text, ...) else shiny::tags$a(text, href = url, target = "_blank", rel = "noopener noreferrer",
    title = interpretation_source_note(item$source), ...)
}
interpretation_summary_items <- function(summary, items) {
  if (!is.null(items)) return(items)
  # Compatibility for callers supplying a plain synthetic/legacy summary. No
  # attempt is made to infer provenance or roles from its English sentences.
  lapply(summary$evidence, function(x) interp_evidence_items(strsplit(x, " | ", fixed = TRUE)[[1L]]))
}
interpretation_source_badges <- function(items, identifiers, fallback = "") {
  items <- unique(items[items$role != "notice", c("source", "gene_id", "gene", "clinvar_id"), drop = FALSE])
  if (!nrow(items) || all(items$source == "unattributed")) return(if (nzchar(fallback)) shiny::span(paste("Sources:", fallback), shiny::tags$small(" · local")) else NULL)
  shiny::tagList(lapply(seq_len(nrow(items)), function(j) {
    item <- items[j, ]; info <- interpretation_source_info(item$source)
    shiny::div(class = "regshiny-source-badge", interpretation_evidence_link(info$label, item, identifiers),
      if (info$kind != "external") shiny::tags$small(paste0(" · ", info$kind)))
  }))
}
interpretation_summary_table <- function(summary, shown = interp_summary_visible_evidence(), items = NULL,
    states = interp_request_states(), baseline = NULL, identifiers = list()) {
  items <- interpretation_summary_items(summary, items)
  shiny::div(class = "regshiny-evidence", lapply(seq_len(nrow(summary)), function(i) {
    x <- items[[i]]; eligible <- which(x$role == "evidence")
    if (!length(eligible)) eligible <- which(x$role == "context")
    visible <- head(eligible, shown); hidden <- setdiff(seq_len(nrow(x)), visible)
    line <- function(j, preview = FALSE) shiny::tags$li(`data-role` = x$role[[j]], `data-source` = x$source[[j]], `data-step` = x$step[[j]],
      interpretation_evidence_link(x$text[[j]], x[j, ], identifiers, class = if (preview) "regshiny-evidence-line" else NULL))
    prior <- if (!is.null(baseline)) baseline$level[match(summary$question[[i]], baseline$question)] else NA_character_
    shiny::div(class = "regshiny-evidence-row", `data-question` = summary$question[[i]],
      shiny::div(class = "regshiny-evidence-head",
        shiny::strong(summary$question[[i]]),
        shiny::div(shiny::span(class = paste("regshiny-level", interpretation_level_class(summary$level[[i]])), summary$level[[i]]),
          if (length(prior) == 1L && !is.na(prior) && prior != summary$level[[i]]) shiny::tags$small(class = "regshiny-prior-level", paste(prior, "before the on-request steps"))),
        shiny::div(class = "regshiny-row-steps", interpretation_step_tags(states, interp_summary_scope()[[summary$question[[i]]]]))),
      shiny::div(class = "regshiny-evidence-body",
        shiny::p(class = "regshiny-finding", summary$finding[[i]]),
        shiny::div(class = "regshiny-evidence-preview", if (length(visible)) shiny::tags$ul(lapply(visible, line, preview = TRUE))),
        shiny::tags$details(class = "regshiny-evidence-more",
          shiny::tags$summary("Full details"), if (length(hidden)) shiny::tags$ul(lapply(hidden, line)))),
      shiny::div(class = "regshiny-evidence-source", interpretation_source_badges(x, identifiers, summary$sources[[i]])))
  }))
}

interpretation_summary_map <- function(summary, items = NULL, states = interp_request_states(), identifiers = list()) {
  items <- interpretation_summary_items(summary, items)
  nodes <- list(); edges <- list()
  add <- function(i, key, step, item, pending = FALSE) {
    id <- paste(step, key, item$gene_id, item$gene, item$clinvar_id, sep = ":")
    if (is.null(nodes[[id]])) nodes[[id]] <<- list(key = key, step = step, item = item, pending = pending)
    edges[[length(edges) + 1L]] <<- list(question = i, node = id, pending = pending)
  }
  for (i in seq_len(nrow(summary))) {
    x <- unique(items[[i]][items[[i]]$role != "notice", c("source", "step", "gene_id", "gene", "clinvar_id"), drop = FALSE])
    for (j in seq_len(nrow(x))) add(i, x$source[[j]], x$step[[j]], x[j, ])
    for (step in interp_summary_scope()[[summary$question[[i]]]]) if (states[[step]] != "after")
      add(i, paste0("pending_", step), step, list(source = "", gene_id = "", gene = "", clinvar_id = ""), TRUE)
  }
  groups <- c(local = "Local", public = "Public query", spliceai = "SpliceAI", ld = "LD scan")
  node_y <- numeric(); headers <- list(); y <- 36
  for (step in names(groups)) {
    ids <- names(nodes)[vapply(nodes, function(n) n$step == step, logical(1))]
    if (!length(ids)) next
    headers[[step]] <- y; y <- y + 25
    for (id in ids) { node_y[[id]] <- y; y <- y + 38 }
    y <- y + 14
  }
  source_ids <- stats::setNames(paste0("s", seq_along(nodes)), names(nodes))
  height <- max(540, y + 12); qy <- seq(65, height - 45, length.out = 7L)
  S <- function(name, ...) htmltools::tag(name, list(...))
  title <- function(text) S("title", text)
  links <- lapply(edges, function(e) S("path", class = if (e$pending) "regshiny-map-edge pending" else "regshiny-map-edge contributed",
    `data-question` = summary$question[[e$question]], `data-step` = nodes[[e$node]]$step,
    `data-question-node` = paste0("q", e$question), `data-source-node` = source_ids[[e$node]],
    d = sprintf("M 330 %.1f C 490 %.1f, 550 %.1f, 710 %.1f", qy[[e$question]], qy[[e$question]], node_y[[e$node]], node_y[[e$node]]),
    title(if (e$pending) "Step not completed" else "Source contributes evidence or context")))
  left <- lapply(seq_len(nrow(summary)), function(i) S("g", class = paste("regshiny-map-question", interpretation_level_class(summary$level[[i]])),
    `data-question-node` = paste0("q", i), tabindex = "0",
    title(paste(summary$question[[i]], summary$level[[i]], sep = ": ")), S("rect", x = 14, y = qy[[i]] - 24, width = 316, height = 48, rx = 5),
    S("text", x = 26, y = qy[[i]] - 4, summary$question[[i]]), S("text", x = 26, y = qy[[i]] + 15, class = "regshiny-map-level", summary$level[[i]])))
  right <- lapply(names(nodes), function(id) {
    n <- nodes[[id]]; info <- interpretation_source_info(n$key)
    label <- if (n$pending) interp_step_label(n$step, states[[n$step]]) else paste0(info$label, if (info$kind == "external") "" else paste0(" · ", info$kind))
    lines <- head(strwrap(label, width = 37), 2L)
    node <- S("g", class = "regshiny-map-source", `data-source` = n$key, `data-step` = n$step, `data-source-node` = source_ids[[id]], tabindex = "0",
      title(paste(label, interpretation_source_note(n$key))), S("rect", x = 710, y = node_y[[id]] - 16, width = 390, height = 32, rx = 4),
      lapply(seq_along(lines), function(k) S("text", x = 722, y = node_y[[id]] - 3 + (k - 1) * 15, lines[[k]])))
    ids <- interpretation_item_identifiers(n$item, identifiers)
    url <- if (!n$pending) interp_source_url(n$key, ids) else NULL
    if (is.null(url)) node else S("a", href = url, target = "_blank", rel = "noopener noreferrer", title = interpretation_source_note(n$key), node)
  })
  # Per-node selectors use generated numeric IDs, never public identifiers or prose.
  # :has() keeps hover/focus tracing entirely in CSS (also for linked source nodes).
  highlight <- function(kind, id) {
    node <- if (kind == "question") ".regshiny-map-question" else ".regshiny-map-source"
    attr <- paste0("data-", kind, "-node")
    paste0('.regshiny-evidence-map:has(', node, '[', attr, '="', id, '"]:is(:hover,:focus-visible,:focus-within)) ',
      '.regshiny-map-edge[', attr, '="', id, '"] {opacity:1;stroke:#0b5cad;stroke-width:3;}')
  }
  highlight_css <- paste(c(vapply(paste0("q", seq_len(nrow(summary))), function(id) highlight("question", id), character(1)),
    vapply(unname(source_ids), function(id) highlight("source", id), character(1))), collapse = "\n")
  S("svg", class = "regshiny-evidence-map", xmlns = "http://www.w3.org/2000/svg", viewBox = paste(0, 0, 1120, height), role = "img",
    `aria-label` = "Evidence questions and their sources, grouped by requested step", title("Solid: source contribution. Dashed: step not completed."),
    S("style", htmltools::HTML(highlight_css)), links, left, lapply(names(headers), function(step) S("text", class = "regshiny-map-group", x = 710, y = headers[[step]], groups[[step]])), right)
}
