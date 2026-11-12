# Variant interpretation tab: what the selected DNA variant does and what is
# known about it. Local lookups run on their own; public databases and SpliceAI
# run only when the user asks. Nothing here changes RNA counts, Model 1 or a
# saved result.

interpretation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("header"), role = "status", `aria-live` = "polite"),
    shiny::div(class = "regshiny-actions",
      bslib::input_task_button(ns("public"), "Query public databases"),
      bslib::input_task_button(ns("spliceai"), "Run SpliceAI (about 1-2 min)"),
      shiny::downloadButton(ns("download_json"), "Interpretation JSON"),
      shiny::downloadButton(ns("download_summary"), "Evidence summary CSV")),
    shiny::helpText("Local sources are read automatically. A public query sends only the variant's chromosome, position and alleles, and the Ensembl IDs of the genes it overlaps, to Open Targets, GTEx and MyVariant.info. No sample, genotype or cohort information is sent."),
    shiny::uiOutput(ns("public_status"), role = "status", `aria-live` = "polite"),
    shiny::h5("Evidence summary"),
    shiny::uiOutput(ns("summary")),
    shiny::tags$details(class = "regshiny-definitions", shiny::tags$summary("How the levels are assigned"),
      shiny::tags$ul(
        shiny::tags$li(shiny::strong("Strong"), ": direct or fine-mapped evidence that this variant is responsible (a protein-truncating or essential splice-site change on the reference transcript; an RNA association together with a prediction or a posterior probability of at least 0.5; ClinVar pathogenic without a benign assertion; a GWAS credible set with posterior probability of at least 0.5)."),
        shiny::tags$li(shiny::strong("Moderate"), ": an association or a prediction, but not both. A QTL alone can be a neighbour in linkage disequilibrium."),
        shiny::tags$li(shiny::strong("Limited"), ": weak or indirect support only."),
        shiny::tags$li(shiny::strong("None found"), ": the sources were consulted and report nothing. ", shiny::strong("Not assessed"), ": the sources needed were unavailable or not queried."),
        shiny::tags$li(shiny::strong("Context"), ": describes the gene or tissue, not evidence about the variant.")),
      shiny::p("The levels rank evidence for follow-up. They are not a clinical classification, and a QTL p-value is not a causal claim. Full rules: docs/METRICS.md.")),
    bslib::navset_card_tab(id = ns("detail"),
      bslib::nav_panel("Transcripts and protein",
        shiny::h5("Most severe consequence per gene"),
        DT::DTOutput(ns("gene_summary")),
        shiny::h5("Every transcript"),
        DT::DTOutput(ns("consequences")),
        shiny::downloadButton(ns("download_consequences"), "Consequences CSV"),
        shiny::p("Consequences are computed here from the selected GENCODE release and the reference sequence, with Sequence Ontology terms and Ensembl VEP's severity order. donor_distance and acceptor_distance count into the intron from its 5' and 3' ends (1 = first or last intronic base); negative values are exonic bases counted from the exon edge. Insertions, deletions and multi-base changes are classified by what they overlap and get no protein change."),
        shiny::uiOutput(ns("predictors"))),
      bslib::nav_panel("Splicing",
        shiny::uiOutput(ns("spliceai_status"), role = "status", `aria-live` = "polite"),
        DT::DTOutput(ns("spliceai_table")),
        shiny::p("SpliceAI gives the change in predicted splice-site probability (0 to 1) caused by the ALT allele, for the strongest acceptor gain, acceptor loss, donor gain and donor loss within 500 bp. Common reading: 0.2 permissive, 0.5 recommended, 0.8 high precision. site_position is the first or last intronic base next to that site, the coordinate used for sites in Genotype association; observed_in_cohort joins the bin-ratio row of that site when this result has one."),
        shiny::h5("RNA evidence in this result"),
        shiny::uiOutput(ns("rna_evidence"))),
      bslib::nav_panel("QTL by tissue",
        shiny::h5("TOPMed cis-QTL in six tissues and cell types"),
        shiny::uiOutput(ns("topmed_status")),
        shiny::plotOutput(ns("context_plot"), height = 420),
        DT::DTOutput(ns("heterogeneity")),
        shiny::p("One row per gene (eQTL) or intron (sQTL) that is significant in at least one context. Slopes are per ALT allele on normalised phenotypes, so they compare across contexts but not with Model 1. A pair is significant when its nominal p is below that gene's permutation threshold (FDR 5%). q_p_value and i_squared describe how much the slopes differ between contexts (Cochran's Q); the contexts share donors, so read them as descriptive. A context with few samples can miss an effect that is there."),
        shiny::h5("Tested pairs"),
        DT::DTOutput(ns("topmed_pairs")),
        shiny::p("Kept: every significant pair, every pair in a SuSiE credible set (pip), pairs with nominal p below 0.001, the three smallest p per context, and every context's row of a phenotype that is significant somewhere. n_estimated is derived from the minor-allele count and frequency."),
        shiny::h5("FHS whole blood"),
        shiny::uiOutput(ns("fhs_status")),
        DT::DTOutput(ns("fhs_table")),
        shiny::h5("GTEx v10 (public query)"),
        DT::DTOutput(ns("gtex_qtl")),
        shiny::plotOutput(ns("expression_plot"), height = 430)),
      bslib::nav_panel("Disease and traits",
        shiny::h5("ClinVar"),
        DT::DTOutput(ns("clinvar")),
        shiny::h5("GWAS credible sets containing this variant (Open Targets)"),
        shiny::uiOutput(ns("gwas_status")),
        DT::DTOutput(ns("gwas")),
        shiny::p("posterior_probability is this variant's probability of being the causal one within that credible set; a variant can sit in a credible set with a small probability. l2g_gene is Open Targets' top locus-to-gene prediction for the signal."),
        shiny::h5("Other clinical and curated evidence"),
        DT::DTOutput(ns("clinical")),
        shiny::h5("Pharmacogenomics"),
        DT::DTOutput(ns("pgx"))),
      bslib::nav_panel("Protein, genes and pathways",
        shiny::h5("Protein and molecular QTL credible sets (Open Targets)"),
        DT::DTOutput(ns("ot_qtl")),
        shiny::h5("Genes the variant falls in"),
        DT::DTOutput(ns("genes")),
        shiny::h5("Reactome pathways"),
        DT::DTOutput(ns("pathways")),
        shiny::h5("Gene Ontology"),
        DT::DTOutput(ns("go")),
        shiny::h5("Diseases most associated with these genes"),
        DT::DTOutput(ns("gene_diseases")),
        shiny::p("Pathways, GO terms and gene-level disease associations describe the gene. They do not show that this variant changes them.")),
      bslib::nav_panel("Sources",
        shiny::h5("Local sources"),
        DT::DTOutput(ns("local_sources")),
        shiny::h5("Population frequency of the ALT allele"),
        DT::DTOutput(ns("frequencies")),
        shiny::verbatimTextOutput(ns("provenance")))))
}

interpretation_level_class <- function(level) switch(level, Strong = "regshiny-level-strong", Moderate = "regshiny-level-moderate",
  Limited = "regshiny-level-limited", `None found` = "regshiny-level-none", Context = "regshiny-level-context", "regshiny-level-unknown")

interpretation_summary_cards <- function(summary) {
  shiny::div(class = "regshiny-evidence", lapply(seq_len(nrow(summary)), function(i) {
    evidence <- strsplit(summary$evidence[[i]], " | ", fixed = TRUE)[[1L]]
    shiny::div(class = "regshiny-evidence-row",
      shiny::div(class = "regshiny-evidence-head",
        shiny::strong(summary$question[[i]]),
        shiny::span(class = paste("regshiny-level", interpretation_level_class(summary$level[[i]])), summary$level[[i]])),
      shiny::p(summary$finding[[i]]),
      if (length(evidence)) shiny::tags$ul(lapply(evidence, shiny::tags$li)),
      if (nzchar(summary$sources[[i]])) shiny::p(class = "regshiny-evidence-source", paste("Sources:", summary$sources[[i]])))
  }))
}

# Slopes of the leading phenotypes in every context, with 95% intervals.
interpretation_plot_contexts <- function(pairs, heterogeneity, max_panels = 4L) {
  if (is.null(pairs) || !nrow(pairs) || is.null(heterogeneity) || !nrow(heterogeneity)) {
    graphics::plot.new(); graphics::text(.5, .5, "No significant TOPMed cis-QTL for this variant"); return(invisible(NULL))
  }
  h <- utils::head(heterogeneity, max_panels)
  contexts <- unique(pairs$context)
  old <- graphics::par(mfrow = c(1L, nrow(h)), mar = c(4.2, if (nrow(h) > 2L) 6.6 else 8.4, 3.6, 1), cex = .85)
  on.exit(graphics::par(old))
  key <- paste(pairs$qtl_type, sub("^([^:]+:[0-9]+:[0-9]+):.*:(ENSG[^:]*)$", "\\1:\\2", pairs$phenotype))
  for (i in seq_len(nrow(h))) {
    x <- pairs[key == paste(h$qtl_type[[i]], h$phenotype[[i]]), , drop = FALSE]
    x <- x[match(contexts, x$context), , drop = FALSE]; y <- rev(seq_along(contexts))
    lo <- x$slope_alt - 1.96 * x$slope_se; hi <- x$slope_alt + 1.96 * x$slope_se
    lim <- range(c(0, lo, hi), na.rm = TRUE); if (!all(is.finite(lim))) lim <- c(-1, 1)
    graphics::plot(NA, xlim = lim + c(-.05, .05) * diff(lim), ylim = c(.5, length(contexts) + .5), yaxt = "n", ylab = "",
      xlab = "Slope per ALT allele (95% CI)")
    graphics::axis(2, at = y, labels = contexts, las = 1, cex.axis = .85)
    graphics::abline(v = 0, lty = 2, col = "#888888")
    graphics::segments(lo, y, hi, y, col = "#1F4E79", lwd = 1.6)
    significant <- x$significant %in% TRUE
    graphics::points(x$slope_alt, y, pch = ifelse(significant, 19, 21), col = "#1F4E79", bg = "white", cex = 1.3)
    missing <- is.na(x$slope_alt)
    if (any(missing)) graphics::text(mean(lim), y[missing], "not tested", cex = .8, col = "#777777")
    gene <- if (nzchar(h$gene[[i]])) h$gene[[i]] else h$gene_id[[i]]
    intron <- if (h$qtl_type[[i]] == "sQTL") sub("^[^:]+:([0-9]+):([0-9]+):.*$", "intron \\1-\\2", h$phenotype[[i]]) else "expression"
    graphics::title(main = paste0(h$qtl_type[[i]], " · ", gene), cex.main = .95, line = 2.1)
    graphics::mtext(intron, side = 3, line = .9, cex = .68)
    graphics::mtext("filled: significant", side = 3, line = .1, cex = .6, col = "#555555")
  }
  invisible(h)
}

interpretation_plot_expression <- function(expression, gene) {
  x <- expression[expression$gene == gene & is.finite(expression$median_tpm), , drop = FALSE]
  if (!nrow(x)) { graphics::plot.new(); graphics::text(.5, .5, "Run the public query to see GTEx expression by tissue"); return(invisible(NULL)) }
  x <- x[order(x$median_tpm), , drop = FALSE]
  old <- graphics::par(mar = c(4.2, 15, 2.6, 1), cex = .72); on.exit(graphics::par(old))
  graphics::barplot(x$median_tpm, names.arg = x$tissue, horiz = TRUE, las = 1, border = NA,
    col = ifelse(x$tissue == "Whole Blood", "#B7410E", "#7A9CC6"), xlab = "Median TPM (GTEx v10)",
    main = paste0(gene, " expression by tissue · whole blood in orange"))
  invisible(x)
}

# variant: reactive list(chrom, pos1, ref, alt, build) or NULL.
# bin_ratio / junctions: reactives returning the open result's Model 1 tables (or NULL).
interpretation_server <- function(id, variant, annotation_resources, annotation_file, reference, splice_file,
                                  interpretation_file, bin_ratio = NULL, junctions = NULL, is_active = shiny::reactive(TRUE)) {
  shiny::moduleServer(id, function(input, output, session) {
    current_annotation <- if (is.function(annotation_resources)) annotation_resources else {
      fixed <- annotation_resources; function() fixed }
    alleles <- shiny::reactive({
      shiny::req(is_active())
      v <- variant(); shiny::req(!is.null(v))
      alt <- strsplit(as.character(v$alt[[1L]]), ",", fixed = TRUE)[[1L]]
      list(chrom = as.character(v$chrom[[1L]]), pos1 = as.integer(v$pos1[[1L]]), ref = as.character(v$ref[[1L]]), alt = alt,
        build = if (is.null(v$build)) "GRCh38" else as.character(v$build[[1L]]))
    })
    selected_alt <- shiny::reactive({
      a <- alleles()
      if (length(a$alt) > 1L && length(input$alt) == 1L && input$alt %in% a$alt) input$alt else a$alt[[1L]]
    })
    key <- shiny::reactive({
      a <- alleles(); res <- current_annotation()
      release <- if (is.list(res) && !inherits(res, "error") && length(res$release) == 1L) res$release else "none"
      paste(a$build, a$chrom, a$pos1, a$ref, selected_alt(), release, sep = ":")
    })
    variant_key <- shiny::reactive(sub(":[^:]*$", "", key()))

    local_store <- shiny::reactiveVal(NULL)
    local_task <- shiny::ExtendedTask$new(function(key, allele, annotation, annotation_file, reference, splice_file, interpretation_file) {
      promises::future_promise({
        value <- tryCatch({
          e <- new.env(parent = globalenv()); sys.source(interpretation_file, envir = e)
          v <- e$interp_variant(allele$chrom, allele$pos1, allele$ref, allele$alt)
          out <- list(variant = v, consequences = e$interp_empty_consequences(), annotation_error = NULL, reference_base = NA_character_,
            release = NULL, introns = NULL)
          fetch <- function(chrom, start1, end1) "N"
          if (!is.null(reference) && !inherits(reference, "error") && !is.null(splice_file)) {
            s <- new.env(parent = globalenv()); sys.source(splice_file, envir = s)
            fetch <- function(chrom, start1, end1) s$splice_fetch_reference(reference, chrom, start1, end1)$sequence
            out$reference_base <- tryCatch(fetch(v$chrom, v$pos1, v$pos1 + nchar(v$ref) - 1L), error = function(err) NA_character_)
          }
          names <- NULL
          if (is.null(annotation) || inherits(annotation, "error") || is.null(annotation_file)) {
            out$annotation_error <- if (inherits(annotation, "error")) conditionMessage(annotation) else "GENCODE annotation is unavailable."
          } else {
            a <- new.env(parent = globalenv()); sys.source(annotation_file, envir = a)
            q <- tryCatch(a$query_annotation(annotation, v$chrom, max(1L, v$pos1 - 5000L), v$pos1 + 5000L, allele$build), error = function(err) err)
            if (inherits(q, "error")) out$annotation_error <- conditionMessage(q) else {
              out$consequences <- e$interp_consequences(v, q, fetch)
              out$release <- q$metadata$release
              out$introns <- unique(q$introns[, c("strand", "donor1", "acceptor1")])
            }
          }
          out$gene_summary <- e$interp_gene_summary(out$consequences)
          out$local <- if (identical(allele$build, "GRCh38")) e$interp_local_evidence(v) else NULL
          if (!is.null(out$local) && !inherits(out$local$topmed, "error") && nrow(out$local$topmed$pairs) &&
              !is.null(annotation) && !inherits(annotation, "error") && !is.null(annotation_file)) {
            a <- new.env(parent = globalenv()); sys.source(annotation_file, envir = a)
            names <- tryCatch(a$annotation_gene_names(annotation, out$local$topmed$pairs$gene_id), error = function(err) NULL)
            out$local$topmed$pairs <- e$interp_name_genes(out$local$topmed$pairs, names)
          }
          out
        }, error = function(err) err)
        list(key = key, value = value)
      }, seed = TRUE)
    })
    shiny::observeEvent(key(), {
      a <- alleles()
      allele <- list(chrom = a$chrom, pos1 = a$pos1, ref = a$ref, alt = selected_alt(), build = a$build)
      local_task$invoke(key(), allele, current_annotation(), annotation_file, reference, splice_file, interpretation_file)
    })
    shiny::observeEvent(local_task$status(), {
      if (local_task$status() == "success") local_store(local_task$result())
    }, ignoreInit = TRUE)
    local_now <- shiny::reactive({
      x <- local_store()
      shiny::validate(shiny::need(!is.null(x) && identical(x$key, key()), "Reading annotation and local sources for this variant..."))
      shiny::validate(shiny::need(!inherits(x$value, "error"), paste("Interpretation failed:", if (inherits(x$value, "error")) conditionMessage(x$value) else "")))
      x$value
    })

    # Public databases: only on request, kept per variant for this session.
    public_store <- shiny::reactiveVal(list())
    public_task <- shiny::ExtendedTask$new(function(key, variant, gene_ids, interpretation_file) {
      promises::future_promise({
        e <- new.env(parent = globalenv()); sys.source(interpretation_file, envir = e)
        list(key = key, value = tryCatch(e$interp_public_evidence(variant, gene_ids), error = function(err) err))
      }, seed = TRUE)
    })
    bslib::bind_task_button(public_task, "public")
    shiny::observeEvent(input$public, {
      x <- local_now()
      genic <- x$consequences[!x$consequences$most_severe %in% c("upstream_gene_variant", "downstream_gene_variant"), , drop = FALSE]
      genes <- unique(if (nrow(genic)) genic$gene_id else utils::head(x$consequences$gene_id, 3L))
      public_task$invoke(variant_key(), x$variant, genes, interpretation_file)
    })
    shiny::observeEvent(public_task$status(), {
      if (public_task$status() != "success") return()
      r <- public_task$result(); store <- public_store(); store[[r$key]] <- r$value; public_store(store)
    }, ignoreInit = TRUE)
    public_now <- shiny::reactive({
      x <- public_store()[[variant_key()]]
      if (is.null(x) || inherits(x, "error")) NULL else x
    })

    spliceai_store <- shiny::reactiveVal(list())
    spliceai_task <- shiny::ExtendedTask$new(function(key, variant, fasta, interpretation_file) {
      promises::future_promise({
        e <- new.env(parent = globalenv()); sys.source(interpretation_file, envir = e)
        list(key = key, value = tryCatch(e$interp_spliceai(variant, fasta), error = function(err) err))
      }, seed = TRUE)
    })
    bslib::bind_task_button(spliceai_task, "spliceai")
    shiny::observeEvent(input$spliceai, {
      x <- local_now()
      fasta <- if (!is.null(reference) && !inherits(reference, "error")) reference$fasta else ""
      spliceai_task$invoke(variant_key(), x$variant, fasta, interpretation_file)
    })
    shiny::observeEvent(spliceai_task$status(), {
      if (spliceai_task$status() != "success") return()
      r <- spliceai_task$result(); store <- spliceai_store(); store[[r$key]] <- r$value; spliceai_store(store)
    }, ignoreInit = TRUE)
    spliceai_now <- shiny::reactive(spliceai_store()[[variant_key()]])

    backend <- new.env(parent = globalenv()); sys.source(interpretation_file, envir = backend)
    rna_now <- shiny::reactive({
      get_table <- function(f) if (is.null(f)) NULL else tryCatch(f(), error = function(e) NULL)
      backend$interp_rna_evidence(get_table(bin_ratio), get_table(junctions))
    })
    summary_now <- shiny::reactive({
      x <- local_now()
      backend$interp_evidence_summary(x$consequences, x$local, public_now(), spliceai_now(), rna_now())
    })

    table_options <- list(pageLength = 10, scrollX = TRUE, dom = "ftip")
    show_table <- function(expr, digits = 3L, empty = "No rows.") DT::renderDT({
      x <- expr()
      shiny::validate(shiny::need(is.data.frame(x) && nrow(x) > 0L, empty))
      numeric <- vapply(x, function(col) is.numeric(col) && !is.integer(col), logical(1))
      x[numeric] <- lapply(x[numeric], function(col) signif(col, digits))
      x
    }, rownames = FALSE, selection = "none", options = table_options)

    output$header <- shiny::renderUI({
      a <- alleles(); x <- local_now(); v <- x$variant
      known <- if (!is.null(x$local) && !inherits(x$local$known, "error")) x$local$known else NULL
      cadd <- if (!is.null(x$local) && !inherits(x$local$cadd, "error")) x$local$cadd else NULL
      af <- if (!is.null(known) && nrow(known$frequencies)) known$frequencies$allele_frequency[match("gnomADg", known$frequencies$population)] else NA_real_
      mismatch <- !is.na(x$reference_base) && !identical(x$reference_base, v$ref)
      shiny::tagList(
        if (length(a$alt) > 1L) shiny::selectInput(session$ns("alt"), "ALT allele to interpret", choices = a$alt, selected = selected_alt()),
        shiny::div(class = "regshiny-context",
          shiny::strong(paste0(v$id, if (!is.null(known) && nzchar(known$rs_id)) paste0(" · ", known$rs_id) else "")),
          shiny::p(paste0(v$type, if (!is.null(x$release)) paste0(" · GENCODE v", x$release) else "",
            if (is.finite(af)) paste0(" · gnomAD genomes ALT frequency ", signif(af, 3)) else "",
            if (!is.null(cadd) && isTRUE(cadd$found)) paste0(" · CADD PHRED ", cadd$phred) else "")),
          if (mismatch) shiny::p(class = "text-danger", paste0("REF allele ", v$ref, " does not match the reference sequence (", x$reference_base,
            ") at this position. Check the assembly before reading the consequences.")),
          if (!is.null(x$annotation_error)) shiny::p(class = "text-danger", paste("Transcript consequences unavailable:", x$annotation_error)),
          if (!identical(a$build, "GRCh38")) shiny::p("Local and public sources are GRCh38 only and were not queried for this assembly.")))
    })
    output$public_status <- shiny::renderUI({
      s <- public_task$status()
      if (s == "running") return(shiny::div(class = "alert alert-info", "Querying Open Targets, GTEx and MyVariant.info..."))
      raw <- public_store()[[variant_key()]]
      if (s == "error" || inherits(raw, "error")) return(shiny::div(class = "alert alert-danger", "Public query failed: ",
        if (inherits(raw, "error")) conditionMessage(raw) else "the request did not complete."))
      if (is.null(raw)) return(NULL)
      failed <- Filter(function(n) inherits(raw[[n]], "error"), c("open_targets", "genes", "gtex", "myvariant"))
      labels <- c(open_targets = "Open Targets (variant)", genes = "Open Targets (genes)", gtex = "GTEx", myvariant = "MyVariant.info")
      shiny::div(class = if (length(failed)) "alert alert-warning" else "alert alert-success",
        paste0("Public sources queried at ", raw$queried_at, "."),
        if (length(failed)) paste0(" Unavailable: ", paste(sprintf("%s (%s)", labels[failed],
          vapply(failed, function(n) conditionMessage(raw[[n]]), character(1))), collapse = "; "), ".") else "",
        if (!inherits(raw$open_targets, "error") && !isTRUE(raw$open_targets$found)) " Open Targets has no record of this variant." else "",
        if (!inherits(raw$gtex, "error") && identical(raw$gtex$variant_in_gtex, FALSE)) " GTEx did not test this variant." else "")
    })
    output$summary <- shiny::renderUI(interpretation_summary_cards(summary_now()))
    output$gene_summary <- show_table(function() local_now()$gene_summary, empty = "No transcript within 5 kb of the variant.")
    output$consequences <- show_table(function() {
      x <- local_now()$consequences
      x[, setdiff(names(x), "most_severe"), drop = FALSE]
    }, empty = "No transcript within 5 kb of the variant.")
    output$predictors <- shiny::renderUI({
      p <- public_now()
      if (is.null(p) || inherits(p$open_targets, "error") || !nrow(p$open_targets$predictors)) return(NULL)
      x <- p$open_targets$predictors
      shiny::tagList(shiny::h5("In-silico predictors (Open Targets)"),
        shiny::tags$ul(lapply(seq_len(nrow(x)), function(i) shiny::tags$li(paste0(x$method[[i]], ": ",
          if (nzchar(x$assessment[[i]])) paste0(x$assessment[[i]], ", ") else "", "score ", signif(x$score[[i]], 3),
          if (nzchar(x$gene[[i]])) paste0(" (", x$gene[[i]], ")") else "")))))
    })

    spliceai_sites <- shiny::reactive({
      s <- spliceai_now(); x <- local_now()
      shiny::req(!is.null(s), !inherits(s, "error"))
      strand <- x$consequences$strand[match(s$scores$gene, x$consequences$gene_name)]
      rows <- backend$interp_spliceai_sites(s$scores, strand, x$introns)
      b <- if (is.null(bin_ratio)) NULL else tryCatch(bin_ratio(), error = function(e) NULL)
      rows$observed_in_cohort <- ""
      if (is.data.frame(b) && nrow(b) && nrow(rows)) {
        k <- match(paste(rows$site_type, rows$site_position), paste(b$splice_site_type, b$splice_site_position))
        rows$observed_in_cohort <- ifelse(is.na(k), "", sprintf("bin ratio %+.3f per ALT allele, permutation p %s (%s)", b$beta[k],
          format(signif(b$p_permutation[k], 2)), b$status[k]))
      }
      rows
    })
    output$spliceai_status <- shiny::renderUI({
      s <- spliceai_task$status(); r <- spliceai_now()
      if (s == "running") return(shiny::div(class = "alert alert-info", "SpliceAI is running. Loading its models takes most of the time; the rest of the app stays usable."))
      if (inherits(r, "error")) return(shiny::div(class = "alert alert-danger", conditionMessage(r)))
      if (is.null(r)) return(shiny::p("SpliceAI has not been run for this variant. Use the button above; it needs no network."))
      if (!nrow(r$scores)) return(shiny::div(class = "alert alert-warning", "SpliceAI returned no score: the variant is outside the genes in its annotation."))
      shiny::div(class = "alert alert-success", paste0(r$model, ", window ", r$distance, " bp, run at ", r$ran_at, "."))
    })
    output$spliceai_table <- show_table(spliceai_sites, empty = "No SpliceAI result.")
    output$rna_evidence <- shiny::renderUI({
      r <- rna_now()
      if (!isTRUE(r$available)) return(shiny::p("No genotype comparison is open, so there is no cohort RNA evidence for this variant yet."))
      if (!length(r$details)) return(shiny::p("Model 1 fitted no site or junction with a permutation p below the Bonferroni cut-off for this result."))
      shiny::tagList(shiny::tags$ul(lapply(r$details, shiny::tags$li)),
        shiny::p("Listed: fitted rows whose permutation p is below 0.05 divided by the number of fitted rows of that table. Sites flagged NO_DEPTH_DROP_AT_SITE are set aside. See Genotype association for every row."))
    })

    topmed_now <- shiny::reactive({
      x <- local_now()
      shiny::validate(shiny::need(!is.null(x$local) && !inherits(x$local$topmed, "error") && isTRUE(x$local$topmed$available),
        "The TOPMed cis-QTL files are not readable from this session."))
      x$local$topmed
    })
    heterogeneity_now <- shiny::reactive(backend$interp_context_heterogeneity(topmed_now()$pairs))
    output$topmed_status <- shiny::renderUI({
      t <- topmed_now(); ctx <- t$contexts
      untested <- unique(ctx$context[!ctx$tested])
      shiny::p(paste0("Significant pairs: ", sum(ctx$significant_pairs), " across ", sum(ctx$significant_pairs > 0), " of ", sum(ctx$tested),
        " tested context-and-type combinations.", if (length(untested)) paste0(" Not tested (absent or below the frequency cut-off): ",
          paste(untested, collapse = ", "), ".") else ""))
    })
    output$context_plot <- shiny::renderPlot(interpretation_plot_contexts(topmed_now()$pairs, heterogeneity_now()))
    output$heterogeneity <- show_table(heterogeneity_now, empty = "This variant is not a significant cis-QTL in any TOPMed context.")
    output$topmed_pairs <- show_table(function() topmed_now()$pairs, empty = "This variant was not tested in the TOPMed cis-QTL runs.")
    fhs_now <- shiny::reactive({
      x <- local_now()
      shiny::validate(shiny::need(!is.null(x$local) && !inherits(x$local$fhs, "error"), "The FHS QTL files are not readable from this session."))
      x$local$fhs
    })
    output$fhs_status <- shiny::renderUI({
      f <- fhs_now()
      shiny::tagList(shiny::p(if (f$sqtl_available) f$sqtl_note else "Isoform-ratio sQTL files are not available."),
        shiny::p(if (f$eqtl_available) f$eqtl_note else "The position-sorted gene eQTL copy is not available (see scripts/build_qtl_index.sh)."))
    })
    output$fhs_table <- show_table(function() {
      f <- fhs_now()
      rbind(if (nrow(f$sqtl)) cbind(qtl_type = "sQTL (isoform ratio)", f$sqtl) else NULL,
            if (nrow(f$eqtl)) cbind(qtl_type = "eQTL (gene)", f$eqtl) else NULL)
    }, empty = "No FHS whole-blood QTL row for this variant at the stored thresholds.")
    public_part <- function(source, field, message = "Run the public query to fill this table.") function() {
      p <- public_now()
      shiny::validate(shiny::need(!is.null(p), message))
      shiny::validate(shiny::need(!inherits(p[[source]], "error"), "This source did not answer."))
      p[[source]][[field]]
    }
    output$gtex_qtl <- show_table(public_part("gtex", "qtl"), empty = "No significant GTEx eQTL or sQTL for this variant.")
    expression_gene <- shiny::reactive({
      p <- public_now(); shiny::req(!is.null(p), !inherits(p$gtex, "error"), nrow(p$gtex$expression) > 0L)
      p$gtex$expression$gene[[1L]]
    })
    output$expression_plot <- shiny::renderPlot({
      p <- public_now()
      if (is.null(p) || inherits(p$gtex, "error") || !nrow(p$gtex$expression)) return(interpretation_plot_expression(backend$interp_gtex_empty()$expression, ""))
      interpretation_plot_expression(p$gtex$expression, expression_gene())
    })
    output$clinvar <- show_table(function() {
      p <- public_now(); x <- local_now()
      if (!is.null(p) && !inherits(p$myvariant, "error") && nrow(p$myvariant$clinvar)) return(p$myvariant$clinvar)
      known <- if (!is.null(x$local) && !inherits(x$local$known, "error")) x$local$known else NULL
      if (!is.null(known) && nzchar(known$clinvar_significance))
        return(data.frame(source = "ClinVar 2023-10 (Ensembl VEP 112 cache)", clinical_significance = known$clinvar_significance, stringsAsFactors = FALSE))
      NULL
    }, empty = "No ClinVar record for this allele in the sources consulted.")
    output$gwas_status <- shiny::renderUI({
      p <- public_now()
      if (is.null(p) || inherits(p$open_targets, "error") || !isTRUE(p$open_targets$found)) return(NULL)
      total <- p$open_targets$summary$gwas_total; shown <- nrow(p$open_targets$gwas)
      shiny::p(paste0(total, " credible set(s)", if (shown < total) paste0("; the first ", shown, " are listed") else "", "."))
    })
    output$gwas <- show_table(public_part("open_targets", "gwas"), empty = "Not in any GWAS credible set in Open Targets.")
    output$clinical <- show_table(public_part("open_targets", "clinical"), empty = "No curated clinical evidence in Open Targets.")
    output$pgx <- show_table(public_part("open_targets", "pharmacogenomics"), empty = "No pharmacogenomic annotation.")
    output$ot_qtl <- show_table(public_part("open_targets", "qtl"), empty = "Not in any molecular-QTL credible set in Open Targets.")
    output$genes <- show_table(public_part("genes", "genes"), empty = "No gene record.")
    output$pathways <- show_table(public_part("genes", "pathways"), empty = "No Reactome pathway.")
    output$go <- show_table(public_part("genes", "go"), empty = "No GO annotation.")
    output$gene_diseases <- show_table(public_part("genes", "diseases"), empty = "No associated disease.")
    output$local_sources <- DT::renderDT({
      x <- local_now()
      if (is.null(x$local)) backend$interp_local_status() else x$local$status
    }, rownames = FALSE, selection = "none", options = list(dom = "t", paging = FALSE, scrollX = TRUE))
    output$frequencies <- show_table(function() {
      x <- local_now()
      if (is.null(x$local) || inherits(x$local$known, "error")) NULL else x$local$known$frequencies
    }, empty = "No frequency record in the Ensembl VEP cache for this allele.")
    output$provenance <- shiny::renderPrint({
      x <- local_now(); p <- public_store()[[variant_key()]]; s <- spliceai_now()
      cat("Variant:", x$variant$id, "\nGENCODE release:", if (is.null(x$release)) "unavailable" else x$release, "\n")
      cat("Reference base(s) at the variant:", x$reference_base, "\n")
      if (!is.null(x$local)) for (n in c("cadd", "known", "fhs", "topmed")) if (inherits(x$local[[n]], "error"))
        cat("Local source", n, "failed:", conditionMessage(x$local[[n]]), "\n")
      cat("Public query:", if (is.null(p)) "not run" else if (inherits(p, "error")) conditionMessage(p) else
        paste0("run at ", p$queried_at, "; sent variant ", p$variant, " and gene IDs ", paste(p$gene_ids, collapse = ", ")), "\n")
      cat("SpliceAI:", if (is.null(s)) "not run" else if (inherits(s, "error")) conditionMessage(s) else paste(s$model, "at", s$ran_at), "\n")
    })
    export <- function() {
      x <- local_now(); p <- public_now(); s <- spliceai_now()
      plain <- function(v) if (inherits(v, "error")) list(error = conditionMessage(v)) else v
      list(variant = x$variant, gencode_release = x$release, summary = summary_now(), gene_summary = x$gene_summary,
        consequences = x$consequences, local = if (is.null(x$local)) NULL else lapply(x$local, plain),
        public = if (is.null(p)) NULL else lapply(p, plain), spliceai = if (is.null(s)) NULL else plain(s), rna_evidence = rna_now())
    }
    output$download_json <- shiny::downloadHandler(function() paste0("interpretation_", gsub(":", "_", local_now()$variant$id), ".json"),
      function(file) writeLines(jsonlite::toJSON(export(), dataframe = "rows", auto_unbox = TRUE, na = "null", null = "null", digits = NA, pretty = TRUE), file))
    output$download_summary <- shiny::downloadHandler("evidence_summary.csv", function(file)
      utils::write.csv(summary_now(), file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"))
    output$download_consequences <- shiny::downloadHandler("transcript_consequences.csv", function(file)
      utils::write.csv(local_now()$consequences, file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"))
    invisible(list(summary = summary_now, local = local_now, public = public_now, spliceai = spliceai_now,
      tasks = list(local = local_task, public = public_task, spliceai = spliceai_task)))
  })
}
