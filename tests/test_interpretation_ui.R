# Synthetic-only UI tests for the Variant interpretation tab. No controlled data,
# no network: annotation, reference, public sources and SpliceAI are fixtures.
needed <- c("shiny", "bslib", "DT", "future", "promises", "later", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("INTERPRETATION_UI_DEPENDENCIES_MISSING: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(shiny))
source("src/interpretation_module.R")
future::plan(future::sequential)
# Local sources must not be read from a test: point every one at nothing.
Sys.setenv(REGSHINY_CADD_SNV = "/nonexistent/cadd", REGSHINY_CADD_INDEL = "/nonexistent/cadd", REGSHINY_VEP_CACHE = "/nonexistent/vep",
  REGSHINY_TOPMED_QTL = "/nonexistent/topmed", REGSHINY_FHS_SQTL = "/nonexistent/a;/nonexistent/b", REGSHINY_FHS_EQTL = "/nonexistent/e",
  REGSHINY_SPLICEAI = "/nonexistent/spliceai", REGSHINY_GWAS_CATALOG = "/nonexistent/gwas", REGSHINY_CLINVAR_VCF = "/nonexistent/clinvar",
  REGSHINY_REACTOME = "/nonexistent/reactome", REGSHINY_MSIGDB = "/nonexistent/msigdb", REGSHINY_TISSUE_MODEL_EQTL = "/nonexistent/model",
  REGSHINY_TISSUE_MODEL_SQTL = "/nonexistent/model", REGSHINY_PQTL_UKBPPP = "/nonexistent/pqtl", REGSHINY_PQTL_FHS = "/nonexistent/pqtl",
  REGSHINY_PQTL_PROTEIN_MAP = "/nonexistent/map", REGSHINY_GTEX_EQTL = "/nonexistent/gtex", REGSHINY_GTEX_EGENES = "/nonexistent/genes",
  REGSHINY_GENE_DISEASE = "/nonexistent/gene_disease")
scratch <- tempfile("interpretation-ui-"); dir.create(scratch)
backend_file <- normalizePath("src/interpretation_backend.R")
fixture_backend <- file.path(scratch, "backend.R")
annotation_source <- file.path(scratch, "annotation.R")
splice_source <- file.path(scratch, "splice.R")
request_log <- file.path(scratch, "requests.txt"); writeLines(character(), request_log)
writeLines(c(
  sprintf('sys.source(%s, envir = environment())', deparse(backend_file)),
  sprintf('request_log <- %s', deparse(request_log)),
  'interp_http_json <- function(url, body = NULL, timeout = 45) {',
  '  cat(url, if (is.null(body)) "" else body, "\\n", file = request_log, append = TRUE)',
  '  page <- function(data) list(data = data, paging_info = list(numberOfPages = 1L))',
  '  if (grepl("myvariant", url, fixed = TRUE)) return(list(dbsnp = list(rsid = "rs77"), clinvar = list(variant_id = 5L,',
  '    rcv = list(list(accession = "RCV1", clinical_significance = "Pathogenic", conditions = list(name = "Synthetic disease"), review_status = "criteria provided")))))',
  '  if (grepl("gtexportal", url, fixed = TRUE)) {',
  '    if (grepl("dataset/variant", url, fixed = TRUE)) return(page(list(list(variantId = "x"))))',
  '    if (grepl("singleTissueSqtl", url, fixed = TRUE)) return(page(list(list(tissueSiteDetailId = "Whole_Blood", geneSymbol = "PLUS",',
  '      gencodeId = "ENSGPLUS.1", phenotypeId = "chrT:161:200:clu_1:ENSGPLUS.1", pValue = 1e-15, nes = -1))))',
  '    if (grepl("singleTissueEqtl", url, fixed = TRUE)) return(page(list()))',
  '    if (grepl("reference/gene", url, fixed = TRUE)) return(page(list(list(gencodeId = "ENSG00000000001.5"))))',
  '    return(page(list(list(geneSymbol = "PLUS", tissueSiteDetailId = "Whole_Blood", median = 12), list(geneSymbol = "PLUS", tissueSiteDetailId = "Liver", median = 3))))',
  '  }',
  '  req <- jsonlite::fromJSON(body, simplifyVector = FALSE)',
  '  if (!is.null(req$variables$ids)) return(list(data = list(targets = list(list(id = "ENSG00000000001", approvedSymbol = "PLUS", approvedName = "plus gene",',
  '    biotype = "protein_coding", pathways = list(list(pathway = "Synthetic pathway", topLevelTerm = "Metabolism", pathwayId = "R-HSA-1")),',
  '    geneOntology = list(list(aspect = "P", evidence = "IDA", term = list(id = "GO:1", label = "synthetic process"))),',
  '    associatedDiseases = list(count = 1L, rows = list(list(score = .7, disease = list(id = "D1", name = "Synthetic disease"))))))))) ',
  '  if (is.null(req$variables$index)) return(list(data = list(variant = list(id = req$variables$id, rsIds = list("rs77"),',
  '    variantDescription = "Synthetic", mostSevereConsequence = list(label = "splice_donor_variant")))))',
  '  gwas <- grepl("studyTypes:[gwas]", req$query, fixed = TRUE)',
  '  rows <- if (gwas) list(list(studyType = "gwas", pValueMantissa = 2, pValueExponent = -20, beta = .3,',
  '    locus = list(rows = list(list(posteriorProbability = .8, is95CredibleSet = TRUE))),',
  '    study = list(id = "S1", traitFromSource = "Synthetic trait", diseases = list(list(name = "Synthetic disease"))))) else list()',
  '  list(data = list(variant = list(credibleSets = list(count = length(rows), rows = rows))))',
  '}',
  'interp_gwas_ld <- function(variant, vcf, ...) {',
  '  if (!nzchar(vcf)) stop("The dataset\'s VCF for this chromosome is not readable.")',
  '  ld <- interp_gwas_ld_empty()[1, ]; ld$rs_id <- "rs88"; ld$position <- variant$pos1 + 40L; ld$distance <- 40L; ld$ref <- "C"; ld$alt <- "T"',
  '  ld$r2 <- .93; ld$r <- .964; ld$n <- 100L; ld$risk_allele <- "rs88-T"; ld$risk_allele_travels_with <- "ALT"; ld$trait <- "Protein Z levels"',
  '  ld$mapped_trait <- "level of protein Z in blood"; ld$p_value <- 1e-20; ld$effect <- .2; ld$pubmed <- "1"; ld$study_accession <- "GCST1"',
  '  list(available = TRUE, found = TRUE, ld = ld, window = 500000L, keep_r2 = .2, positions_in_window = 3L, positions_examined = 3L, truncated = FALSE,',
  '    samples = 100L, release = vcf, ran_at = "now")',
  '}',
  '# Synthetic cis-pQTL scan and GTEx copy: the variant at 161 is in both, any other position in neither.',
  'interp_pqtl <- function(variant, sources = NULL) {',
  '  hit <- variant$pos1 == 161L',
  '  a <- interp_pqtl_empty()[1, ]; a$study <- "UKB-PPP"; a$protein <- "PLUS"; a$gene <- "PLUS"; a$gene_id <- "ENSG00000000001"; a$assay_id <- ""; a$rs_id <- "rs77"',
  '  a$beta_alt <- -.6; a$se <- .02; a$p_value <- 1e-30; a$log10_p <- -30; a$significance <- "experiment-wide"; a$alt_frequency <- .1; a$n <- 30000L; a$info <- 1; a$gene_distance <- 0L',
  '  list(available = TRUE, found = hit, associations = if (hit) a else a[0, ], studies = data.frame(study = "UKB-PPP", label = "synthetic scan", available = TRUE,',
  '    variant_tested = hit, proteins_tested = as.integer(hit), min_p = if (hit) 1e-30 else NA_real_, stringsAsFactors = FALSE), cut = interp_pqtl_cut(), source = "synthetic")',
  '}',
  'interp_gtex_local <- function(variant, sources = NULL) {',
  '  out <- c(list(available = TRUE, found = FALSE, release = "v99", gene_table = TRUE), interp_gtex_local_empty(), list(source = "synthetic GTEx"))',
  '  if (variant$pos1 != 161L) return(out)',
  '  out$pairs <- data.frame(gene = "PLUS", gene_id = "ENSG00000000001.1", tissue = c("Whole Blood", "Liver"), slope_alt = c(-.4, -.2), slope_se = .05, p_nominal = c(1e-12, 1e-5),',
  '    p_threshold = 1e-4, alt_frequency = .1, tss_distance = 60L, lead_variant = c(TRUE, FALSE), stringsAsFactors = FALSE)',
  '  out$genes <- data.frame(gene = "PLUS", gene_id = "ENSG00000000001.1", biotype = "protein_coding", tissues_significant = 2L, tissues_tested = 40L,',
  '    direction = "ALT lowers expression in every significant tissue", strongest_tissue = "Whole Blood", strongest_slope = -.4, strongest_p = 1e-12, lead_variant_in = 1L,',
  '    whole_blood = "significant", tissues = "Whole Blood, Liver", stringsAsFactors = FALSE)',
  '  out$found <- TRUE; out',
  '}',
  '# Two molecular traits at 161 (the protein above and a synthetic eQTL), compared by a canned colocalisation.',
  'interp_coloc_traits <- function(local, ...) {',
  '  hit <- !is.null(local$pqtl) && isTRUE(local$pqtl$found)',
  '  data.frame(kind = c("eQTL", "pQTL"), label = c("eQTL PLUS (Lung)", "pQTL PLUS (UKB-PPP)"), gene = "PLUS", source = c("Lung", "UKB-PPP"),',
  '    phenotype = c("ENSG00000000001.1", "PLUS"), p_value = c(1e-12, 1e-30), stringsAsFactors = FALSE)[if (hit) 1:2 else 0L, ]',
  '}',
  'interp_colocalise <- function(variant, local, ...) {',
  '  traits <- interp_coloc_traits(local); traits$variants <- 500L; traits$this_variant_posterior <- c(.7, .9)',
  '  pairs <- interp_coloc_empty()[1, ]; pairs$trait_1 <- traits$label[[1L]]; pairs$trait_2 <- traits$label[[2L]]; pairs$variants <- 480L; pairs$shared_cause <- .93',
  '  pairs$distinct_causes <- .05; pairs$only_trait_1 <- .01; pairs$only_trait_2 <- .01; pairs$neither <- 0; pairs$this_variant_share <- .81',
  '  pairs$leading_variant <- "this variant"; pairs$reading <- "one shared causal variant"',
  '  list(available = TRUE, traits = traits, pairs = pairs, window = 1000000L, priors = c(p1 = 1e-4), ran_at = "now")',
  '}',
  'interp_spliceai <- function(variant, fasta, ...) {',
  '  if (!nzchar(fasta)) stop("SpliceAI needs the reference FASTA of the connected dataset.")',
  '  list(scores = interp_parse_spliceai(c("#h", paste0("chrT\\t", variant$pos1, "\\t.\\tA\\tG\\t.\\t.\\tSpliceAI=G|PLUS|0.00|0.00|0.05|0.93|3|4|9|-1")), variant),',
  '    distance = 500L, model = "SYNTHETIC SpliceAI", ran_at = "now")',
  '}'), fixture_backend)
writeLines(c(
  'query_annotation <- function(resources, chrom, start1, end1, build) {',
  '  if (chrom == "chrBad") stop("Chromosome is not present in this annotation")',
  '  id <- "ENSG00000000001.1"',
  '  list(transcripts = data.frame(transcript_id = "TP", gene_id = id, gene_name = "PLUS", transcript_name = "TP", transcript_type = "protein_coding",',
  '      transcript_support_level = "1", tags = "basic;Ensembl_canonical;MANE_Select", chrom = chrom, start1 = 101L, end1 = 360L, strand = "+", source = "S", stringsAsFactors = FALSE),',
  '    exons = data.frame(transcript_id = "TP", gene_id = id, exon_id = paste0("E", 1:3), exon_number = as.character(1:3), chrom = chrom,',
  '      start1 = c(101L, 201L, 301L), end1 = c(160L, 260L, 360L), strand = "+", stringsAsFactors = FALSE),',
  '    features = data.frame(transcript_id = "TP", gene_id = id, feature = c("CDS", "CDS", "CDS", "stop_codon"), chrom = chrom,',
  '      start1 = c(121L, 201L, 301L, 318L), end1 = c(160L, 260L, 317L, 320L), strand = "+", frame = c("0", "2", "2", "0"), stringsAsFactors = FALSE),',
  '    introns = data.frame(transcript_id = "TP", strand = "+", donor1 = c(161L, 261L), acceptor1 = c(200L, 300L), stringsAsFactors = FALSE),',
  '    metadata = list(release = resources$release))',
  '}',
  'annotation_gene_names <- function(resources, gene_ids) data.frame(gene_id = character(), gene_name = character(), gene_type = character())'), annotation_source)
writeLines(c(
  'splice_fetch_reference <- function(reference, chrom, start1, end1) {',
  '  genome <- rep("A", 400L); genome[121:129] <- strsplit("ATGGCTTGG", "")[[1L]]; genome[318:320] <- c("T", "A", "A")',
  '  list(sequence = paste(genome[start1:end1], collapse = ""))',
  '}'), splice_source)

# Static pieces.
ui <- as.character(interpretation_ui("x"))
for (id in c("x-public", "x-spliceai", "x-summary", "x-consequences", "x-context_plot", "x-gwas_slot", "x-pathways_slot", "x-local_sources", "x-pqtl_table",
  "x-pqtl_status", "x-gtex_local_plot_box", "x-gtex_local_genes", "x-gtex_local_pairs", "x-splice_protein", "x-splice_protein_status", "x-coloc_table",
  "x-coloc_status", "x-gene_disease", "x-gene_disease_status"))
  stopifnot(grepl(paste0('id="', id, '"'), ui, fixed = TRUE))
# No deferred DT containers in the initial UI; local tables are still present.
deferred <- interpretation_deferred_tables()
stopifnot(length(unlist(deferred)) == 12L)
for (id in unlist(deferred)) stopifnot(!grepl(paste0('id="x-', id, '"'), ui, fixed = TRUE),
  grepl(paste0('id="x-', id, '_slot"'), ui, fixed = TRUE))
for (id in c("clinvar", "catalog", "gtex_local_pairs", "linked_genes"))
  stopifnot(grepl(paste0('id="x-', id, '"'), ui, fixed = TRUE))
stopifnot(grepl("sends only the variant", ui, fixed = TRUE), grepl("not a clinical classification", ui, fixed = TRUE))
cards <- as.character(interpretation_summary_cards(data.frame(question = c("Splicing", "Tissue"), level = c("Strong", "Not assessed"),
  finding = c("f1", "f2"), evidence = c("a | b", ""), sources = c("S", ""), stringsAsFactors = FALSE)))
stopifnot(grepl("regshiny-level-strong", cards, fixed = TRUE), grepl("regshiny-level-unknown", cards, fixed = TRUE),
  lengths(regmatches(cards, gregexpr("<li>", cards, fixed = TRUE))) == 2L, grepl("Sources: S", cards, fixed = TRUE))
stopifnot(interpretation_level_class("Moderate") == "regshiny-level-moderate", interpretation_level_class("Context") == "regshiny-level-context")
# A long evidence list shows its first four lines; the rest open on request and nothing is dropped.
long_card <- as.character(interpretation_summary_cards(data.frame(question = "Splicing", level = "Strong", finding = "f", evidence = paste(paste0("line", 1:6), collapse = " | "),
  sources = "S", stringsAsFactors = FALSE)))
stopifnot(lengths(regmatches(long_card, gregexpr("<li>", long_card, fixed = TRUE))) == 6L, grepl("<summary>2 more evidence lines</summary>", long_card, fixed = TRUE),
  regexpr("line4", long_card, fixed = TRUE) < regexpr("<details", long_card, fixed = TRUE), regexpr("line5", long_card, fixed = TRUE) > regexpr("<details", long_card, fixed = TRUE),
  !grepl("<details", cards, fixed = TRUE), grepl(".regshiny-evidence-more", css <- paste(readLines("www/frontend.css", warn = FALSE), collapse = "\n"), fixed = TRUE))
css <- paste(readLines("www/frontend.css", warn = FALSE), collapse = "\n")
for (cls in c("regshiny-level-strong", "regshiny-level-moderate", "regshiny-level-limited", "regshiny-level-none", "regshiny-level-context",
  "regshiny-level-unknown", "regshiny-evidence-row", "regshiny-cell-wrap", "regshiny-site-list")) stopifnot(grepl(paste0(".", cls), css, fixed = TRUE))
# Table helper: p-value columns keep their stored value and get the scientific renderer, long text gets the wrapping
# renderer, other numbers are rounded; rows stay on one line.
stopifnot(identical(interpretation_p_columns(c("p_value", "p_nominal", "p_threshold", "min_p", "q_p_value", "p_permutation", "pip",
  "posterior_probability", "log10_p", "phenotype", "pooled_slope")), c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)))
tab <- interpretation_datatable(data.frame(trait = "x", p_value = 0, pip = .123456, min_p = 1.23456e-9, log10_p = -3.21987, n = 5L, stringsAsFactors = FALSE))
rendered <- Filter(function(d) !is.null(d$render), tab$x$options$columnDefs)
stopifnot(inherits(tab, "datatables"), length(rendered) == 2L, identical(as.numeric(rendered[[1L]]$targets), c(1, 3)), grepl("1e-300", rendered[[1L]]$render, fixed = TRUE),
  identical(as.numeric(rendered[[2L]]$targets), 0), grepl("regshiny-cell-wrap", rendered[[2L]]$render, fixed = TRUE),
  tab$x$data$p_value == 0, tab$x$data$min_p == 1.23456e-9, tab$x$data$pip == .123, tab$x$data$log10_p == -3.22, tab$x$data$n == 5L,
  grepl("nowrap", tab$x$container, fixed = TRUE) || grepl("nowrap", paste(unlist(tab$x), collapse = " "), fixed = TRUE),
  # A short table needs no search box or pager.
  tab$x$options$dom == "t", interpretation_datatable(data.frame(n = 1:11))$x$options$dom == "ftip")
# Plots draw with and without data.
pair <- function(type, context, phenotype, slope, significant) data.frame(qtl_type = type, context = context, phenotype = phenotype,
  gene_id = "ENSG1.1", gene = "ONE", alt_frequency = .1, n_estimated = 500L, slope_alt = slope, slope_se = .1, p_nominal = 1e-9,
  p_threshold = 1e-5, significant = significant, pip = NA_real_, credible_set = "", tss_distance = 1L, stringsAsFactors = FALSE)
pairs <- rbind(pair("sQTL", "Lung", "chr1:10:90:clu_5_+:ENSG1.1", 2, TRUE), pair("sQTL", "PBMC", "chr1:10:90:clu_9_+:ENSG1.1", 1.8, FALSE),
  pair("eQTL", "Lung", "ENSG1.1", .4, TRUE))
e <- new.env(); sys.source(backend_file, envir = e)
png <- file.path(scratch, "contexts.png"); grDevices::png(png, width = 1400, height = 600, res = 130)
shown <- interpretation_plot_contexts(pairs, e$interp_context_heterogeneity(pairs)); grDevices::dev.off()
# With model estimates the same panels are drawn, the estimates as diamonds.
grDevices::png(png, width = 1400, height = 600, res = 130)
with_model <- interpretation_plot_contexts(pairs, e$interp_context_heterogeneity(pairs), model = data.frame(qtl_type = "sQTL", phenotype = "chr1:10:90:ENSG1.1",
  context = c("Lung", "PBMC"), posterior_slope = c(1.9, 1.85), call = c("effect", "no clear effect"), stringsAsFactors = FALSE)); grDevices::dev.off()
stopifnot(nrow(with_model) == 2L)
stopifnot(nrow(shown) == 2L, file.info(png)$size > 1000)
grDevices::png(png, width = 900, height = 500); stopifnot(is.null(interpretation_plot_contexts(pairs[0, ], NULL))); grDevices::dev.off()
grDevices::png(png, width = 900, height = 700)
bars <- interpretation_plot_expression(data.frame(gene = "G", tissue = c("Whole Blood", "Liver"), median_tpm = c(5, 50)), "G")
stopifnot(identical(bars$tissue, c("Whole Blood", "Liver")), is.null(interpretation_plot_expression(data.frame(gene = character(), tissue = character(), median_tpm = numeric()), "G")))
grDevices::dev.off()
gtex_pairs <- data.frame(gene = "G", gene_id = "ENSG1.1", tissue = c("Whole Blood", "Liver", "Lung"), slope_alt = c(.3, -.2, .1), slope_se = .05,
  lead_variant = c(TRUE, FALSE, NA), stringsAsFactors = FALSE)
grDevices::png(png, width = 900, height = 700)
drawn <- interpretation_plot_gtex(gtex_pairs, "ENSG1.1", "v99")
stopifnot(identical(drawn$tissue, c("Liver", "Lung", "Whole Blood")), is.null(interpretation_plot_gtex(gtex_pairs, "ENSG2.1")), is.null(interpretation_plot_gtex(NULL, "")))
grDevices::dev.off()

await <- function(session, task, ready) {
  for (i in seq_len(200L)) {
    later::run_now(.05); session$flushReact()
    if (shiny::isolate(task$status()) == "error") stop(shiny::isolate(task$result()))
    if (shiny::isolate(task$status()) == "success" && isTRUE(tryCatch(shiny::isolate(ready()), error = function(e) FALSE))) {
      later::run_now(.01); session$flushReact(); return(invisible(TRUE))
    }
  }
  stop("Interpretation task timed out")
}
html_text <- function(x) gsub("[[:space:]]+", " ", if (is.list(x) && !is.null(x$html)) as.character(x$html) else as.character(x))
level <- function(s, question) s$level[s$question == question]
record <- function(chrom = "chrT", pos = 161L, ref = "A", alt = "G", build = "GRCh38")
  data.frame(chrom = chrom, pos1 = pos, ref = ref, alt = alt, build = build, vcf = "/synthetic/chrT.vcf.gz", stringsAsFactors = FALSE)
current <- reactiveVal(record())
bins <- reactiveVal(data.frame(site_id = "donor:+:chrT:161", splice_site_type = "donor", splice_site_position = 161L, status = "OK", beta = -.4,
  p_permutation = 1e-4, n = 800L, bin_check = "OK", stringsAsFactors = FALSE))
release <- reactiveVal("48")
shiny::testServer(interpretation_server, args = list(variant = reactive(current()),
  annotation_resources = function() list(release = release()), annotation_file = annotation_source, reference = list(fasta = "/synthetic.fa"),
  splice_file = splice_source, interpretation_file = fixture_backend, bin_ratio = reactive(bins()), junctions = NULL), {
  session$flushReact(); await(session, local_task, function() !is.null(local_now()))
  check_slots <- function(ids, ready) {
    for (id in ids) {
      html <- html_text(output[[paste0(id, "_slot")]])
      if (ready) stopifnot(grepl("datatables", html, fixed = TRUE), grepl(paste0('id="', session$ns(id), '"'), html, fixed = TRUE)) else
        stopifnot(!grepl("datatables", html, fixed = TRUE), grepl("Use .* to show this table", html))
    }
  }
  check_slots(unlist(deferred), FALSE)
  x <- isolate(local_now())
  stopifnot(x$variant$id == "chrT:161:A:G", x$consequences$consequence == "splice_donor_variant", x$consequences$hgvs_c == "c.40+1A>G",
    x$release == "48", x$reference_base == "A", is.null(x$annotation_error), !any(x$local$status$available))
  s <- isolate(summary_now())
  # Essential site plus the cohort's own bin-ratio association: Strong without any external source.
  stopifnot(nrow(s) == 7L, level(s, "Splicing") == "Strong", level(s, "Protein sequence") == "Strong", level(s, "Disease and traits") == "Not assessed",
    level(s, "Pathway") == "Not assessed", grepl("bin ratio at donor:+:chrT:161", s$evidence[s$question == "Splicing"], fixed = TRUE),
    # The local cis-pQTL scan and GTEx copy are read with the other local sources: a cis-pQTL of the gene whose splice site
    # the variant destroys is Strong, and the gene joins the linked genes through both.
    level(s, "Protein abundance (pQTL)") == "Strong", grepl("PLUS in UKB-PPP: -0.600 per ALT allele", s$evidence[s$question == "Protein abundance (pQTL)"], fixed = TRUE),
    level(s, "Gene expression (eQTL)") == "Moderate", grepl("GTEx v99 eQTL (local): PLUS in 2 of 40 tested tissue(s)", s$evidence[s$question == "Gene expression (eQTL)"], fixed = TRUE),
    level(s, "Tissue") == "Context", grepl("the gene's lead variant in 1 tissue(s)", s$evidence[s$question == "Tissue"], fixed = TRUE),
    x$local$pathways$genes$relation == "overlaps; eQTL target (GTEx); cis-pQTL protein")
  stopifnot(grepl("this variant was tested against 1 nearby protein(s); smallest p 1.0e-30", html_text(output$pqtl_status), fixed = TRUE),
    grepl("significant cis-eQTL of 1 gene(s) in 2 tissue(s)", html_text(output$gtex_local_status), fixed = TRUE), isolate(gtex_gene()) == "ENSG00000000001.1",
    is.null(output$gtex_gene_choice) || !nzchar(html_text(output$gtex_gene_choice)),
    # Two tissues: the plot box is sized for the three-row minimum.
    grepl("height:222px", html_text(output$gtex_local_plot_box), fixed = TRUE))
  header <- html_text(output$header)
  stopifnot(grepl("chrT:161:A:G", header, fixed = TRUE), grepl("GENCODE v48", header, fixed = TRUE), !grepl("does not match", header, fixed = TRUE))
  # The splice outcome is computed on its own once the local sources are in: an essential donor with nothing
  # to show what replaces it gives the assumed outcome (first exon: retention only), read with the ALT allele.
  await(session, protein_task, function() !is.null(protein_now()))
  outcome <- isolate(splice_protein_table())
  stopifnot(nrow(outcome) == 1L, outcome$transcript_id == "TP", outcome$basis == "assumed outcome", outcome$event == "if the intron is retained (40 nt)",
    # The open result's bin ratio at that donor falls with ALT: the cohort's own line joins the variant's.
    outcome$support == "the variant changes the essential donor site; this cohort: donor site at 161 used less with ALT (bin ratio -0.400 per ALT allele, permutation p 1.0e-04)",
    outcome$protein_length_ref == 39L, nzchar(outcome$coding_effect),
    grepl("1 splice site(s) of the open result whose bin ratio moves with ALT, 0 SpliceAI site(s) (SpliceAI not run)", html_text(output$splice_protein_status), fixed = TRUE),
    grepl("1 outcome(s) on coding transcripts", html_text(output$splice_protein_status), fixed = TRUE),
    grepl("through splicing (assumed outcome), TP: if the intron is retained (40 nt); 40 intronic nt kept", isolate(summary_now())$evidence[[1L]], fixed = TRUE),
    level(isolate(summary_now()), "Protein sequence") == "Strong")
  # The variant has two molecular traits: they are compared by themselves, and the shared cause is named under both.
  await(session, coloc_task, function() !is.null(coloc_now()))
  s <- isolate(summary_now())
  stopifnot(nrow(isolate(coloc_now())$pairs) == 1L,
    grepl("2 trait(s), with this variant's posterior probability within each trait's own signal: eQTL PLUS (Lung) 0.70; pQTL PLUS (UKB-PPP) 0.90", html_text(output$coloc_status), fixed = TRUE),
    grepl("1 with one shared causal variant, 0 with different ones", html_text(output$coloc_status), fixed = TRUE),
    grepl("colocalisation, eQTL PLUS (Lung) and pQTL PLUS (UKB-PPP): one shared causal variant (posterior 0.93; this variant's share of the shared signal 0.81)",
      s$evidence[s$question == "Protein abundance (pQTL)"], fixed = TRUE),
    grepl("colocalisation, eQTL PLUS (Lung)", s$evidence[s$question == "Gene expression (eQTL)"], fixed = TRUE),
    # The reading box ends with what fine-mapping and colocalisation add; the assumed splice outcome is left out of it.
    grepl("What fine-mapping, the splice prediction and colocalisation add", html_text(output$reading), fixed = TRUE),
    # 0.90 within the protein's signal counts as leading it; 0.70 within the eQTL does not.
    grepl("likely cause of the plasma level of PLUS (posterior 0.90).", html_text(output$reading), fixed = TRUE),
    !grepl("the expression of PLUS (posterior 0.70)", html_text(output$reading), fixed = TRUE),
    grepl("One causal variant is shared by the expression of PLUS and the plasma level of PLUS (posterior 0.93).", html_text(output$reading), fixed = TRUE),
    !grepl("The splice change it is tied to", html_text(output$reading), fixed = TRUE), length(isolate(mechanism_now())) == 2L)
  stopifnot(grepl("regshiny-level-strong", html_text(output$summary), fixed = TRUE), !length(readLines(request_log)))
  # The reading names the best-supported effect and the steps that have not been run.
  stopifnot(grepl("Best-supported molecular effect: a change to the protein sequence and protein abundance and splicing (Strong)", html_text(output$reading), fixed = TRUE),
    grepl("Not run yet: SpliceAI", html_text(output$reading), fixed = TRUE))
  stopifnot(grepl("has not been run", html_text(output$spliceai_status), fixed = TRUE))
  # Nothing is sent until the user asks.
  session$setInputs(public = 1L); await(session, public_task, function() !is.null(public_now()))
  check_slots(deferred$public, TRUE); check_slots(c(deferred$spliceai, deferred$ld), FALSE)
  p <- isolate(public_now()); s <- isolate(summary_now()); sent <- paste(readLines(request_log), collapse = "\n")
  stopifnot(p$open_targets$found, nrow(p$open_targets$gwas) == 1L, level(s, "Disease and traits") == "Strong", level(s, "Pathway") == "Context",
    level(s, "Tissue") == "Context", level(s, "Protein abundance (pQTL)") == "Strong", grepl("GTEx sQTL in 1 tissue", s$evidence[s$question == "Splicing"], fixed = TRUE),
    grepl("T_161_A_G", sent, fixed = TRUE), grepl("ENSG00000000001", sent, fixed = TRUE),
    # Only the variant and gene identifiers leave: no site table, sample count or cohort label.
    !grepl("bin ratio|donor:|800|FHS", sent))
  stopifnot(grepl("Public sources queried", html_text(output$public_status), fixed = TRUE))
  session$setInputs(spliceai = 1L); await(session, spliceai_task, function() !is.null(spliceai_now()))
  check_slots(deferred$spliceai, TRUE)
  sites <- isolate(spliceai_sites())
  stopifnot(nrow(sites) == 1L, sites$effect == "donor loss", sites$site_position == 161L, sites$annotated_site == "yes",
    grepl("bin ratio -0.400", sites$observed_in_cohort, fixed = TRUE), grepl("SYNTHETIC SpliceAI", html_text(output$spliceai_status), fixed = TRUE))
  stopifnot(grepl("donor loss 0.93 at exon base 160", isolate(summary_now())$evidence[3L], fixed = TRUE))
  # The SpliceAI site joins the splice changes and the outcome is computed again with it.
  await(session, protein_task, function() !is.null(protein_now()) && nrow(protein_now()$events) == 2L)
  stopifnot(identical(isolate(protein_now())$events$source, c("this cohort", "SpliceAI")), all(isolate(protein_now())$events$effect == "loss"), nrow(isolate(splice_protein_table())) == 1L,
    grepl("1 SpliceAI site(s);", html_text(output$splice_protein_status), fixed = TRUE))
  # The LD scan runs on request with the record's own VCF and feeds Disease and the pQTL question.
  stopifnot(is.null(isolate(ld_now())), grepl("Not computed for this variant yet", html_text(output$ld_status), fixed = TRUE))
  session$setInputs(ld = 1L); await(session, ld_task, function() !is.null(ld_now()))
  check_slots(deferred$ld, TRUE)
  l <- isolate(ld_now()); s <- isolate(summary_now())
  stopifnot(l$release == "/synthetic/chrT.vcf.gz", grepl("3 GWAS Catalog position(s) within 500 kb", html_text(output$ld_status), fixed = TRUE),
    grepl("1 variant(s) with r2 of at least 0.2, 1 of them at 0.8 or more", html_text(output$ld_status), fixed = TRUE),
    grepl("r2 >= 0.8 with 1 other reported variant(s)", s$evidence[s$question == "Disease and traits"], fixed = TRUE),
    grepl("through LD (r2 >= 0.8 in this cohort): protein Z (rs88, r2 0.93)", s$evidence[s$question == "Protein abundance (pQTL)"], fixed = TRUE))
  # A new variant starts clean: earlier public and SpliceAI answers are not shown for it.
  current(record(pos = 124L, ref = "G", alt = "A")); session$flushReact(); await(session, local_task, function() identical(local_now()$variant$pos1, 124L))
  stopifnot(is.null(isolate(public_now())), is.null(isolate(spliceai_now())), is.null(isolate(ld_now())))
  check_slots(unlist(deferred), FALSE)
  s <- isolate(summary_now())
  stopifnot(isolate(local_now())$consequences$hgvs_p == "p.Ala2Thr", level(s, "Protein sequence") == "Moderate", level(s, "Disease and traits") == "Not assessed",
    # This position is in neither the scan nor GTEx.
    level(s, "Protein abundance (pQTL)") == "Not assessed", level(s, "Gene expression (eQTL)") == "None found",
    grepl("this variant is not in the scan", html_text(output$pqtl_status), fixed = TRUE), grepl("not a significant cis-eQTL of any gene", html_text(output$gtex_local_status), fixed = TRUE),
    # With fewer than two traits nothing is compared.
    is.null(isolate(coloc_now())), grepl("Not applicable: the variant is associated with 0 molecular trait(s)", html_text(output$coloc_status), fixed = TRUE),
    !length(isolate(mechanism_now())), !grepl("What fine-mapping", html_text(output$reading), fixed = TRUE))
  # Going back shows the answers already fetched for that variant.
  current(record()); session$flushReact(); await(session, local_task, function() identical(local_now()$variant$pos1, 161L))
  stopifnot(!is.null(isolate(public_now())), !is.null(isolate(spliceai_now())), !is.null(isolate(ld_now())))
  check_slots(unlist(deferred), TRUE)
  # Changing the GENCODE release recomputes consequences and keeps the public answer.
  release("50"); session$flushReact(); await(session, local_task, function() identical(local_now()$release, "50"))
  stopifnot(!is.null(isolate(public_now())))
  # Multi-allelic record: the first ALT by default, the chosen one after selection.
  current(record(pos = 124L, ref = "G", alt = "A,T")); session$flushReact(); await(session, local_task, function() identical(local_now()$variant$alt, "A"))
  stopifnot(grepl("ALT allele to interpret", html_text(output$header), fixed = TRUE))
  session$setInputs(alt = "T"); await(session, local_task, function() identical(local_now()$variant$alt, "T"))
  stopifnot(isolate(local_now())$consequences$hgvs_p == "p.Ala2Ser")
  # REF that disagrees with the reference is flagged at the top.
  current(record(pos = 124L, ref = "C", alt = "A")); session$flushReact(); await(session, local_task, function() identical(local_now()$variant$ref, "C"))
  stopifnot(grepl("does not match the reference sequence (G)", html_text(output$header), fixed = TRUE))
  # Annotation failure and a non-GRCh38 assembly degrade without a crash.
  current(record(chrom = "chrBad")); session$flushReact(); await(session, local_task, function() identical(local_now()$variant$chrom, "chrBad"))
  stopifnot(grepl("not present in this annotation", isolate(local_now())$annotation_error, fixed = TRUE), nrow(isolate(summary_now())) == 7L,
    grepl("Transcript consequences unavailable", html_text(output$header), fixed = TRUE))
  current(record(build = "GRCh37")); session$flushReact(); await(session, local_task, function() is.null(local_now()$local) && identical(local_now()$variant$chrom, "chrT"))
  stopifnot(grepl("GRCh38 only", html_text(output$header), fixed = TRUE), nrow(isolate(summary_now())) == 7L)
})
# A junction of the open result that rises with ALT is re-spliced into the transcript: exon 2 left out, in frame.
rising <- reactiveVal(data.frame(chrom = "chrT", intron_start1 = 161L, intron_end1 = 300L, status = "OK", beta = 3, p_permutation = 1e-4, n = 800L, stringsAsFactors = FALSE))
shiny::testServer(interpretation_server, args = list(variant = reactive(record(pos = 175L, ref = "A", alt = "G")), annotation_resources = list(release = "48"),
  annotation_file = annotation_source, reference = list(fasta = "/synthetic.fa"), splice_file = splice_source, interpretation_file = fixture_backend,
  junctions = reactive(rising())), {
  session$flushReact(); await(session, local_task, function() !is.null(local_now()))
  await(session, protein_task, function() !is.null(protein_now()))
  outcome <- isolate(splice_protein_table()); s <- isolate(summary_now())
  stopifnot(nrow(outcome) == 1L, outcome$event == "exon skipping: 1 exon(s) left out (60 nt)", outcome$basis == "observed", outcome$coding_effect == "in-frame deletion",
    outcome$protein_length_new == 19L, grepl("this cohort: junction used more with ALT (+3.00 reads per 100", outcome$support, fixed = TRUE),
    level(s, "Protein sequence") == "Moderate", grepl("The variant itself does not change the protein sequence. Through the splice change (observed): residues ", s$finding[[1L]], fixed = TRUE),
    grepl("(20 residues) are removed; the rest of the protein is unchanged (19 instead of 39 residues).", s$finding[[1L]], fixed = TRUE),
    grepl("1 intron(s) used more and 0 used less with ALT (0 from TOPMed sQTL, 1 from the open result)", html_text(output$splice_protein_status), fixed = TRUE))
})
# Without a reference sequence the outcome is not computed and says why.
shiny::testServer(interpretation_server, args = list(variant = reactive(record()), annotation_resources = list(release = "48"), annotation_file = annotation_source,
  reference = simpleError("no FASTA"), splice_file = splice_source, interpretation_file = fixture_backend), {
  session$flushReact(); await(session, local_task, function() !is.null(local_now()))
  await(session, protein_task, function() !is.null(protein_now()))
  stopifnot(inherits(isolate(protein_now()), "error"), is.null(isolate(splice_protein_table())),
    grepl("Not computed: the reference sequence is unavailable.", html_text(output$splice_protein_status), fixed = TRUE), nrow(isolate(summary_now())) == 7L)
})
# No reference: consequences still classify by position, SpliceAI reports why it cannot run.
shiny::testServer(interpretation_server, args = list(variant = reactive(record(pos = 124L, ref = "G", alt = "A")),
  annotation_resources = list(release = "48"), annotation_file = annotation_source, reference = simpleError("no FASTA"),
  splice_file = splice_source, interpretation_file = fixture_backend), {
  session$flushReact(); await(session, local_task, function() !is.null(local_now()))
  x <- isolate(local_now())
  stopifnot(x$consequences$consequence == "coding_sequence_variant", grepl("could not be read", x$consequences$notes), is.na(x$reference_base),
    !isolate(rna_now())$available)
  session$setInputs(spliceai = 1L); await(session, spliceai_task, function() !is.null(spliceai_now()))
  stopifnot(inherits(isolate(spliceai_now()), "error"), grepl("reference FASTA", html_text(output$spliceai_status), fixed = TRUE))
})
# A click made before the local lookup has come back is kept and runs once it has.
writeLines(character(), request_log)
shiny::testServer(interpretation_server, args = list(variant = reactive(record()), annotation_resources = list(release = "48"),
  annotation_file = annotation_source, reference = list(fasta = "/synthetic.fa"), splice_file = splice_source, interpretation_file = fixture_backend), {
  session$flushReact()
  session$setInputs(spliceai = 1L, public = 1L)
  stopifnot(is.null(isolate(local_ready())), isolate(pending$spliceai) == isolate(variant_key()), !length(readLines(request_log)),
    grepl("will start as soon as they are in", html_text(output$public_status), fixed = TRUE))
  await(session, spliceai_task, function() !is.null(spliceai_now())); await(session, public_task, function() !is.null(public_now()))
  stopifnot(!inherits(isolate(spliceai_now()), "error"), isolate(public_now())$open_targets$found, is.null(isolate(pending$spliceai)),
    is.null(isolate(pending$public)), grepl("T_161_A_G", paste(readLines(request_log), collapse = "\n"), fixed = TRUE))
})
# A click whose variant was replaced before the lookup came back is dropped: nothing is sent for either variant.
writeLines(character(), request_log)
moving <- reactiveVal(record())
shiny::testServer(interpretation_server, args = list(variant = reactive(moving()), annotation_resources = list(release = "48"),
  annotation_file = annotation_source, reference = list(fasta = "/synthetic.fa"), splice_file = splice_source, interpretation_file = fixture_backend), {
  session$flushReact()
  session$setInputs(public = 1L)
  moving(record(pos = 124L, ref = "G", alt = "A")); session$flushReact()
  await(session, local_task, function() identical(local_now()$variant$pos1, 124L))
  for (i in 1:5) { later::run_now(.05); session$flushReact() }
  stopifnot(is.null(isolate(pending$public)), is.null(isolate(public_now())), !length(readLines(request_log)))
})
unlink(scratch, recursive = TRUE)
cat("interpretation UI tests passed\n")
