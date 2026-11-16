# Public-copy links use public identifiers only; no lookup runs while rendering.
source("src/interpretation_backend.R")
source("src/interpretation_module.R")
ids <- list(chrom = "22", pos = 19963748L, ref = "G", alt = "A", rs_id = "rs4680", gene_id = "ENSG00000093010", gene = "COMT", clinvar_id = "17591")
expected <- c(open_targets = "https://platform.opentargets.org/variant/22_19963748_G_A",
  open_targets_gene = "https://platform.opentargets.org/target/ENSG00000093010",
  gtex = "https://gtexportal.org/home/snp/rs4680", gtex_gene = "https://gtexportal.org/home/gene/COMT",
  gtex_local = "https://gtexportal.org/home/snp/rs4680", gwas_catalog = "https://www.ebi.ac.uk/gwas/variants/rs4680",
  reactome_local = "https://reactome.org/content/query?q=COMT",
  clinvar_local = "https://www.ncbi.nlm.nih.gov/clinvar/variation/VCV000017591",
  clinvar_api = "https://www.ncbi.nlm.nih.gov/clinvar/variation/VCV000017591",
  cadd = "https://cadd.gs.washington.edu/snv/GRCh38-v1.7/22%3A19963748_G_A")
stopifnot(identical(names(expected), interp_source_links()$key))
for (key in names(expected)) {
  stopifnot(identical(interp_source_url(key, ids), unname(expected[[key]])),
    identical(interp_source_url(key, c(ids, list(sample = "private", genotype = "private", cohort = "private"))), unname(expected[[key]])))
  item <- data.frame(source = key, role = "evidence", gene_id = "", gene = "", clinvar_id = "")
  html <- as.character(interpretation_evidence_link("a source", item, ids))
  stopifnot(grepl(htmltools::htmlEscape(expected[[key]], attribute = TRUE), html, fixed = TRUE),
    grepl('target="_blank"', html, fixed = TRUE), grepl('rel="noopener noreferrer"', html, fixed = TRUE))
  item$role <- "notice"
  stopifnot(!grepl("href=", as.character(interpretation_evidence_link("not consulted", item, ids)), fixed = TRUE))
}
stopifnot(grepl("Release V10", interpretation_source_note("gtex_local"), fixed = TRUE),
  grepl("may differ", interpretation_source_note("gtex_local"), fixed = TRUE),
  interpretation_source_info("msigdb")$kind == "local copy", interpretation_source_info("myvariant")$kind == "via API")
for (key in c("topmed", "fhs", "cohort_ld", "cohort_rna", "colocalisation", "splice_outcome", "pqtl_finemap", "msigdb", "myvariant"))
  stopifnot(is.null(interp_source_url(key, ids)))
stopifnot(is.null(interp_source_url("clinvar_local", modifyList(ids, list(clinvar_id = "17591?sample=x")))),
  is.null(interp_source_url("cadd", modifyList(ids, list(alt = "AT")))),
  identical(interp_source_url("clinvar_local", modifyList(ids, list(clinvar_id = "VCV000017591"))), expected[["clinvar_local"]]))

# Exercise the real local VCF parser on a synthetic row (no native tool or dataset).
e <- new.env(parent = globalenv()); sys.source("src/interpretation_backend.R", envir = e)
f <- tempfile("synthetic-clinvar-"); writeLines("22\t19963748\t17591\tG\tA\t.\t.\tCLNSIG=Benign;GENEINFO=COMT:1312;RS=4680", f)
file.create(paste0(f, ".tbi"))
e$interp_run_lines <- function(...) "##fileDate=20990101"
e$interp_tabix <- function(...) readLines(f)
v <- e$interp_variant("chr22", 19963748L, "G", "A")
cv <- e$interp_clinvar_local(v, list(tabix = unname(Sys.which("true")), clinvar_vcf = f))
stopifnot(cv$found, identical(cv$records$variation_id, "17591"))
local_summary <- e$interp_evidence_summary(e$interp_empty_consequences(), local = list(clinvar = cv), structured = TRUE)
x <- local_summary$items[[5L]]; row <- x[x$source == "clinvar_local", ]
stopifnot(nrow(row) == 1L, row$clinvar_id == "17591", grepl("VCV000017591", as.character(interpretation_evidence_link(row$text, row, list())), fixed = TRUE))
myv <- e$interp_parse_myvariant(list(dbsnp = list(rsid = "rs4680"), clinvar = list(variant_id = 17591, rcv = list(accession = "RCV1", clinical_significance = "Benign"))))
pub_summary <- e$interp_evidence_summary(e$interp_empty_consequences(), public = list(myvariant = myv), structured = TRUE)
x <- pub_summary$items[[5L]]; row <- x[x$source == "clinvar_api", ]
stopifnot(nrow(row) == 1L, row$clinvar_id == "17591", grepl("VCV000017591", as.character(interpretation_evidence_link(row$text, row, list())), fixed = TRUE))
unlink(c(f, paste0(f, ".tbi")))
cat("SUMMARY_SOURCE_LINKS_OK: every template; notice exclusion; local VCF and API ClinVar identity\n")

source("tests/fixtures/summary_fixture.R")
dedup <- summary_source_dedup_fixture(); original_summary <- dedup$summary
x <- dedup$items[[6L]]
displays <- lapply(seq_len(nrow(x)), function(i) interpretation_source_display(x[i, ], dedup$identifiers))
stopifnot(identical(vapply(displays, `[[`, character(1), "label"), c("Reactome (COMT)", "Reactome (TXNRD2)", "MSigDB", "MSigDB")),
  displays[[1L]]$identity != displays[[2L]]$identity, displays[[3L]]$identity == displays[[4L]]$identity)
# Duplicate linked entries also collapse when their displayed text and URL match.
badges <- as.character(interpretation_source_badges(rbind(x, x[1L, ]), dedup$identifiers))
count <- function(text, pattern) sum(gregexpr(pattern, text, fixed = TRUE)[[1L]] > 0L)
stopifnot(count(badges, 'class="regshiny-source-badge"') == 3L, count(badges, 'Reactome (') == 2L, count(badges, '>MSigDB<') == 1L,
  grepl('href="https://reactome.org/content/query?q=COMT"', badges, fixed = TRUE),
  grepl('href="https://reactome.org/content/query?q=TXNRD2"', badges, fixed = TRUE))
# Evidence entries retain their original text and individual destinations.
for (i in 1:2) stopifnot(grepl(paste0('href="https://reactome.org/content/query?q=', x$gene[[i]], '"'),
  as.character(interpretation_evidence_link(x$text[[i]], x[i, ], list())), fixed = TRUE))
graph <- as.character(interpretation_summary_map(dedup$summary, dedup$items, interp_request_states(list(), list(), list())))
stopifnot(count(graph, 'class="regshiny-map-source"') == 3L, count(graph, 'class="regshiny-map-edge contributed"') == 4L,
  count(graph, 'data-contributions="2"') == 2L, count(graph, 'data-contributions="1"') == 2L,
  count(graph, 'data-question="Pathway"') == 3L, count(graph, 'data-question="Tissue"') == 1L,
  identical(dedup$summary, original_summary))
cat("SUMMARY_SOURCE_DEDUP_OK: Table Reactome=2 distinct labels, MSigDB=1; Map sources=3, question-source edges=4, contributions=6 (none lost)\n")
