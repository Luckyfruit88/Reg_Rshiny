# Pure presentation metadata shared by the backend and UI. No network or file data.
interp_request_steps <- function() data.frame(
  key = c("public", "spliceai", "ld"),
  label = c("Query public databases", "Run SpliceAI", "LD with GWAS variants"),
  duration = c("", " (about 1-2 min)", " (about 30 s)"), stringsAsFactors = FALSE)

interp_summary_scope <- function() list(
  `Protein sequence` = c("public", "spliceai"),
  `Protein abundance (pQTL)` = c("public", "spliceai", "ld"),
  Splicing = c("public", "spliceai"), `Gene expression (eQTL)` = "public",
  `Disease and traits` = c("public", "ld"), Pathway = "public", Tissue = "public")

interp_summary_visible_evidence <- function() 2L
interp_request_states <- function(public = NULL, spliceai = NULL, ld = NULL) {
  vapply(list(public = public, spliceai = spliceai, ld = ld), function(x)
    if (inherits(x, "error")) "failed" else if (is.null(x)) "before" else "after", character(1))
}
interp_step_label <- function(key, state) {
  label <- interp_request_steps()$label[match(key, interp_request_steps()$key)]
  ifelse(state == "failed", paste(label, "failed"), paste(ifelse(state == "after", "After", "Before"), label))
}
interp_summary_csv <- function(summary, states) {
  out <- summary
  out$on_request_steps <- vapply(summary$question, function(q) {
    steps <- interp_summary_scope()[[q]]; paste(paste0(steps, ": ", states[steps]), collapse = "; ")
  }, character(1))
  out
}

interp_empty_evidence_items <- function() data.frame(text = character(), source = character(), step = character(), role = character(),
  gene_id = character(), gene = character(), clinvar_id = character(), stringsAsFactors = FALSE)
interp_evidence_tag <- function(text, source, step = "local", role = "evidence", gene_id = "", gene = "", clinvar_id = "") {
  text <- as.character(text)
  if (!length(text)) return(text)
  stopifnot(all(step %in% c("local", "public", "spliceai", "ld")), all(role %in% c("evidence", "notice", "context")))
  attr(text, "evidence_items") <- data.frame(text = text, source = rep_len(source, length(text)), step = rep_len(step, length(text)),
    role = rep_len(role, length(text)), gene_id = rep_len(gene_id, length(text)), gene = rep_len(gene, length(text)), clinvar_id = rep_len(clinvar_id, length(text)), stringsAsFactors = FALSE)
  text
}
interp_evidence_items <- function(text) {
  x <- attr(text, "evidence_items")
  if (is.null(x)) x <- if (!length(text)) interp_empty_evidence_items() else
    attr(interp_evidence_tag(text, "unattributed"), "evidence_items")
  x[nzchar(x$text), , drop = FALSE]
}
interp_evidence_combine <- function(...) {
  args <- list(...)
  # Preserve empty strings until the legacy row formatter removes them.
  out <- unlist(lapply(args, as.character), use.names = FALSE)
  items <- lapply(args, function(x) {
    m <- attr(x, "evidence_items")
    if (!is.null(m)) m else if (length(x)) attr(interp_evidence_tag(x, "unattributed"), "evidence_items") else interp_empty_evidence_items()
  })
  attr(out, "evidence_items") <- do.call(rbind, items)
  out
}

# Public-copy provenance is independent of whether its public entity page is linked.
# Notice lines are filtered by the renderer before URL construction.
interp_source_links <- function() data.frame(
  key = c("open_targets", "open_targets_gene", "gtex", "gtex_gene", "gtex_local", "gwas_catalog", "reactome_local", "clinvar_local", "clinvar_api", "cadd"),
  template = c("https://platform.opentargets.org/variant/{variant}", "https://platform.opentargets.org/target/{gene_id}",
    "https://gtexportal.org/home/snp/{rs_id}", "https://gtexportal.org/home/gene/{gene}",
    "https://gtexportal.org/home/snp/{rs_id}", "https://www.ebi.ac.uk/gwas/variants/{rs_id}",
    "https://reactome.org/content/query?q={gene}", "https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}",
    "https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}", "https://cadd.gs.washington.edu/snv/GRCh38-v1.7/{cadd_variant}"),
  token = c("variant", "gene_id", "rs_id", "gene", "rs_id", "rs_id", "gene", "clinvar_id", "clinvar_id", "cadd_variant"),
  spa = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  note = c("", "", "", "", "The GTEx portal shows Release V10; its version may differ from the local GTEx files.", "", "", "", "", ""),
  stringsAsFactors = FALSE)
interp_summary_sources <- function() data.frame(
  key = c("gencode", "cadd", "splice_outcome", "pqtl", "pqtl_finemap", "gwas_catalog", "cohort_ld", "cohort_rna", "topmed", "fhs",
    "gtex_local", "clinvar_local", "gene_disease", "local_pathways", "colocalisation", "spliceai", "open_targets", "open_targets_gene", "gtex", "gtex_gene", "myvariant", "unattributed",
    "reactome_local", "msigdb", "clinvar_api"),
  label = c("GENCODE + reference", "CADD", "Splice outcome", "cis-pQTL scans", "pQTL fine-mapping", "GWAS Catalog", "Cohort LD", "Cohort RNA", "TOPMed", "FHS",
    "GTEx", "ClinVar", "ClinVar gene context", "Pathway annotations", "Colocalisation", "SpliceAI", "Open Targets", "Open Targets genes", "GTEx", "GTEx expression", "MyVariant.info", "Evidence",
    "Reactome", "MSigDB", "ClinVar"),
  kind = c("computed here", "local copy", "computed here", "local", "computed here", "local copy", "computed here", "computed here", "local", "local",
    "local copy", "local copy", "local copy", "local copy", "computed here", "computed here", "external", "external", "external", "external", "via API", "local",
    "local copy", "local copy", "via MyVariant.info API"), stringsAsFactors = FALSE)
interp_single_identifier <- function(x) {
  x <- unique(as.character(x)); x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 1L) x else ""
}
interp_source_url <- function(source, identifiers = list()) {
  registry <- interp_source_links(); i <- match(source, registry$key)
  if (is.na(i)) return(NULL)
  # Allowlisted scalar public identifiers only. Caller fields such as sample,
  # genotype, cohort or arbitrary URLs are never copied into a destination.
  scalar <- function(name) { x <- identifiers[[name]]; if (is.null(x) || length(x) != 1L || is.na(x)) "" else as.character(x) }
  chrom <- sub("^chr", "", scalar("chrom")); pos <- scalar("pos"); ref <- toupper(scalar("ref")); alt <- toupper(scalar("alt"))
  variant <- if (grepl("^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$", chrom) && grepl("^[1-9][0-9]*$", pos) &&
    grepl("^[ACGT]+$", ref) && grepl("^[ACGT]+$", alt)) paste(chrom, pos, ref, alt, sep = "_") else ""
  token <- registry$token[[i]]
  clinvar <- scalar("clinvar_id")
  if (grepl("^[0-9]{1,9}$", clinvar) && as.numeric(clinvar) > 0) clinvar <- sprintf("VCV%09.0f", as.numeric(clinvar))
  cadd <- if (nzchar(variant) && nchar(ref) == 1L && nchar(alt) == 1L) paste0(chrom, ":", pos, "_", ref, "_", alt) else ""
  value <- switch(token, variant = variant, gene_id = sub("\\.[0-9]+$", "", scalar("gene_id")), rs_id = scalar("rs_id"), gene = scalar("gene"), clinvar_id = clinvar, cadd_variant = cadd)
  valid <- switch(token, variant = nzchar(value), gene_id = grepl("^ENSG[0-9]+$", value), rs_id = grepl("^rs[0-9]+$", value), gene = grepl("^[A-Za-z0-9_.-]+$", value), clinvar_id = grepl("^VCV[0-9]{9}$", value) && value != "VCV000000000", cadd_variant = nzchar(value))
  if (!valid) return(NULL)
  sub(paste0("{", token, "}"), utils::URLencode(value, reserved = TRUE), registry$template[[i]], fixed = TRUE)
}
