# Rscript tests/test_interpretation.R [--native]
# Synthetic transcripts, a synthetic reference and recorded-shape API answers:
# no controlled data and no network. --native also builds tiny indexed files and
# needs tabix, bgzip and look.
source("src/interpretation_backend.R", encoding = "UTF-8")
native <- "--native" %in% commandArgs(trailingOnly = TRUE)

expect_error <- function(expr, text = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(e, "error"))
  if (!is.null(text)) stopifnot(grepl(text, conditionMessage(e), fixed = TRUE))
}

# ---- Alleles ----------------------------------------------------------------
v <- interp_variant("chr1", 100L, "g", "a")
stopifnot(v$type == "SNV", v$start1 == 100L, v$end1 == 100L, v$ref_changed == "G", v$alt_changed == "A", v$id == "chr1:100:G:A")
d <- interp_variant("chr1", 100L, "ATCT", "A")
stopifnot(d$type == "deletion", d$start1 == 101L, d$end1 == 103L, d$ref_changed == "TCT", d$alt_changed == "")
i <- interp_variant("chr1", 100L, "A", "ATG")
stopifnot(i$type == "insertion", i$start1 == 101L, i$end1 == 100L, i$alt_changed == "TG")
m <- interp_variant("chr1", 100L, "AT", "AC")
stopifnot(m$type == "SNV", m$start1 == 101L, m$ref_changed == "T", m$alt_changed == "C")
stopifnot(interp_variant("chr1", 100L, "TA", "CA")$start1 == 100L, interp_variant("chr1", 100L, "AC", "GT")$type == "MNV",
  interp_variant("chr1", 100L, "AA", "A")$type == "deletion")
expect_error(interp_variant("chr1", 100L, "A", "A"), "identical")
expect_error(interp_variant("chr1", 100L, "A", "C,G"), "one ALT")
expect_error(interp_variant("chr1", 100L, "A", "<DEL>"), "symbolic")
stopifnot(interp_hgvs_g(v) == "chr1:g.100G>A", interp_hgvs_g(d) == "chr1:g.101_103del", interp_hgvs_g(i) == "chr1:g.100_101insTG",
  interp_hgvs_g(interp_variant("1", 5L, "AT", "A")) == "chr1:g.6del", interp_hgvs_g(interp_variant("chrM", 5L, "A", "G")) == "chrMT:g.5A>G")

# ---- Genetic code -----------------------------------------------------------
std <- interp_codon_table(); mito <- interp_codon_table(TRUE)
stopifnot(length(std) == 64L, std[["ATG"]] == "M", std[["TGG"]] == "W", std[["TAA"]] == "*", std[["TAG"]] == "*", std[["TGA"]] == "*",
  std[["GCT"]] == "A", std[["AGA"]] == "R", std[["ATA"]] == "I", sum(std == "*") == 3L,
  mito[["AGA"]] == "*", mito[["AGG"]] == "*", mito[["ATA"]] == "M", mito[["TGA"]] == "W", sum(mito == "*") == 4L)
stopifnot(interp_most_severe(c("intron_variant", "missense_variant&splice_region_variant")) == "missense_variant",
  identical(interp_so_impact(c("stop_gained", "missense_variant", "synonymous_variant", "intron_variant")), c("HIGH", "MODERATE", "LOW", "MODIFIER")))

# ---- Synthetic annotation ---------------------------------------------------
# Plus-strand coding transcript on chrT: exons 101-160, 201-260, 301-360;
# CDS 121-160, 201-260, 301-317; stop codon 318-320.
genome <- list(chrT = rep("A", 400L), chrU = rep("A", 400L), chrM = rep("A", 60L))
put <- function(chrom, start, text) genome[[chrom]][start + seq_len(nchar(text)) - 1L] <<- strsplit(text, "")[[1L]]
put("chrT", 121L, "ATGGCTTGG")          # Met Ala Trp
put("chrT", 157L, "GCGA"); put("chrT", 201L, "T")   # codon 13 = GCG; codon 14 = A|TA (split by intron 1)
put("chrT", 318L, "TAA")
# Minus-strand transcript on chrU with the same exons; start codon at 338-340.
put("chrU", 338L, "CAT"); put("chrU", 335L, "AGC"); put("chrU", 141L, "TTA")
put("chrM", 1L, "ATAAGCTAA")            # Met Ser stop under the mitochondrial code
fetch <- function(chrom, start1, end1) paste(genome[[chrom]][start1:end1], collapse = "")
tx_row <- function(id, gene, chrom, strand, start, end, tags = "basic", type = "protein_coding")
  data.frame(transcript_id = id, gene_id = paste0("ENSG", gene, ".1"), gene_name = gene, transcript_name = id, transcript_type = type,
    transcript_support_level = "1", tags = tags, chrom = chrom, start1 = start, end1 = end, strand = strand, source = "SYNTHETIC",
    stringsAsFactors = FALSE)
exon_rows <- function(id, gene, chrom, strand, starts, ends) data.frame(transcript_id = id, gene_id = paste0("ENSG", gene, ".1"),
  exon_id = paste0(id, ".", seq_along(starts)), exon_number = as.character(seq_along(starts)), chrom = chrom, start1 = starts, end1 = ends,
  strand = strand, stringsAsFactors = FALSE)
feature_rows <- function(id, gene, chrom, strand, feature, starts, ends, frames) data.frame(transcript_id = id,
  gene_id = paste0("ENSG", gene, ".1"), feature = feature, chrom = chrom, start1 = starts, end1 = ends, strand = strand,
  frame = as.character(frames), stringsAsFactors = FALSE)
annotation_T <- list(
  transcripts = rbind(tx_row("TP", "PLUS", "chrT", "+", 101L, 360L, "basic;Ensembl_canonical;MANE_Select"),
    tx_row("TNC", "PLUS", "chrT", "+", 101L, 260L, "basic", "lncRNA"),
    tx_row("TNF", "PARTIAL", "chrT", "+", 121L, 160L, "cds_start_NF;cds_end_NF")),
  exons = rbind(exon_rows("TP", "PLUS", "chrT", "+", c(101L, 201L, 301L), c(160L, 260L, 360L)),
    exon_rows("TNC", "PLUS", "chrT", "+", c(101L, 201L), c(160L, 260L)), exon_rows("TNF", "PARTIAL", "chrT", "+", 121L, 160L)),
  features = rbind(feature_rows("TP", "PLUS", "chrT", "+", c("CDS", "CDS", "CDS", "stop_codon"), c(121L, 201L, 301L, 318L),
      c(160L, 260L, 317L, 320L), c(0, 2, 2, 0)),
    feature_rows("TNF", "PARTIAL", "chrT", "+", "CDS", 121L, 160L, 2)))
annotation_U <- list(transcripts = tx_row("TM", "MINUS", "chrU", "-", 101L, 360L, "basic;Ensembl_canonical"),
  exons = exon_rows("TM", "MINUS", "chrU", "-", c(101L, 201L, 301L), c(160L, 260L, 360L)),
  features = feature_rows("TM", "MINUS", "chrU", "-", c("CDS", "CDS", "CDS", "stop_codon"), c(301L, 201L, 144L, 141L),
    c(340L, 260L, 160L, 143L), c(0, 2, 2, 0)))
annotation_M <- list(transcripts = tx_row("TMT", "MT-X", "chrM", "+", 1L, 9L, "basic;Ensembl_canonical"),
  exons = exon_rows("TMT", "MT-X", "chrM", "+", 1L, 9L),
  features = feature_rows("TMT", "MT-X", "chrM", "+", c("CDS", "stop_codon"), c(1L, 7L), c(6L, 9L), c(0, 0)))
csq <- function(annotation, chrom, pos, ref, alt, tx = NULL) {
  out <- interp_consequences(interp_variant(chrom, pos, ref, alt), annotation, fetch)
  if (is.null(tx)) out else out[out$transcript_id == tx, , drop = FALSE]
}
one <- function(pos, ref, alt, tx = "TP") csq(annotation_T, "chrT", pos, ref, alt, tx)

# Coding changes on the plus strand.
x <- one(124L, "G", "A")
stopifnot(x$consequence == "missense_variant", x$hgvs_c == "c.4G>A", x$hgvs_p == "p.Ala2Thr", x$codon_change == "Gct/Act",
  x$amino_acids == "A/T", x$protein_position == 2L, x$impact == "MODERATE", x$location == "exon 1/3")
x <- one(126L, "T", "C"); stopifnot(x$consequence == "synonymous_variant", x$hgvs_p == "p.Ala2=", x$amino_acids == "A", x$impact == "LOW")
x <- one(129L, "G", "A"); stopifnot(x$consequence == "stop_gained", x$hgvs_p == "p.Trp3Ter", x$impact == "HIGH")
x <- one(121L, "A", "G"); stopifnot(x$consequence == "start_lost", x$hgvs_p == "p.Met1?", x$hgvs_c == "c.1A>G")
x <- one(318L, "T", "C"); stopifnot(x$consequence == "stop_lost", x$hgvs_p == "p.Ter40Glnext*?", x$hgvs_c == "c.118T>C")
x <- one(320L, "A", "G"); stopifnot(x$consequence == "stop_retained_variant", x$hgvs_p == "p.Ter40=")
# A codon split by an intron: A|TA (Ile) -> A|AA (Lys).
x <- one(201L, "T", "A")
stopifnot(x$consequence == "missense_variant&splice_region_variant", x$hgvs_c == "c.41T>A", x$hgvs_p == "p.Ile14Lys",
  x$codon_change == "aTa/aAa", x$acceptor_distance == -1L)
x <- one(159L, "G", "A"); stopifnot(x$consequence == "synonymous_variant&splice_region_variant", x$donor_distance == -2L, x$hgvs_c == "c.39G>A")
x <- one(160L, "A", "G"); stopifnot(x$consequence == "missense_variant&splice_region_variant", x$hgvs_p == "p.Ile14Val", x$codon_change == "Ata/Gta")
# UTRs and flanks.
x <- one(110L, "A", "G"); stopifnot(x$consequence == "5_prime_UTR_variant", x$hgvs_c == "c.-11A>G")
x <- one(330L, "A", "G"); stopifnot(x$consequence == "3_prime_UTR_variant", x$hgvs_c == "c.*10A>G")
x <- one(50L, "A", "G"); stopifnot(x$consequence == "upstream_gene_variant", x$distance_to_transcript == 51L)
x <- one(400L, "A", "G"); stopifnot(x$consequence == "downstream_gene_variant", x$distance_to_transcript == 40L)
stopifnot(!nrow(interp_consequences(interp_variant("chrT", 6000L, "A", "G"), annotation_T, fetch)))
# Intron 1 is 161-200: every splice term at its defined distance.
terms <- function(pos) one(pos, "A", "G")$consequence
stopifnot(terms(161L) == "splice_donor_variant", terms(162L) == "splice_donor_variant",
  terms(163L) == "splice_donor_region_variant&intron_variant", terms(164L) == "splice_donor_region_variant&intron_variant",
  terms(165L) == "splice_donor_5th_base_variant&intron_variant", terms(166L) == "splice_donor_region_variant&intron_variant",
  terms(167L) == "splice_region_variant&intron_variant", terms(168L) == "splice_region_variant&intron_variant",
  terms(169L) == "intron_variant", terms(183L) == "intron_variant",
  terms(184L) == "splice_polypyrimidine_tract_variant&intron_variant",
  terms(193L) == "splice_region_variant&splice_polypyrimidine_tract_variant&intron_variant",
  terms(198L) == "splice_region_variant&splice_polypyrimidine_tract_variant&intron_variant",
  terms(199L) == "splice_acceptor_variant", terms(200L) == "splice_acceptor_variant")
x <- one(161L, "A", "G"); stopifnot(x$hgvs_c == "c.40+1A>G", x$donor_distance == 1L, x$acceptor_distance == 40L, x$location == "intron 1/2")
x <- one(200L, "A", "G"); stopifnot(x$hgvs_c == "c.41-1A>G", x$acceptor_distance == 1L)
stopifnot(identical(interp_intron_terms(1L, 2L), c("splice_acceptor_variant", "splice_donor_variant", "intron_variant")))
# Non-coding and partial transcripts.
x <- one(124L, "G", "A", "TNC"); stopifnot(x$consequence == "non_coding_transcript_exon_variant", x$hgvs_c == "n.24G>A", x$hgvs_p == "")
x <- one(121L, "A", "G", "TNF"); stopifnot(x$consequence == "coding_sequence_variant", grepl("cds_start_NF", x$notes), grepl("first complete codon", x$notes))
x <- one(124L, "G", "A", "TNF"); stopifnot(x$consequence == "missense_variant", x$protein_position == 1L, x$hgvs_p == "p.Gly1Asp")
# Ordering: most severe first, MANE first within a consequence; gene summary keeps the reference transcript.
all <- one(124L, "G", "A", NULL)
stopifnot(all$transcript_id[[1L]] == "TP", nrow(all) == 3L)
g <- interp_gene_summary(all)
stopifnot(g$gene_name[[1L]] == "PLUS", g$most_severe[[1L]] == "missense_variant", g$transcripts[[1L]] == 2L, g$transcripts_most_severe[[1L]] == 1L,
  grepl("TP (", g$reference_transcript[[1L]], fixed = TRUE), g$reference_hgvs[[1L]] == "c.4G>A p.Ala2Thr")
# REF that disagrees with the reference is reported, not silently used.
x <- one(124L, "C", "A"); stopifnot(grepl("differs from the reference", x$notes))
# Minus strand: alleles are complemented, numbering runs from the right.
x <- csq(annotation_U, "chrU", 337L, "C", "T")
stopifnot(x$consequence == "missense_variant", x$hgvs_c == "c.4G>A", x$hgvs_p == "p.Ala2Thr", x$codon_change == "Gct/Act", x$strand == "-")
x <- csq(annotation_U, "chrU", 300L, "A", "G"); stopifnot(x$consequence == "splice_donor_variant", x$hgvs_c == "c.40+1T>C", x$donor_distance == 1L)
x <- csq(annotation_U, "chrU", 261L, "A", "G"); stopifnot(x$consequence == "splice_acceptor_variant", x$hgvs_c == "c.41-1T>C")
x <- csq(annotation_U, "chrU", 143L, "A", "G"); stopifnot(x$consequence == "stop_lost", x$hgvs_c == "c.118T>C")
x <- csq(annotation_U, "chrU", 400L, "A", "G"); stopifnot(x$consequence == "upstream_gene_variant")
x <- csq(annotation_U, "chrU", 350L, "A", "G"); stopifnot(x$consequence == "5_prime_UTR_variant", x$hgvs_c == "c.-10T>C")
# Mitochondrial code: AGC -> AGA is a stop there, an arginine elsewhere; no flank genes on chrM.
x <- csq(annotation_M, "chrM", 6L, "C", "A"); stopifnot(x$consequence == "stop_gained", x$hgvs_p == "p.Ser2Ter")
x <- csq(annotation_M, "chrM", 3L, "A", "G"); stopifnot(x$consequence == "start_retained_variant")
stopifnot(!nrow(csq(annotation_M, "chrM", 20L, "A", "G")))
# Indels are classified by overlap only.
x <- csq(annotation_T, "chrT", 123L, "GG", "G", "TP"); stopifnot(x$consequence == "frameshift_variant", grepl("non-SNV", x$notes), x$hgvs_p == "")
x <- csq(annotation_T, "chrT", 123L, "GGCT", "G", "TP"); stopifnot(x$consequence == "inframe_deletion")
x <- csq(annotation_T, "chrT", 124L, "G", "GAAA", "TP"); stopifnot(x$consequence == "inframe_insertion")
x <- csq(annotation_T, "chrT", 160L, "AA", "A", "TP"); stopifnot(grepl("splice_donor_variant", x$consequence))
x <- csq(annotation_T, "chrT", 170L, "A", "AT", "TP"); stopifnot(x$consequence == "intron_variant")
x <- csq(annotation_T, "chrT", 109L, "AA", "A", "TP"); stopifnot(x$consequence == "5_prime_UTR_variant")
stopifnot(!nrow(interp_consequences(v, list(transcripts = NULL), fetch)), !nrow(interp_gene_summary(interp_empty_consequences())))

# ---- Nested JSON helpers ----------------------------------------------------
stopifnot(interp_get(list(a = list(b = list(1, 2))), "a", "b", 2L) == 2, is.na(interp_get(list(a = NULL), "a", "b")),
  interp_get(list(), "x", default = "d") == "d", interp_join(c("a", "", NA, "a", "b")) == "a; b",
  interp_join(letters[1:5], ", ", max = 2) == "a, b, +3 more", interp_join(NULL) == "")
p <- interp_p_value(5, -400); stopifnot(p[["p"]] == 0, abs(p[["log10"]] - (log10(5) - 400)) < 1e-9, is.na(interp_p_value(NULL, 1)[["p"]]))
stopifnot(identical(interp_fmt_p(c(0, 1e-320, 3.2e-15, NA)), c("< 1e-300", "< 1e-300", "3.2e-15", "NA")))
stopifnot(interp_bare_chrom("chrM") == "MT", interp_bare_chrom("chr7") == "7", interp_bare_gene("ENSG01.12") == "ENSG01")

# ---- Open Targets: variant --------------------------------------------------
credible <- function(type, mant, exp, pip, study) list(studyLocusId = "x", studyType = type, isTransQtl = isTRUE(study$trans), pValueMantissa = mant,
  pValueExponent = exp, beta = .5, finemappingMethod = "SuSiE", locus = list(rows = list(list(posteriorProbability = pip, is95CredibleSet = TRUE))),
  study = study, l2GPredictions = list(rows = list(list(score = .9, target = list(approvedSymbol = "GENE1")))))
ot_core <- list(data = list(variant = list(id = "1_100_G_A", rsIds = list("rs1"), hgvsId = "1:g.100G>A", variantDescription = "Synthetic",
  mostSevereConsequence = list(label = "missense_variant"),
  alleleFrequencies = list(list(populationName = "nfe_adj", alleleFrequency = .1)),
  variantEffect = list(list(method = "AlphaMissense", assessment = "likely_pathogenic", score = .9, normalisedScore = .8, target = list(approvedSymbol = "GENE1"))),
  dbXrefs = list(list(id = "VCV1", source = "clinvar"), list(id = "rs1", source = "ensembl_variation")),
  evidences = list(count = 2L, rows = list(
    list(datasourceId = "eva", score = .9, clinicalSignificances = list("pathogenic"), confidence = "criteria provided", disease = list(name = "Disease A"),
      target = list(approvedSymbol = "GENE1"), literature = list("1", "2")),
    list(datasourceId = "gwas_credible_sets", score = .5, disease = list(name = "Trait"), target = list(approvedSymbol = "GENE1")))),
  pharmacogenomics = list(list(genotype = "AA", phenotypeText = "slow", pgxCategory = "metabolism", evidenceLevel = "1A",
    drugs = list(list(drugFromSource = "drugA"), list(drugFromSource = "drugB")), target = list(approvedSymbol = "GENE1"))))))
gwas_rows <- list(credible("gwas", 2, -9, .05, list(id = "S1", traitFromSource = "Trait low", nSamples = 100L, diseases = list(list(name = "D low")))),
  credible("gwas", 1, -400, .97, list(id = "S2", traitFromSource = "Trait high", nSamples = 5000L, publicationFirstAuthor = "A", pubmedId = "9",
    diseases = list(list(name = "D high"), list(name = "D other")))))
qtl_rows <- list(credible("pqtl", 1, -30, .95, list(id = "Q1", projectId = "UKB_PPP", target = list(approvedSymbol = "GENE1", id = "ENSG1"),
    biosample = list(biosampleName = "blood plasma"))),
  credible("eqtl", 3, -8, .2, list(id = "Q2", projectId = "GTEx", target = list(approvedSymbol = "GENE2", id = "ENSG2"), biosample = list(biosampleName = "lung"))),
  credible("pqtl", 1, -12, .99, list(id = "Q3", trans = TRUE, target = list(approvedSymbol = "FAR", id = "ENSG3"), biosample = list(biosampleName = "blood plasma"))))
calls <- 0L; pages <- list()
http_ot <- function(url, body = NULL, timeout = 45) {
  calls <<- calls + 1L
  req <- jsonlite::fromJSON(body, simplifyVector = FALSE)
  stopifnot(identical(url, interp_open_targets_url()), identical(req$variables$id, "1_100_G_A"))
  if (is.null(req$variables$index)) return(ot_core)
  gwas <- grepl("studyTypes:[gwas]", req$query, fixed = TRUE)
  pages[[length(pages) + 1L]] <<- c(kind = if (gwas) "gwas" else "qtl", index = req$variables$index)
  rows <- if (gwas) gwas_rows else qtl_rows
  # GWAS arrives one row per page to exercise paging.
  take <- if (gwas) rows[req$variables$index + 1L] else rows
  list(data = list(variant = list(credibleSets = list(count = length(rows), rows = take))))
}
ot <- interp_open_targets_variant(v, http_ot, page_size = 1L)
stopifnot(ot$found, ot$summary$rs_ids == "rs1", ot$summary$clinvar_ids == "VCV1", ot$summary$gwas_total == 2L, ot$summary$qtl_total == 3L,
  nrow(ot$gwas) == 2L, ot$gwas$trait[[1L]] == "Trait high", ot$gwas$diseases[[1L]] == "D high; D other", ot$gwas$p_value[[1L]] == 0,
  abs(ot$gwas$log10_p[[1L]] + 400) < 1e-9, ot$gwas$posterior_probability[[1L]] == .97, ot$gwas$l2g_gene[[1L]] == "GENE1",
  identical(ot$qtl$qtl_type, c("pqtl", "pqtl", "eqtl")), identical(ot$qtl$trans, c(TRUE, FALSE, FALSE)), ot$qtl$biosample[[3L]] == "lung",
  nrow(ot$clinical) == 1L, ot$clinical$source == "eva", ot$clinical$clinical_significance == "pathogenic", ot$clinical$literature == "1; 2",
  ot$predictors$method == "AlphaMissense", ot$frequencies$allele_frequency == .1, ot$pharmacogenomics$drugs == "drugA; drugB",
  sum(vapply(pages, function(p) p[["kind"]] == "gwas", logical(1))) == 2L)
absent <- interp_open_targets_variant(v, function(url, body = NULL, timeout = 45) list(data = list(variant = NULL)))
stopifnot(!absent$found, !nrow(absent$gwas), !nrow(absent$qtl), identical(names(absent$gwas), names(interp_ot_empty()$gwas)))
expect_error(interp_open_targets_variant(v, function(url, body = NULL, timeout = 45) list(errors = list(list(message = "too expensive")))), "too expensive")

# ---- Open Targets: genes ----------------------------------------------------
targets <- list(data = list(targets = list(list(id = "ENSG1", approvedSymbol = "GENE1", approvedName = "gene one", biotype = "protein_coding",
  functionDescriptions = list("Does things."), proteinIds = list(list(id = "P1", source = "uniprot_swissprot"), list(id = "X", source = "ensembl_PRO")),
  subcellularLocations = list(list(location = "Membrane")),
  pathways = list(list(pathway = "Lipid transport", topLevelTerm = "Transport", pathwayId = "R-HSA-1"),
    list(pathway = "Lipid transport", topLevelTerm = "Transport", pathwayId = "R-HSA-1")),
  geneOntology = list(list(aspect = "P", evidence = "IDA", term = list(id = "GO:1", label = "transport")),
    list(aspect = "F", evidence = "IEA", term = list(id = "GO:2", label = "binding"))),
  geneticConstraint = list(list(constraintType = "syn", oeUpper = 1), list(constraintType = "lof", oeUpper = .3, upperBin6 = 0)),
  associatedDiseases = list(count = 12L, rows = list(list(score = .8, disease = list(id = "D1", name = "Disease A"))))))))
seen <- NULL
genes <- interp_open_targets_genes(c("ENSG1.4", "ENSG1.5", "not-a-gene"), function(url, body = NULL, timeout = 45) {
  seen <<- jsonlite::fromJSON(body, simplifyVector = FALSE)$variables$ids; targets })
stopifnot(identical(unlist(seen), "ENSG1"), genes$genes$uniprot == "P1", genes$genes$lof_oe_upper == .3, genes$genes$associated_diseases == 12L,
  nrow(genes$pathways) == 1L, genes$pathways$top_level == "Transport", identical(sort(genes$go$aspect), c("Biological process", "Molecular function")),
  genes$diseases$disease == "Disease A")
stopifnot(!nrow(interp_open_targets_genes("nonsense", function(...) stop("must not be called"))$genes))

# ---- GTEx -------------------------------------------------------------------
gtex_urls <- character()
http_gtex <- function(url, body = NULL, timeout = 45) {
  gtex_urls <<- c(gtex_urls, url)
  page <- function(data, pages = 1L) list(data = data, paging_info = list(numberOfPages = pages))
  if (grepl("dataset/variant", url, fixed = TRUE)) return(page(list(list(variantId = "chr1_100_G_A_b38"))))
  if (grepl("singleTissueEqtl", url, fixed = TRUE)) {
    first <- grepl("page=0", url, fixed = TRUE)
    return(page(list(list(tissueSiteDetailId = if (first) "Whole_Blood" else "Lung", geneSymbol = "GENE1", gencodeId = "ENSG1.9",
      pValue = if (first) 1e-9 else 1e-20, nes = .4)), 2L))
  }
  if (grepl("singleTissueSqtl", url, fixed = TRUE)) return(page(list(list(tissueSiteDetailId = "Liver", geneSymbol = "GENE1", gencodeId = "ENSG1.9",
    phenotypeId = "chr1:90:120:clu_1:ENSG1.9", pValue = 1e-12, nes = -.9))))
  if (grepl("reference/gene", url, fixed = TRUE)) return(page(list(list(gencodeId = "ENSG1.9"))))
  if (grepl("medianGeneExpression", url, fixed = TRUE)) return(page(list(
    list(geneSymbol = "GENE1", tissueSiteDetailId = "Whole_Blood", median = 1), list(geneSymbol = "GENE1", tissueSiteDetailId = "Liver", median = 99),
    list(geneSymbol = "GENE1", tissueSiteDetailId = "Lung", median = 9))))
  stop("unexpected URL ", url)
}
gtex <- interp_gtex(v, "ENSG1.4", http_gtex)
stopifnot(gtex$variant_in_gtex, gtex$dataset == "gtex_v10", nrow(gtex$qtl) == 3L, gtex$qtl$tissue[[1L]] == "Lung", gtex$qtl$qtl_type[[2L]] == "sQTL",
  "Whole Blood" %in% gtex$qtl$tissue, nrow(gtex$expression) == 3L, gtex$genes_mapped == "ENSG1.9",
  any(grepl("variantId=chr1_100_G_A_b38", gtex_urls, fixed = TRUE)), any(grepl("genomeBuild=GRCh38%2Fhg38", gtex_urls, fixed = TRUE)))
spec <- interp_tissue_specificity(gtex$expression)
stopifnot(spec$top_tissue == "Liver", spec$whole_blood_rank == 3L, spec$tissues == 3L, spec$tau > 0, spec$tau < 1,
  !nrow(interp_tissue_specificity(NULL)))
uniform <- data.frame(gene = "G", tissue = c("a", "b", "Whole Blood"), median_tpm = 5)
stopifnot(abs(interp_tissue_specificity(uniform)$tau) < 1e-12)
not_in <- interp_gtex(v, character(), function(url, body = NULL, timeout = 45) list(data = list(), paging_info = list(numberOfPages = 0L)))
stopifnot(identical(not_in$variant_in_gtex, FALSE), !nrow(not_in$qtl))

# ---- MyVariant.info ---------------------------------------------------------
myv <- interp_parse_myvariant(list(dbsnp = list(rsid = "rs1"), clinvar = list(variant_id = 7L, rcv = list(
  list(accession = "RCV1", clinical_significance = "Pathogenic", conditions = list(name = "Disease A"), review_status = "criteria provided"),
  list(accession = "RCV2", clinical_significance = "Benign", conditions = list(list(name = "X"), list(name = "Y")), review_status = "none"))),
  gnomad_genome = list(af = list(af = .012)), dbnsfp = list(revel = list(score = list(.81, .5)))))
stopifnot(myv$found, myv$rs_id == "rs1", myv$clinvar_variant_id == "7", myv$gnomad_af == .012, myv$revel == .81, nrow(myv$clinvar) == 2L,
  myv$clinvar$condition[[2L]] == "X; Y")
single <- interp_parse_myvariant(list(clinvar = list(rcv = list(accession = "RCV9", clinical_significance = "risk factor", conditions = list(name = "Z")))))
stopifnot(nrow(single$clinvar) == 1L, single$clinvar$accession == "RCV9", !interp_parse_myvariant(NULL)$found)
url_seen <- NULL
invisible(interp_myvariant(v, function(url, body = NULL, timeout = 45) { url_seen <<- url; NULL }))
stopifnot(grepl("chr1%3Ag.100G%3EA", url_seen, fixed = TRUE), grepl("assembly=hg38", url_seen, fixed = TRUE))
# One failing source does not take the others down.
public <- interp_public_evidence(v, "ENSG1.4", function(url, body = NULL, timeout = 45) {
  if (grepl("gtexportal", url, fixed = TRUE)) stop("GTEx is down")
  if (grepl("myvariant", url, fixed = TRUE)) return(NULL)
  req <- jsonlite::fromJSON(body, simplifyVector = FALSE)
  if (!is.null(req$variables$ids)) return(targets)
  http_ot(url, body) })
stopifnot(inherits(public$gtex, "error"), public$open_targets$found, nrow(public$genes$genes) == 1L, !public$myvariant$found,
  public$variant == "chr1:100:G:A", identical(public$gene_ids, "ENSG1"))

# ---- SpliceAI parsing -------------------------------------------------------
vcf <- c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  "chr1\t100\t.\tG\tA\t.\t.\tSpliceAI=A|GENE1|0.00|0.02|0.42|0.37|-31|496|24|-1,A|GENE2|0.91|0.00|0.00|0.00|5|1|2|3")
sai <- interp_parse_spliceai(vcf, v)
stopifnot(nrow(sai) == 2L, sai$donor_gain[[1L]] == .42, sai$donor_gain_pos1[[1L]] == 124L, sai$donor_loss_pos1[[1L]] == 99L,
  sai$max_delta[[2L]] == .91, sai$acceptor_gain_pos1[[2L]] == 105L)
stopifnot(!nrow(interp_parse_spliceai(c("#h", "chr1\t100\t.\tG\tA\t.\t.\t."), v)), !nrow(interp_parse_spliceai("#h", v)))
introns <- data.frame(strand = c("+", "-"), donor1 = c(100L, 104L), acceptor1 = c(300L, 104L))
sites <- interp_spliceai_sites(sai, c("+", "-"), introns)
stopifnot(nrow(sites) == 3L, sites$effect[[1L]] == "acceptor gain", sites$site_position[[1L]] == 106L, sites$annotated_site[[1L]] == "no",
  # Plus-strand donor: the intron starts one base after the exon end.
  sites$site_position[sites$effect == "donor gain"] == 125L, sites$site_position[sites$effect == "donor loss"] == 100L,
  sites$annotated_site[sites$effect == "donor loss"] == "yes", sites$annotated_site[sites$effect == "donor gain"] == "no")
stopifnot(!nrow(interp_spliceai_sites(interp_spliceai_empty(), character())), interp_spliceai_sites(sai, c(NA, NA))$annotated_site[[1L]] == "unknown")
expect_error(interp_spliceai(v, "/nonexistent.fa", list(spliceai = "/nonexistent")), "not available")

# ---- Local sources ----------------------------------------------------------
old_env <- Sys.getenv(c("REGSHINY_CADD_SNV", "REGSHINY_FHS_SQTL"), unset = NA)
Sys.setenv(REGSHINY_CADD_SNV = "/tmp/custom_cadd.tsv.gz", REGSHINY_FHS_SQTL = "/a.csv;/b.csv")
src <- interp_local_sources()
stopifnot(src$cadd_snv == "/tmp/custom_cadd.tsv.gz", identical(unname(src$fhs_sqtl), c("/a.csv", "/b.csv")),
  identical(names(src$fhs_sqtl), c("discovery", "replication")))
for (n in names(old_env)) if (is.na(old_env[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(old_env[[n]]), n))
none <- list(tabix = "", look = "", cadd_snv = "/missing", cadd_indel = "/missing", vep_cache = "/missing", topmed_qtl = "/missing",
  fhs_sqtl = c(discovery = "/missing", replication = "/missing"), fhs_eqtl = "/missing", gwas_catalog = "/missing",
  clinvar_vcf = "/missing", reactome = "/missing", msigdb = "/missing", tissue_model = c(eQTL = "/missing", sQTL = "/missing"), spliceai = "/missing")
status <- interp_local_status(none)
stopifnot(nrow(status) == 11L, !any(status$available))
# A source list written before a source existed simply reports it as unavailable.
stopifnot(!any(interp_local_status(none[c("tabix", "look", "cadd_snv", "cadd_indel", "vep_cache", "topmed_qtl", "fhs_sqtl", "fhs_eqtl", "spliceai")])$available))
gone <- interp_local_evidence(v, none)
stopifnot(!gone$cadd$available, !gone$known$available, !gone$fhs$sqtl_available, !gone$fhs$eqtl_available, !gone$topmed$available,
  !nrow(gone$topmed$pairs), !nrow(gone$fhs$sqtl), !gone$gwas_catalog$available, !gone$clinvar$available,
  !nrow(gone$gwas_catalog$associations), !nrow(gone$clinvar$records))
no_paths <- interp_local_pathways(data.frame(gene_id = "ENSG1", gene = "G", relation = "overlaps", stringsAsFactors = FALSE), none)
stopifnot(!no_paths$available, !nrow(no_paths$reactome), !nrow(no_paths$gene_sets))

# Replication needs the same transcript and direction in both sets.
sq <- data.frame(set = c("discovery", "replication", "discovery", "replication", "discovery"), feature = c("T1", "T1", "T2", "T2", "T3"),
  gene = "G", cis = TRUE, beta_alt = c(.5, .4, .2, -.2, .9), t_value = 1, r_squared = .1, log10_p = c(-30, -12, -8, -9, -50),
  alt_frequency = .1, n = 100L, stringsAsFactors = FALSE)
rep <- interp_fhs_replicated(sq)
stopifnot(nrow(rep) == 1L, rep$feature == "T1", rep$log10_p_replication == -12, !nrow(interp_fhs_replicated(sq[sq$set == "discovery", ])))

# Cross-context comparison.
pair <- function(type, context, phenotype, slope, se, p, significant, gene_id = "ENSG1.1") data.frame(qtl_type = type, context = context,
  phenotype = phenotype, gene_id = gene_id, gene = "", alt_frequency = .1, n_estimated = 500L, slope_alt = slope, slope_se = se,
  p_nominal = p, p_threshold = 1e-5, significant = significant, pip = NA_real_, credible_set = "", tss_distance = 10L, stringsAsFactors = FALSE)
pairs <- rbind(
  pair("sQTL", "Lung", "chr1:10:90:clu_5_+:ENSG1.1", 2, .1, 1e-50, TRUE), pair("sQTL", "PBMC", "chr1:10:90:clu_77_+:ENSG1.1", 2.1, .2, 1e-20, TRUE),
  pair("eQTL", "Lung", "ENSG2.1", .5, .05, 1e-20, TRUE, "ENSG2.1"), pair("eQTL", "PBMC", "ENSG2.1", .45, .3, .13, FALSE, "ENSG2.1"),
  pair("eQTL", "Lung", "ENSG3.1", .6, .05, 1e-25, TRUE, "ENSG3.1"), pair("eQTL", "PBMC", "ENSG3.1", .0, .05, .9, FALSE, "ENSG3.1"),
  pair("eQTL", "Lung", "ENSG4.1", .6, .05, 1e-25, TRUE, "ENSG4.1"), pair("eQTL", "PBMC", "ENSG4.1", -.6, .05, 1e-25, TRUE, "ENSG4.1"),
  pair("eQTL", "Lung", "ENSG5.1", .6, .05, 1e-25, TRUE, "ENSG5.1"), pair("eQTL", "Lung", "ENSG6.1", .1, .05, .04, FALSE, "ENSG6.1"))
stopifnot(interp_phenotype_key("sQTL", "chr1:10:90:clu_5_+:ENSG1.1") == "chr1:10:90:ENSG1.1", interp_phenotype_key("eQTL", "ENSG2.1") == "ENSG2.1")
h <- interp_context_heterogeneity(pairs)
pattern <- stats::setNames(h$pattern, h$gene_id)
stopifnot(nrow(h) == 5L, !"ENSG6.1" %in% h$gene_id, h$contexts_tested[h$gene_id == "ENSG1.1"] == 2L,
  pattern[["ENSG1.1"]] == "shared: significant in every tested context",
  pattern[["ENSG2.1"]] == "consistent effect size; significant in 1 of 2 (the others may lack power)",
  pattern[["ENSG3.1"]] == "context-specific: significant only in Lung, effect differs elsewhere",
  pattern[["ENSG4.1"]] == "opposite directions between contexts", pattern[["ENSG5.1"]] == "tested in one context only",
  is.na(h$q_p_value[h$gene_id == "ENSG5.1"]), h$q_p_value[h$gene_id == "ENSG3.1"] < 1e-10, h$i_squared[h$gene_id == "ENSG3.1"] > .9)
named <- interp_name_genes(pairs, data.frame(gene_id = c("ENSG1.7", "ENSG3.1"), gene_name = c("ONE", "THREE")))
stopifnot(named$gene[[1L]] == "ONE", named$gene[[5L]] == "THREE", named$gene[[3L]] == "", !nrow(interp_context_heterogeneity(NULL)))

# ---- Cohort RNA evidence ----------------------------------------------------
bins <- data.frame(site_id = c("donor:+:chr1:100", "donor:+:chr1:115", "donor:+:chr1:150", "donor:+:chr1:170"), status = c("OK", "OK", "OK", "LOW_MAC"),
  beta = c(-.3, -.29, .01, NA), p_permutation = c(1e-4, 1e-4, .6, NA), n = c(900L, 900L, 880L, 3L),
  bin_check = c("OK", "EXON_BIN_IN_OTHER_TRANSCRIPT_INTRON;NO_DEPTH_DROP_AT_SITE", "OK", "OK"), stringsAsFactors = FALSE)
junc <- data.frame(chrom = "chr1", intron_start1 = c(100L, 100L), intron_end1 = c(300L, 400L), status = "OK", beta = c(-1.2, .1),
  p_permutation = c(2e-4, .2), n = 900L, stringsAsFactors = FALSE)
rna <- interp_rna_evidence(bins, junc)
stopifnot(rna$available, rna$supported, identical(rna$sites, "donor:+:chr1:100"), identical(rna$junctions, "chr1:100-300"),
  any(grepl("set aside", rna$details)), length(rna$details) == 3L)
quiet <- interp_rna_evidence(bins[3:4, ], NULL)
stopifnot(quiet$available, !quiet$supported, !interp_rna_evidence(NULL, NULL)$available, !interp_rna_evidence(simpleError("x"), NULL)$available)

# ---- Evidence summary rules -------------------------------------------------
level <- function(s, question) s$level[s$question == question]
missense <- one(124L, "G", "A", NULL)
s0 <- interp_evidence_summary(missense)
stopifnot(identical(s0$question, c("Protein sequence", "Protein abundance (pQTL)", "Splicing", "Gene expression (eQTL)", "Disease and traits", "Pathway", "Tissue")),
  level(s0, "Protein sequence") == "Moderate", level(s0, "Protein abundance (pQTL)") == "Not assessed", level(s0, "Splicing") == "None found",
  level(s0, "Gene expression (eQTL)") == "Not assessed", level(s0, "Disease and traits") == "Not assessed", level(s0, "Pathway") == "Not assessed",
  level(s0, "Tissue") == "Not assessed", all(s0$level %in% c(interp_levels(), "Context")))
stopifnot(level(interp_evidence_summary(one(129L, "G", "A", NULL)), "Protein sequence") == "Strong",
  level(interp_evidence_summary(one(110L, "A", "G", NULL)), "Protein sequence") == "None found",
  level(interp_evidence_summary(interp_empty_consequences()), "Protein sequence") == "None found")
# An essential splice site in a UTR intron is a splicing finding, not a protein one.
utr <- transform(one(161L, "A", "G", NULL), hgvs_c = "c.*2+1A>G")
stopifnot(level(interp_evidence_summary(utr), "Protein sequence") == "Limited", level(interp_evidence_summary(utr), "Splicing") == "Moderate",
  level(interp_evidence_summary(one(161L, "A", "G", NULL)), "Protein sequence") == "Strong")
# A truncating change only off the reference transcript is Moderate.
off <- one(129L, "G", "A", NULL); off$transcript_flags[off$transcript_id == "TP"] <- "basic"
off <- rbind(off, transform(off[off$transcript_id == "TP", ], transcript_id = "REF", transcript_flags = "MANE_Select", consequence = "intron_variant",
  most_severe = "intron_variant", impact = "MODIFIER"))
stopifnot(level(interp_evidence_summary(off), "Protein sequence") == "Moderate")
pub <- list(open_targets = ot, genes = genes, gtex = gtex, myvariant = myv)
s1 <- interp_evidence_summary(missense, NULL, pub)
stopifnot(level(s1, "Protein abundance (pQTL)") == "Strong", grepl("GENE1 in blood plasma (cis)", s1$finding[[2L]], fixed = TRUE),
  grepl("predictors call it damaging", s1$finding[[1L]], fixed = TRUE), grepl("REVEL 0.81", s1$evidence[[1L]], fixed = TRUE),
  # A GTEx sQTL alone is an association without a pointer to this variant.
  level(s1, "Splicing") == "Moderate", level(s1, "Gene expression (eQTL)") == "Moderate", level(s1, "Disease and traits") == "Strong",
  level(s1, "Pathway") == "Context", grepl("GENE1: Transport", s1$finding[[6L]], fixed = TRUE), level(s1, "Tissue") == "Context",
  grepl("highest in Liver", s1$evidence[[7L]], fixed = TRUE))
# Trans-only pQTL sets cannot be Strong; a set with small PIP is Limited.
trans_only <- pub; trans_only$open_targets$qtl <- ot$qtl[ot$qtl$trans, , drop = FALSE]
stopifnot(level(interp_evidence_summary(missense, NULL, trans_only), "Protein abundance (pQTL)") == "Moderate")
weak <- pub; weak$open_targets$qtl <- transform(ot$qtl[!ot$qtl$trans & ot$qtl$qtl_type == "pqtl", ], posterior_probability = .02)
stopifnot(level(interp_evidence_summary(missense, NULL, weak), "Protein abundance (pQTL)") == "Limited")
no_pqtl <- pub; no_pqtl$open_targets$qtl <- ot$qtl[ot$qtl$qtl_type == "eqtl", , drop = FALSE]
stopifnot(level(interp_evidence_summary(missense, NULL, no_pqtl), "Protein abundance (pQTL)") == "None found")
# Disease: benign next to pathogenic is not Strong on ClinVar alone; GWAS PIP drives the rest.
conflict <- pub; conflict$open_targets$gwas <- ot$gwas[ot$gwas$posterior_probability < .1, , drop = FALSE]
stopifnot(level(interp_evidence_summary(missense, NULL, conflict), "Disease and traits") == "Moderate")
vus <- conflict; vus$myvariant$clinvar <- transform(myv$clinvar[1L, ], clinical_significance = "Uncertain significance"); vus$open_targets$pharmacogenomics <- ot$pharmacogenomics[0, ]
stopifnot(level(interp_evidence_summary(missense, NULL, vus), "Disease and traits") == "Limited")
clean <- vus; clean$myvariant$clinvar <- myv$clinvar[0, ]; clean$open_targets$gwas <- ot$gwas[0, ]
stopifnot(level(interp_evidence_summary(missense, NULL, clean), "Disease and traits") == "None found")
# Splicing: prediction alone, association alone, both.
intronic <- one(175L, "A", "G", NULL)
sai_strong <- list(scores = transform(sai[1L, ], donor_gain = .8, max_delta = .8))
sai_weak <- list(scores = transform(sai[1L, ], donor_gain = .25, donor_loss = .1, max_delta = .25))
stopifnot(level(interp_evidence_summary(intronic), "Splicing") == "None found",
  level(interp_evidence_summary(intronic, spliceai = sai_weak), "Splicing") == "Limited",
  level(interp_evidence_summary(intronic, spliceai = sai_strong), "Splicing") == "Moderate",
  level(interp_evidence_summary(intronic, rna = rna), "Splicing") == "Moderate",
  level(interp_evidence_summary(intronic, spliceai = sai_strong, rna = rna), "Splicing") == "Strong",
  level(interp_evidence_summary(intronic, spliceai = sai_weak, rna = rna), "Splicing") == "Moderate",
  level(interp_evidence_summary(one(161L, "A", "G", NULL)), "Splicing") == "Moderate",
  level(interp_evidence_summary(one(161L, "A", "G", NULL), rna = rna), "Splicing") == "Strong",
  level(interp_evidence_summary(one(167L, "A", "G", NULL)), "Splicing") == "Limited")
topmed <- list(available = TRUE, pairs = pairs, contexts = data.frame(qtl_type = c("eQTL", "eQTL", "sQTL", "sQTL"), context = c("Lung", "PBMC", "Lung", "PBMC"),
  tested = TRUE, pairs_tested = 3L, significant_pairs = 1L, min_p = 1e-9))
local <- list(cadd = list(available = TRUE, found = TRUE, phred = 31, raw = 4), known = list(available = TRUE, found = TRUE, rs_id = "rs1",
  clinvar_significance = "", frequencies = data.frame()), fhs = list(sqtl = sq, eqtl = interp_fhs_empty()), topmed = topmed)
s2 <- interp_evidence_summary(intronic, local)
stopifnot(level(s2, "Splicing") == "Moderate", level(s2, "Gene expression (eQTL)") == "Moderate", level(s2, "Tissue") == "Context",
  grepl("TOPMed sQTL significant in 2 of 2 contexts", s2$evidence[[3L]], fixed = TRUE), grepl("FHS whole-blood isoform ratio: 1 transcript", s2$evidence[[3L]], fixed = TRUE),
  grepl("intron 10-90", s2$evidence[[7L]], fixed = TRUE), grepl("CADD PHRED 31.0", s2$evidence[[1L]], fixed = TRUE),
  level(s2, "Disease and traits") == "Not assessed")
fine <- local; fine$topmed$pairs$pip[fine$topmed$pairs$qtl_type == "sQTL"][[1L]] <- .9
fine$topmed$pairs$pip[fine$topmed$pairs$gene_id == "ENSG2.1"][[1L]] <- .7
s3 <- interp_evidence_summary(intronic, fine)
stopifnot(level(s3, "Splicing") == "Strong", level(s3, "Gene expression (eQTL)") == "Strong")
# The offline ClinVar snapshot is used when the public query was not run.
offline <- local; offline$known$clinvar_significance <- "pathogenic, likely_pathogenic"
stopifnot(level(interp_evidence_summary(intronic, offline), "Disease and traits") == "Strong")
# Errors from a source read as "not available", never as a finding.
broken <- list(cadd = simpleError("x"), known = simpleError("x"), fhs = simpleError("x"), topmed = simpleError("x"))
s4 <- interp_evidence_summary(intronic, broken, list(open_targets = simpleError("x"), genes = simpleError("x"), gtex = simpleError("x"), myvariant = simpleError("x")))
stopifnot(level(s4, "Protein abundance (pQTL)") == "Not assessed", level(s4, "Gene expression (eQTL)") == "Not assessed",
  level(s4, "Disease and traits") == "Not assessed", level(s4, "Tissue") == "Not assessed")

# Genes the variant can act through: overlapped genes first, QTL targets after, one row per gene.
lg_pairs <- rbind(pair("sQTL", "Lung", "chr1:10:90:clu_5_+:ENSGPLUS.1", 2, .1, 1e-50, TRUE, "ENSGPLUS.1"), pair("eQTL", "Lung", "ENSG9.2", .5, .05, 1e-20, TRUE, "ENSG9.2"),
  pair("eQTL", "PBMC", "ENSG9.2", .5, .05, 1e-20, TRUE, "ENSG9.2"), pair("eQTL", "Lung", "ENSG8.1", .1, .05, .04, FALSE, "ENSG8.1"))
lg_pairs$gene[lg_pairs$gene_id == "ENSG9.2"] <- "NINE"
lg <- interp_linked_genes(missense, list(topmed = list(available = TRUE, pairs = lg_pairs)))
stopifnot(nrow(lg) == 3L, identical(lg$gene, c("PARTIAL", "PLUS", "NINE")), lg$gene_id[[2L]] == "ENSGPLUS",
  lg$relation[[1L]] == "overlaps", lg$relation[[2L]] == "overlaps; sQTL target (TOPMed)", lg$relation[[3L]] == "eQTL target (TOPMed)",
  !"ENSG8" %in% lg$gene_id, nrow(interp_linked_genes(interp_empty_consequences())) == 0L, nrow(interp_linked_genes(missense)) == 2L)

# ---- Cross-context sharing model -----------------------------------------------
set.seed(11)
mctx <- c("A", "B", "C"); mV <- matrix(.3, 3L, 3L); diag(mV) <- 1
canon <- interp_mash_canonical(mctx)
stopifnot(identical(names(canon), c("identity", "equal_effects", "simple_het_0.25", "simple_het_0.5", "simple_het_0.75", "only_A", "only_B", "only_C")),
  canon$only_B["B", "B"] == 1, sum(canon$only_B) == 1, all(canon$equal_effects == 1), canon$simple_het_0.5["A", "C"] == .5,
  identical(interp_mash_class(c("null", "only_B", "equal_effects", "simple_het_0.5", "identity", "learned_all")),
    c("no effect", "one context", "shared", "shared", "independent", "learned")))
expect_error(interp_mash_model(mctx, mV, canon, pi = c(.5, .6)), "length")
# Simulated pairs: 70% no effect, 20% the same effect everywhere, 10% in B only; 15% of the z missing.
simulate <- function(n) {
  kind <- sample(c("null", "shared", "onlyB"), n, TRUE, c(.7, .2, .1)); z <- matrix(0, n, 3L)
  z[kind == "shared", ] <- stats::rnorm(sum(kind == "shared"), 0, 4); z[kind == "onlyB", 2L] <- stats::rnorm(sum(kind == "onlyB"), 0, 5)
  zh <- z + matrix(stats::rnorm(n * 3L), n, 3L) %*% chol(mV); zh[sample(length(zh), round(length(zh) * .15))] <- NA
  list(z = z, zh = zh[rowSums(is.finite(zh)) > 0L, , drop = FALSE])
}
sim <- simulate(1500L)
vh <- interp_mash_null_correlation(sim$zh)
# Estimated among |z| < 2 only, so somewhat below the true 0.3.
stopifnot(all(vh[upper.tri(vh)] > .12 & vh[upper.tri(vh)] < .4), all(abs(diag(vh) - 1) < 1e-12))
learned <- interp_mash_learned(sim$zh[apply(abs(sim$zh), 1L, max, na.rm = TRUE) > 4, ], vh, mctx)
stopifnot("learned_all" %in% names(learned), max(diag(learned$learned_all)) == 1, isSymmetric(unname(learned$learned_all)))
mm <- interp_mash_model(mctx, vh, c(canon, learned))
ll <- interp_mash_loglik(sim$zh, mm)
# The null column is the plain multivariate normal density of the tested contexts.
one_row <- which(rowSums(is.finite(sim$zh)) == 1L)[[1L]]
stopifnot(dim(ll) == c(nrow(sim$zh), length(mm$pi)), all(is.finite(ll)),
  abs(ll[one_row, 1L] - stats::dnorm(sim$zh[one_row, is.finite(sim$zh[one_row, ])], log = TRUE)) < 1e-10)
mm$pi <- interp_mash_fit_pi(ll)
by_class <- tapply(mm$pi, interp_mash_class(interp_mash_components(mm)$pattern), sum)
stopifnot(abs(sum(mm$pi) - 1) < 1e-8, abs(by_class[["no effect"]] - .7) < .1, by_class[["shared"]] + by_class[["learned"]] > .12,
  by_class[["one context"]] > .04, by_class[["independent"]] < .05)
# Same effect in two contexts: confident there, and predicted in the untested third.
shared <- interp_mash_posterior(c(5, NA, 5.5), mm)
stopifnot(all(shared$contexts$lfsr < .01), !shared$contexts$tested[[2L]], shared$contexts$posterior_mean[[2L]] > 3, shared$weights$class[[1L]] %in% c("shared", "learned"),
  shared$log_bf > 5, all(shared$contexts$posterior_mean[c(1L, 3L)] < c(5, 5.5) + .5), identical(shared$contexts$estimate, c(5, NA, 5.5)))
# One context only: the others stay at zero.
single <- interp_mash_posterior(c(.3, 6, -.5), mm)
stopifnot(single$contexts$lfsr[[2L]] < .01, all(single$contexts$lfsr[c(1L, 3L)] > .5), single$weights$class[[1L]] == "one context", single$top_single$context == "B")
# Nothing: no context is called, the Bayes factor favours no effect.
quiet <- interp_mash_posterior(c(1.2, -.8, .5), mm)
stopifnot(all(quiet$contexts$lfsr > .5), quiet$log_bf < 0, all(is.na(interp_mash_posterior(c(NA, NA, NA), mm)$contexts$lfsr)))
# A left-out context is predicted from the others, with an interval that covers about 95% of fresh pairs.
fresh <- simulate(600L)$zh; fresh <- fresh[is.finite(fresh[, 2L]) & rowSums(is.finite(fresh)) > 1L, ]
pred <- t(vapply(seq_len(nrow(fresh)), function(i) interp_mash_predict_left_out(fresh[i, ], mm, 2L), numeric(2)))
stopifnot(stats::cor(pred[, 1L], fresh[, 2L]) > .5, abs(mean(abs(fresh[, 2L] - pred[, 1L]) <= 1.96 * pred[, 2L]) - .95) < .04)
model_file <- tempfile(fileext = ".json"); mm$meta <- list(qtl_type = "eQTL", trained_at = "then", strong_pairs_fitted = 5L)
interp_tissue_model_write(mm, model_file); back <- interp_tissue_model_read(model_file)
stopifnot(identical(back$contexts, mctx), back$scale == "z", all.equal(unname(back$V), unname(mm$V)), all.equal(back$pi, mm$pi), identical(names(back$U), names(mm$U)),
  all.equal(unname(back$U$learned_all), unname(mm$U$learned_all)), back$meta$trained_at == "then",
  all.equal(interp_mash_posterior(c(5, NA, 5.5), back)$contexts$lfsr, shared$contexts$lfsr))
writeLines("{\"format\": \"other\"}", model_file); expect_error(interp_tissue_model_read(model_file), "Not a tissue sharing model")
# Slope scale: contexts share the effect size, and each estimate counts by its own standard error.
sim_slope <- function(n) {
  se <- cbind(stats::runif(n, .03, .06), stats::runif(n, .1, .2), stats::runif(n, .05, .1))
  kind <- sample(c("null", "shared", "onlyB"), n, TRUE, c(.6, .3, .1)); b <- matrix(0, n, 3L)
  b[kind == "shared", ] <- stats::rnorm(sum(kind == "shared"), 0, .5); b[kind == "onlyB", 2L] <- stats::rnorm(sum(kind == "onlyB"), 0, .8)
  bh <- b + (matrix(stats::rnorm(n * 3L), n, 3L) %*% chol(mV)) * se; bh[sample(length(bh), round(length(bh) * .15))] <- NA
  keep <- rowSums(is.finite(bh)) > 0L
  list(b = bh[keep, , drop = FALSE], se = se[keep, , drop = FALSE])
}
ss <- sim_slope(900L)
big <- apply(abs(ss$b / ss$se), 1L, max, na.rm = TRUE) > 4
sm <- interp_mash_model(mctx, vh, c(canon, interp_mash_learned(ss$b[big, ], vh, mctx, SE = ss$se[big, ])), scale = "slope")
stopifnot(identical(sm$grid, interp_mash_grid("slope")), max(sm$grid) < 4, min(sm$grid) < .05)
expect_error(interp_mash_loglik(ss$b, sm), "is.null")
sl <- interp_mash_loglik(ss$b, sm, ss$se)
# One tested context: the null density is a normal with that context's own standard error.
one_row <- which(rowSums(is.finite(ss$b)) == 1L)[[1L]]; col <- which(is.finite(ss$b[one_row, ]))
stopifnot(all(is.finite(sl)), abs(sl[one_row, 1L] - stats::dnorm(ss$b[one_row, col], 0, ss$se[one_row, col], log = TRUE)) < 1e-10)
sm$pi <- interp_mash_fit_pi(sl)
s_class <- tapply(sm$pi, interp_mash_class(interp_mash_components(sm)$pattern), sum)
stopifnot(abs(s_class[["no effect"]] - .6) < .12, s_class[["one context"]] > .03)
# A precise estimate in A and a noisy one in C: B, untested, is predicted on the slope scale with A's size.
sp <- interp_mash_posterior(c(.5, NA, .7), sm, c(.04, NA, .3))
stopifnot(abs(sp$contexts$posterior_mean[[2L]] - .5) < .1, sp$contexts$lfsr[[2L]] < .05, !sp$contexts$tested[[2L]],
  abs(sp$contexts$posterior_mean[[3L]] - .5) < .1, sp$contexts$posterior_sd[[1L]] < .05,
  # Without a standard error a context does not count as tested.
  !interp_mash_posterior(c(.5, .4, .7), sm, c(.04, NA, .3))$contexts$tested[[2L]])
fresh <- sim_slope(500L); ok <- is.finite(fresh$b[, 1L]) & rowSums(is.finite(fresh$b)) > 1L
pred <- t(vapply(which(ok), function(i) interp_mash_predict_left_out(fresh$b[i, ], sm, 1L, fresh$se[i, ]), numeric(2)))
stopifnot(stats::cor(pred[, 1L], fresh$b[ok, 1L]) > .5, abs(mean(abs(fresh$b[ok, 1L] - pred[, 1L]) <= 1.96 * pred[, 2L]) - .95) < .05)
interp_tissue_model_write(sm, model_file); stopifnot(interp_tissue_model_read(model_file)$scale == "slope")
slope_ts <- interp_tissue_sharing(rbind(pair("eQTL", "A", "ENSG1.1", .5, .04, 1e-30, TRUE), pair("eQTL", "C", "ENSG1.1", .7, .3, .02, FALSE)), models = list(eQTL = sm))
stopifnot(identical(slope_ts$contexts$call, c("effect", "predicted effect", "effect")), abs(slope_ts$contexts$posterior_slope[[2L]] - .5) < .1,
  is.na(slope_ts$contexts$posterior_z[[2L]]), slope_ts$models$eQTL$scale == "slope",
  # C alone is not significant, but with A it is called: the small context borrows strength.
  slope_ts$phenotypes$summary == "effect in all 2 tested contexts, similar size; predicted in B")
# Applied to a variant's pairs: one row per context of the model, calls from the lfsr, slopes back on the slope scale.
tp <- rbind(pair("eQTL", "A", "ENSG1.1", .5, .1, 1e-7, TRUE), pair("eQTL", "C", "ENSG1.1", .55, .1, 1e-8, TRUE),
  pair("eQTL", "A", "ENSG2.1", .03, .1, .7, FALSE, "ENSG2.1"), pair("eQTL", "B", "ENSG2.1", .6, .1, 1e-9, TRUE, "ENSG2.1"), pair("eQTL", "C", "ENSG2.1", -.05, .1, .6, FALSE, "ENSG2.1"),
  pair("eQTL", "A", "ENSG3.1", .1, .1, .3, FALSE, "ENSG3.1"), pair("sQTL", "A", "chr1:10:90:clu_1_+:ENSG1.1", 2, .1, 1e-50, TRUE))
ts <- interp_tissue_sharing(tp, models = list(eQTL = mm))
first <- ts$contexts[ts$contexts$phenotype == "ENSG1.1", ]
stopifnot(ts$available, nrow(ts$phenotypes) == 2L, nrow(ts$contexts) == 6L, identical(first$context, mctx), identical(first$call, c("effect", "predicted effect", "effect")),
  is.na(first$posterior_slope[[2L]]), abs(first$posterior_slope[[1L]] - first$posterior_z[[1L]] * .1) < 1e-12, first$z[[1L]] == 5,
  ts$phenotypes$summary[ts$phenotypes$phenotype == "ENSG1.1"] == "effect in all 2 tested contexts, similar size; predicted in B",
  ts$phenotypes$summary[ts$phenotypes$phenotype == "ENSG2.1"] == "effect in 1 of 3 tested contexts (B)",
  ts$phenotypes$leading_pattern[ts$phenotypes$phenotype == "ENSG2.1"] == "one context (B)", ts$phenotypes$no_clear_effect_in[ts$phenotypes$phenotype == "ENSG2.1"] == "A, C",
  # Not significant anywhere, or no model for that QTL type: left out.
  !"ENSG3.1" %in% ts$phenotypes$phenotype, !"sQTL" %in% ts$phenotypes$qtl_type, ts$models$eQTL$trained_at == "then")
none_model <- interp_tissue_sharing(tp, none)
stopifnot(!none_model$available, !nrow(none_model$phenotypes), !nrow(interp_tissue_sharing(tp[0, ], models = list(eQTL = mm))$phenotypes))
# The Tissue line carries the model's reading next to the comparison of the contexts.
tissue_line <- interp_evidence_summary(missense, list(topmed = list(available = TRUE, pairs = tp, contexts = data.frame()), tissue = ts))
stopifnot(grepl("model: effect in 1 of 3 tested contexts (B)", tissue_line$evidence[tissue_line$question == "Tissue"], fixed = TRUE),
  grepl("cross-context model", tissue_line$sources[tissue_line$question == "Tissue"], fixed = TRUE))

# ---- Sorted-file lookup and tabix-backed sources ----------------------------
work <- tempfile("regshiny_interp_"); dir.create(work, mode = "0700")
if (nzchar(Sys.which("look"))) {
  csv <- file.path(work, "sqtl.csv")
  body <- c("chr1:100:G:A:rs1,0.1,T1,GENEA,1,0.5,0.2,5,-30,2622", "chr1:100:G:C:rs1,0.1,T1,GENEA,1,0.5,0.2,5,-30,2622",
    "chr1:1000:G:A:rs2,0.1,T9,GENEB,0,0.5,0.2,5,-30,2622", "chr1:100:G:A:rs1,0.1,T2,FAR,0,0.01,-0.1,-4,-4.2,2622")
  writeLines(c("SNP,EAF,Tx,Tx_Gene,Is_Cis,RSq,Fx,T,log10P,N", sort(body, method = "radix")), csv)
  srcs <- none; srcs$look <- unname(Sys.which("look")); srcs$fhs_sqtl <- c(discovery = csv, replication = csv); srcs$fhs_eqtl <- csv
  f <- interp_fhs_qtl(v, srcs)
  stopifnot(f$sqtl_available, f$eqtl_available, nrow(f$sqtl) == 4L, nrow(f$eqtl) == 2L, f$eqtl$gene[[1L]] == "GENEA", f$eqtl$cis[[1L]],
    !f$eqtl$cis[[2L]], f$eqtl$beta_alt[[1L]] == .2, f$eqtl$alt_frequency[[1L]] == .1, f$eqtl$n[[1L]] == 2622L,
    !nrow(interp_fhs_qtl(interp_variant("chr1", 10L, "G", "A"), srcs)$eqtl), !nrow(interp_fhs_qtl(interp_variant("chr1", 100L, "G", "T"), srcs)$eqtl))
  # Reactome mapping through its build script, and MSigDB gene sets.
  raw <- file.path(work, "Ensembl2Reactome_All_Levels.txt")
  writeLines(c("ENSG00000000002\tR-HSA-9\turl\tSignal Transduction\tIEA\tHomo sapiens", "ENSMUSG00000000001\tR-MMU-1\turl\tMetabolism\tIEA\tMus musculus",
    "ENSG00000000001\tR-HSA-2\turl\tSphingolipid metabolism\tTAS\tHomo sapiens", "ENSG00000000001\tR-HSA-1\turl\tMetabolism\tTAS\tHomo sapiens",
    "ENSG00000000001\tR-HSA-1\turl\tMetabolism\tIEA\tHomo sapiens", "P12345\tR-HSA-1\turl\tMetabolism\tTAS\tHomo sapiens"), raw)
  reactome <- file.path(work, "pathway", "reactome.tsv")
  built <- system2("bash", c("scripts/build_pathway_index.sh", shQuote(raw), shQuote(reactome)), stdout = TRUE, stderr = TRUE)
  stopifnot(file.exists(reactome), length(readLines(reactome)) == 4L, grepl("genes: 2", built[[length(built)]], fixed = TRUE),
    grepl("up to date", system2("bash", c("scripts/build_pathway_index.sh", shQuote(raw), shQuote(reactome)), stdout = TRUE), fixed = TRUE))
  msig <- file.path(work, "msigdb"); dir.create(msig)
  writeLines(c("HALLMARK_ONE\turl\tGENEA\tGENEB\tGENEC", "HALLMARK_TWO\turl\tGENEB"), file.path(msig, "h.all.v2099.1.Hs.symbols.gmt"))
  writeLines(c("WP_SMALL\turl\tGENEA\tX1", "REACTOME_DROPPED\turl\tGENEA", "KEGG_BIG\turl\tGENEA\tX1\tX2\tX3"), file.path(msig, "c2.cp.v2099.1.Hs.symbols.gmt"))
  srcs$reactome <- reactome; srcs$msigdb <- msig
  linked <- data.frame(gene_id = c("ENSG00000000001", "ENSG00000000003", "ENSG00000000002"), gene = c("GENEA", "", "GENEAB"),
    relation = c("overlaps", "eQTL target (TOPMed)", "sQTL target (TOPMed)"), stringsAsFactors = FALSE)
  lp <- interp_local_pathways(linked, srcs)
  stopifnot(lp$available, lp$reactome_available, lp$msigdb_available, nrow(lp$reactome) == 3L,
    # Top-level pathways first; the two evidence codes of one pathway are merged.
    identical(lp$reactome$pathway[lp$reactome$gene == "GENEA"], c("Metabolism", "Sphingolipid metabolism")), lp$reactome$top_level[[1L]],
    !lp$reactome$top_level[[2L]], lp$reactome$evidence[[1L]] == "IEA, TAS", lp$reactome$gene[[3L]] == "GENEAB",
    # Whole-symbol match only, Reactome sets left to the mapping, smaller sets first within a collection.
    identical(lp$gene_sets$gene_set, c("HALLMARK_ONE", "WP_SMALL", "KEGG_BIG")), identical(lp$gene_sets$genes_in_set, c(3L, 2L, 4L)),
    identical(lp$gene_sets$collection, c("Hallmark", "Canonical pathways", "Canonical pathways")), grepl("MSigDB 2099.1", lp$source[[2L]], fixed = TRUE))
  st <- interp_local_status(srcs)
  stopifnot(st$available[st$source == "Reactome gene-to-pathway mapping"], st$available[st$source == "MSigDB gene sets"])
  sp <- interp_evidence_summary(interp_empty_consequences(), list(pathways = lp))
  stopifnot(level(sp, "Pathway") == "Context", grepl("GENEA (overlaps): Metabolism", sp$finding[sp$question == "Pathway"], fixed = TRUE),
    grepl("no pathway annotation: ENSG00000000003", sp$evidence[sp$question == "Pathway"], fixed = TRUE))
} else if (native) stop("--native needs look")
tabix <- interp_local_sources()$tabix; bgzip <- file.path(dirname(tabix), "bgzip")
if (nzchar(tabix) && file.exists(bgzip)) {
  index <- function(lines, path, args) {
    plain <- sub("\\.gz$", "", path); writeLines(lines, plain)
    stopifnot(system2(bgzip, c("-f", shQuote(plain))) == 0L, system2(tabix, c(args, shQuote(path))) == 0L)
  }
  cadd <- file.path(work, "cadd.tsv.gz")
  index(c("## CADD synthetic", "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED", "1\t100\tG\tA\t2.5\t22.1", "1\t100\tG\tC\t0.1\t3.0", "1\t100\tGT\tG\t1.0\t12.0"),
    cadd, c("-s", "1", "-b", "2", "-e", "2"))
  srcs <- none; srcs$tabix <- tabix; srcs$cadd_snv <- cadd; srcs$cadd_indel <- cadd
  c1 <- interp_cadd(v, srcs)
  stopifnot(c1$available, c1$found, c1$phred == 22.1, c1$raw == 2.5, !interp_cadd(interp_variant("chr1", 100L, "G", "T"), srcs)$found,
    interp_cadd(interp_variant("chr1", 100L, "GT", "G"), srcs)$phred == 12, !interp_cadd(interp_variant("chrM", 100L, "G", "A"), srcs)$available)
  cache <- file.path(work, "vep"); dir.create(file.path(cache, "1"), recursive = TRUE)
  writeLines(c("species\thomo_sapiens", "variation_cols\tchr,variation_name,failed,somatic,start,end,allele_string,strand,clin_sig,phenotype_or_disease,clin_sig_allele,pubmed,var_synonyms,AF,EUR,gnomADg"),
    file.path(cache, "info.txt"))
  index(c("1\tCM1\t.\t1\t100\t.\tHGMD_MUTATION\t.\t.\t1\t.\t.\t.\t.\t.\t.",
    "1\trs1\t.\t.\t100\t.\tG/A/C\t.\tpathogenic,benign\t1\tA:pathogenic/likely_pathogenic;C:benign\t11,12\t.\tA:0.01,C:0.2\tA:0.03\tA:0.025,C:0.1",
    "1\trs9\t.\t.\t101\t103\tTCT/-\t.\t.\t.\t.\t.\t.\t.\t.\t-:0.4"), file.path(cache, "1", "all_vars.gz"), c("-s", "1", "-b", "5", "-e", "5"))
  srcs$vep_cache <- cache
  k <- interp_known_variant(v, srcs)
  stopifnot(k$found, k$rs_id == "rs1", k$clinvar_significance == "pathogenic, likely_pathogenic", k$pubmed_n == 2L,
    identical(k$frequencies$population, c("AF", "EUR", "gnomADg")), identical(k$frequencies$allele_frequency, c(.01, .03, .025)))
  kc <- interp_known_variant(interp_variant("chr1", 100L, "G", "C"), srcs)
  stopifnot(kc$clinvar_significance == "benign", kc$frequencies$allele_frequency[[1L]] == .2,
    !interp_known_variant(interp_variant("chr1", 100L, "G", "T"), srcs)$found,
    interp_known_variant(interp_variant("chr1", 100L, "ATCT", "A"), srcs)$rs_id == "rs9")
  # TOPMed layout: two contexts, eQTL and sQTL, thresholds and credible sets.
  root <- file.path(work, "topmed")
  for (d in c("cis-eqtl", "cis-sqtl2")) for (s in c("nominal", "top-hits", "susie")) dir.create(file.path(root, d, s), recursive = TRUE)
  header <- "#chrom\tstart\tend\tphenotype_id\tvariant_id\ttss_distance\taf\tma_samples\tma_count\tpval_nominal\tslope\tslope_se\tnegative_log10_p"
  row <- function(pheno, variant, af, mac, p, slope, se) paste("chr1", 99, 100, pheno, variant, 500, af, mac, mac, p, slope, se, -log10(p), sep = "\t")
  args <- c("-s", "1", "-b", "3", "-e", "3")
  index(c(header, row("ENSG1.1", "chr1_100_G_A", .1, 200, 1e-12, .5, .07), row("ENSG2.1", "chr1_100_G_A", .1, 200, .2, .05, .07),
    row("ENSG1.1", "chr1_100_G_C", .3, 600, 1e-30, 1, .07)), file.path(root, "cis-eqtl", "nominal", "Lung.maf001.tsv.gz"), args)
  index(c(header, row("ENSG1.1", "chr1_100_G_A", .9, 100, .03, .2, .1)), file.path(root, "cis-eqtl", "nominal", "PBMC.maf001.tsv.gz"), args)
  index(c(header, row("chr1:10:90:clu_1_+:ENSG1.1", "chr1_100_G_A", .1, 200, 1e-40, -1.5, .1), row("chr1:10:95:clu_1_+:ENSG1.1", "chr1_100_G_A", .1, 200, .6, .01, .1)),
    file.path(root, "cis-sqtl2", "nominal", "Lung.maf001.tsv.gz"), args)
  index(c(header, row("ENSG1.1", "chr1_100_G_A", .004, 40, 1e-9, .9, .1)), file.path(root, "cis-eqtl", "nominal", "Whole_blood.maf0001.tsv.gz"), args)
  index(c(header, row("ENSG7.1", "chr1_500_G_A", .2, 40, 1e-9, .9, .1)), file.path(root, "cis-eqtl", "nominal", "Whole_blood.maf001.tsv.gz"), args)
  hits <- function(path, ids, q, threshold) {
    con <- gzfile(path, "w"); on.exit(close(con))
    writeLines(c("phenotype_id\tqval\tpval_nominal_threshold", paste(ids, q, threshold, sep = "\t")), con)
  }
  hits(file.path(root, "cis-eqtl", "top-hits", "Lung.75.maf001.permutations.txt.gz"), c("ENSG1.1", "ENSG2.1"), c(.001, .5), c(1e-5, .5))
  hits(file.path(root, "cis-eqtl", "top-hits", "PBMC.30.maf001.permutations.txt.gz"), "ENSG1.1", .01, 1e-4)
  hits(file.path(root, "cis-sqtl2", "top-hits", "Lung.maf001.permutations.txt.gz"), "chr1:10:90:clu_1_+:ENSG1.1", .001, 1e-6)
  hits(file.path(root, "cis-eqtl", "top-hits", "Whole_blood.100.maf0001.permutations.txt.gz"), "ENSG1.1", .001, 1e-5)
  writeLines(c("phenotype_id\tvariant_id\tpip\taf\tcs_id", "ENSG1.1\tchr1_100_G_A\t0.83\t0.1\t1", "ENSG1.1\tchr1_100_G_AA\t0.5\t0.1\t1"),
    file.path(root, "cis-eqtl", "susie", "Lung.maf001.cs.txt"))
  writeLines(c("phenotype_id\tvariant_id\tpip\taf\tcs_id", "chr1:10:90:clu_1_+:ENSG1.1\tchr1_100_G_A\t0.99\t0.1\t2"),
    file.path(root, "cis-sqtl2", "susie", "Lung.maf001.cs-per-intron.txt"))
  srcs$topmed_qtl <- root
  t <- interp_topmed_qtl(v, srcs)
  ctx <- function(type, context) t$contexts[t$contexts$qtl_type == type & t$contexts$context == context, ]
  pr <- function(type, context, pheno) t$pairs[t$pairs$qtl_type == type & t$pairs$context == context & t$pairs$phenotype == pheno, ]
  lung <- pr("eQTL", "Lung", "ENSG1.1")
  stopifnot(t$available, nrow(t$contexts) == 12L, ctx("eQTL", "Lung")$pairs_tested == 2L, ctx("eQTL", "Lung")$significant_pairs == 1L,
    ctx("eQTL", "Monocytes")$tested == FALSE, ctx("sQTL", "PBMC")$tested == FALSE,
    lung$significant, lung$pip == .83, lung$credible_set == "1", lung$n_estimated == 1000L, lung$p_threshold == 1e-5, lung$slope_alt == .5,
    # The ENSG2 gene failed FDR, so it has no threshold and cannot be significant.
    !pr("eQTL", "Lung", "ENSG2.1")$significant, is.na(pr("eQTL", "Lung", "ENSG2.1")$p_threshold),
    # PBMC: same gene, not significant, kept for the comparison; minor allele is REF there.
    !pr("eQTL", "PBMC", "ENSG1.1")$significant, pr("eQTL", "PBMC", "ENSG1.1")$n_estimated == 500L,
    # Whole blood falls back to the 0.1% run when the 1% run has no row for this variant.
    pr("eQTL", "Whole blood", "ENSG1.1")$significant, pr("eQTL", "Whole blood", "ENSG1.1")$alt_frequency == .004,
    pr("sQTL", "Lung", "chr1:10:90:clu_1_+:ENSG1.1")$pip == .99, pr("sQTL", "Lung", "chr1:10:90:clu_1_+:ENSG1.1")$significant,
    # The other allele at this position never leaks in.
    !any(t$pairs$slope_alt == 1))
  hh <- interp_context_heterogeneity(t$pairs)
  stopifnot(hh$contexts_tested[hh$qtl_type == "eQTL"] == 3L, hh$contexts_significant[hh$qtl_type == "eQTL"] == "Lung, Whole blood")
  stopifnot(!nrow(interp_topmed_qtl(interp_variant("chr1", 100L, "G", "T"), srcs)$pairs))
  # GWAS Catalog through its build script: columns found by name, multi-variant rows split, unmapped rows dropped.
  catalog_cols <- c("PUBMEDID", "FIRST AUTHOR", "DATE", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "CHR_ID", "CHR_POS", "MAPPED_GENE",
    "STRONGEST SNP-RISK ALLELE", "SNPS", "RISK ALLELE FREQUENCY", "P-VALUE", "PVALUE_MLOG", "OR or BETA", "95% CI (TEXT)", "MAPPED_TRAIT", "STUDY ACCESSION")
  catalog_row <- function(chr, pos, allele, snp, trait, p, mlog, effect = "") paste("1", "Author A", "2020-01-01", trait, "1,000 cases", chr, pos, "GENEA",
    allele, snp, "0.1", p, mlog, effect, "[1.1-1.3]", tolower(trait), "GCST1", sep = "\t")
  raw <- file.path(work, "gwas.tsv")
  writeLines(c(paste(catalog_cols, collapse = "\t"), catalog_row("1", "100", "rs1-A", "rs1", "Weak trait", "3E-6", "5.52"),
    catalog_row("1", "100", "rs1-G", "rs1", "Strong trait", "1E-400", "400", "1.2"), catalog_row("1", "100", "rs1-?", "rs1", "Middle trait", "2E-9", "8.7"),
    catalog_row("1;2", "100;300", "rs1-A; rs5-T", "rs1; rs5", "Pair trait", "5E-8", "7.3"), catalog_row("", "", "rs7-A", "rs7", "Unmapped", "1E-9", "9"),
    catalog_row("X", "50", "rs8-C", "rs8", "Other place", "1E-9", "9")), raw)
  catalog <- file.path(work, "disease", "gwas_catalog_e999_r2099-01-01.b38.tsv.gz")
  built <- system2("bash", c("scripts/build_gwas_index.sh", shQuote(raw), shQuote(catalog)), stdout = TRUE, stderr = TRUE,
    env = c(paste0("BGZIP=", bgzip), paste0("TABIX=", tabix)))
  stopifnot(file.exists(catalog), file.exists(paste0(catalog, ".tbi")), grepl("rows: 6", built[[length(built)]], fixed = TRUE))
  srcs$gwas_catalog <- catalog
  g <- interp_gwas_catalog(v, srcs)
  stopifnot(g$available, g$found, g$release == "e999_r2099-01-01", nrow(g$associations) == 4L,
    # Strongest first by the stored -log10 p; a p-value below double precision is kept as 0.
    identical(g$associations$trait, c("Strong trait", "Middle trait", "Pair trait", "Weak trait")), g$associations$p_value[[1L]] == 0,
    identical(g$associations$risk_allele_is, c("REF", "not reported", "ALT", "ALT")), g$associations$effect[[1L]] == 1.2,
    g$associations$rs_id[[3L]] == "rs1", g$associations$mapped_gene[[1L]] == "GENEA", g$associations$pubmed[[1L]] == "1")
  far <- interp_gwas_catalog(interp_variant("chr1", 101L, "G", "A"), srcs)
  stopifnot(far$available, !far$found, interp_gwas_catalog(interp_variant("chr2", 300L, "C", "T"), srcs)$associations$rs_id == "rs5",
    interp_gwas_catalog(interp_variant("chrX", 50L, "A", "C"), srcs)$associations$risk_allele_is == "ALT")
  # ClinVar VCF: this exact allele only; other alleles at the position are counted, not shown.
  clinvar <- file.path(work, "disease", "clinvar.vcf.gz")
  index(c("##fileformat=VCFv4.1", "##fileDate=2099-02-03", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t100\t11\tG\tA\t.\t.\tCLNDN=Disease_one|not_provided;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic/Likely_pathogenic;GENEINFO=GENEA:1;MC=SO:0001583|missense_variant;RS=1",
    "1\t100\t12\tG\tC\t.\t.\tCLNDN=Disease_two;CLNREVSTAT=no_assertion_criteria_provided;CLNSIG=Benign",
    "1\t100\t13\tGT\tG\t.\t.\tCLNSIG=Uncertain_significance"), clinvar, c("-p", "vcf"))
  srcs$clinvar_vcf <- clinvar
  cv <- interp_clinvar_local(v, srcs)
  stopifnot(cv$available, cv$found, nrow(cv$records) == 1L, cv$records$variation_id == "11", cv$records$clinical_significance == "Pathogenic/Likely pathogenic",
    cv$records$conditions == "Disease one; not provided", cv$records$review_status == "criteria provided, single submitter", cv$records$gene == "GENEA",
    cv$records$molecular_consequence == "missense variant", cv$records$rs_id == "rs1", cv$other_alleles == 2L, cv$file_date == "2099-02-03",
    cv$source == "ClinVar VCF 2099-02-03 (local)")
  none_here <- interp_clinvar_local(interp_variant("chr1", 100L, "G", "T"), srcs)
  stopifnot(none_here$available, !none_here$found, none_here$other_alleles == 3L, interp_clinvar_local(interp_variant("chr1", 100L, "GT", "G"), srcs)$found)
  # Disease level from local sources alone.
  disease <- function(local) { s <- interp_evidence_summary(missense, local); s[s$question == "Disease and traits", ] }
  d1 <- disease(list(clinvar = cv, gwas_catalog = g))
  stopifnot(d1$level == "Strong", grepl("ClinVar: pathogenic, likely pathogenic for Disease one", d1$evidence, fixed = TRUE),
    grepl("GWAS Catalog: 4 association(s) reported at this position, 3 at p <= 5e-8, 4 distinct trait(s)", d1$evidence, fixed = TRUE),
    grepl("not consulted: public query not run", d1$evidence, fixed = TRUE), grepl("ClinVar VCF 2099-02-03 (local); GWAS Catalog e999_r2099-01-01", d1$sources, fixed = TRUE))
  # A published association at the variant is Moderate; below genome-wide significance, Limited.
  stopifnot(disease(list(clinvar = none_here, gwas_catalog = g))$level == "Moderate")
  weak_only <- g; weak_only$associations <- g$associations[g$associations$trait == "Weak trait", ]
  stopifnot(disease(list(clinvar = none_here, gwas_catalog = weak_only))$level == "Limited")
  # Both local sources consulted and empty: None found, naming what was read. One of them missing: Not assessed.
  d0 <- disease(list(clinvar = none_here, gwas_catalog = far))
  stopifnot(d0$level == "None found", grepl("ClinVar VCF 2099-02-03 (local), GWAS Catalog e999_r2099-01-01", d0$finding, fixed = TRUE),
    disease(list(clinvar = none_here))$level == "Not assessed", disease(list(gwas_catalog = far))$level == "Not assessed",
    grepl("Only GWAS Catalog", disease(list(gwas_catalog = far))$finding, fixed = TRUE))
  st <- interp_local_status(srcs)
  stopifnot(st$available[st$source == "CADD v1.7"], st$available[st$source == "TOPMed multi-tissue cis-QTL"], !st$available[st$source == "SpliceAI 1.3.1"],
    st$available[st$source == "GWAS Catalog"], st$available[st$source == "ClinVar VCF"])
  # LD with nearby reported variants, from a small VCF: twelve samples, the target's ALT dosage 0,0,0,0,1,1,1,1,2,2,2,2.
  stopifnot(identical(interp_gt_dosage(c("0/0", "0|1", "1/1", "./.", "1", "0", "1/2", "2|2", "./1"), 1L), c(0, 1, 2, NA, 1, 0, 1, 0, NA)),
    identical(interp_gt_dosage(c("1/2", "2|2", "0/0"), 2L), c(1, 2, 0)))
  bcftools <- unname(Sys.which("bcftools"))
  if (nzchar(bcftools)) {
    gt <- function(d) c("0/0", "0|1", "1/1")[d + 1L]
    own <- c(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2); vcf_row <- function(pos, id, ref, alt, calls) paste(c("chr1", pos, id, ref, alt, ".", "PASS", ".", "GT", calls), collapse = "\t")
    partial <- gt(own); partial[[3L]] <- "./."
    third <- rep("0/0", 12L); third[own == 1] <- "0/2"; third[own == 2] <- "2/2"; third[[1L]] <- "0/1"
    ld_vcf <- file.path(work, "disease", "cohort.vcf.gz")
    index(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste0("S", 1:12)), collapse = "\t"),
      vcf_row(100, "rs1", "G", "A", gt(own)), vcf_row(150, "rs150", "C", "T", gt(own)), vcf_row(180, "rs180", "A", "G", gt(2 - own)),
      vcf_row(220, "rs220", "T", "C", gt(rep(0:2, 4L))), vcf_row(260, "rs260", "G", "A", partial), vcf_row(260, ".", "G", "T", gt(rep(c(0, 1), 6L))),
      vcf_row(300, "rs300", "A", "C,T", third), vcf_row(350, "rs350", "C", "G", rep("0/0", 12L))), ld_vcf, c("-p", "vcf"))
    raw <- file.path(work, "gwas_ld.tsv")
    writeLines(c(paste(catalog_cols, collapse = "\t"), catalog_row("1", "100", "rs1-A", "rs1", "Own trait", "1E-20", "20"),
      catalog_row("1", "150", "rs150-T", "rs150", "Level of protein X in blood", "1E-30", "30"), catalog_row("1", "180", "rs180-G", "rs180", "Opposite trait", "2E-9", "8.7"),
      catalog_row("1", "220", "rs220-C", "rs220", "Unlinked trait", "1E-12", "12"), catalog_row("1", "260", "rs260-?", "rs260", "Partial trait", "3E-6", "5.5"),
      catalog_row("1", "300", "rs300-T", "rs300", "Third allele trait", "1E-9", "9"), catalog_row("1", "350", "rs350-G", "rs350", "Monomorphic", "1E-9", "9"),
      catalog_row("1", "400", "rs400-A", "rs400", "Absent from the VCF", "1E-9", "9"), catalog_row("1", "90000", "rs9-A", "rs9", "Far away", "1E-9", "9")), raw)
    ld_catalog <- file.path(work, "disease", "gwas_catalog_e998_r2099-01-02.b38.tsv.gz")
    stopifnot(!is.null(system2("bash", c("scripts/build_gwas_index.sh", shQuote(raw), shQuote(ld_catalog)), stdout = TRUE, stderr = TRUE,
      env = c(paste0("BGZIP=", bgzip), paste0("TABIX=", tabix)))), file.exists(ld_catalog))
    ld_src <- srcs; ld_src$gwas_catalog <- ld_catalog
    ld <- interp_gwas_ld(v, ld_vcf, ld_src, window = 1000L)
    row_of <- function(id, alt = NULL) ld$ld[ld$ld$rs_id == id & (is.null(alt) | ld$ld$alt == (if (is.null(alt)) "" else alt)), ]
    stopifnot(ld$available, ld$found, ld$samples == 12L, ld$positions_in_window == 7L, ld$positions_examined == 7L, !ld$truncated, ld$release == "e998_r2099-01-02",
      # The variant's own catalog row is not an LD partner; unlinked, monomorphic and absent variants give no row.
      !any(c("rs1", "rs220", "rs350", "rs400", "rs9") %in% ld$ld$rs_id),
      abs(row_of("rs150")$r2 - 1) < 1e-12, row_of("rs150")$r > 0, row_of("rs150")$risk_allele_travels_with == "ALT", row_of("rs150")$distance == 50L, row_of("rs150")$n == 12L,
      # Opposite phase: r2 is 1 and the reported risk allele is found with this variant's REF.
      abs(row_of("rs180")$r2 - 1) < 1e-12, row_of("rs180")$r < 0, row_of("rs180")$risk_allele_travels_with == "REF",
      # A missing call drops that sample; the unnamed second record at the position is not used for rs260.
      nrow(row_of("rs260")) == 1L, row_of("rs260")$alt == "A", row_of("rs260")$n == 11L, row_of("rs260")$risk_allele_travels_with == "",
      # Multi-allelic record: each ALT has its own dosage; the T allele follows the target apart from one sample.
      row_of("rs300", "T")$r2 > .8, row_of("rs300", "T")$risk_allele_travels_with == "ALT", !nrow(row_of("rs300", "C")),
      # Strongest LD first, then the smaller p-value.
      ld$ld$rs_id[[1L]] == "rs150", all(diff(round(ld$ld$r2, 6)) <= 0))
    near <- interp_gwas_ld(v, ld_vcf, ld_src, window = 1000L, max_positions = 2L)
    stopifnot(near$truncated, near$positions_examined == 2L, near$positions_in_window == 7L, identical(sort(unique(near$ld$rs_id)), c("rs150", "rs180")))
    expect_error(interp_gwas_ld(interp_variant("chr1", 150L, "C", "G"), ld_vcf, ld_src, window = 1000L), "own record was not found")
    expect_error(interp_gwas_ld(interp_variant("chr1", 350L, "C", "G"), ld_vcf, ld_src, window = 1000L), "does not vary")
    expect_error(interp_gwas_ld(v, "/missing.vcf.gz", ld_src), "VCF for this chromosome is not readable")
    expect_error(interp_gwas_ld(v, ld_vcf, none), "GWAS Catalog copy is not readable")
    alone <- interp_gwas_ld(interp_variant("chr1", 100L, "G", "A"), ld_vcf, ld_src, window = 10L)
    stopifnot(alone$available, !alone$found, alone$positions_examined == 0L)
    # Disease: LD with a genome-wide significant reported variant is an association (Moderate), named with its traits.
    dl <- interp_evidence_summary(missense, list(clinvar = none_here, gwas_catalog = far), ld = ld)
    dline <- dl[dl$question == "Disease and traits", ]
    stopifnot(dline$level == "Moderate", grepl("r2 >= 0.8 with 4 other reported variant(s) (4 association(s), 3 at p <= 5e-8)", dline$evidence, fixed = TRUE),
      !any(ld$ld$r2 < .8), !grepl("r2 0.2-0.8 with", dline$evidence, fixed = TRUE), grepl("cohort WGS genotypes", dline$sources, fixed = TRUE))
    quiet_ld <- interp_evidence_summary(missense, list(clinvar = none_here, gwas_catalog = far), ld = alone)
    stopifnot(level(quiet_ld, "Disease and traits") == "None found",
      grepl("no other reported variant with r2 >= 0.2", quiet_ld$evidence[quiet_ld$question == "Disease and traits"], fixed = TRUE))
    # Protein abundance: a published protein-level association at a variant in tight LD counts; clinical traits do not.
    pl <- dl[dl$question == "Protein abundance (pQTL)", ]
    stopifnot(pl$level == "Moderate", grepl("through LD (r2 >= 0.8 in this cohort): protein x (rs150, r2 1.00)", pl$evidence, fixed = TRUE),
      grepl("blood level of 1 protein(s)", pl$finding, fixed = TRUE), level(quiet_ld, "Protein abundance (pQTL)") == "Not assessed",
      identical(interp_is_protein_trait(c("level of liver carboxylesterase 1 in blood", "Level of X in blood serum", "protein measurement",
        "high density lipoprotein cholesterol measurement", "apolipoprotein A 1 measurement")), c(TRUE, TRUE, FALSE, FALSE, FALSE)),
      interp_protein_name("level of liver carboxylesterase 1 in blood") == "liver carboxylesterase 1")
  } else if (native) stop("--native needs bcftools")
} else if (native) stop("--native needs tabix and bgzip") else cat("tabix/bgzip not found: indexed-source checks skipped\n")
unlink(work, recursive = TRUE)
cat("interpretation backend tests passed\n")
