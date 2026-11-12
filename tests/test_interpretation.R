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
  fhs_sqtl = c(discovery = "/missing", replication = "/missing"), fhs_eqtl = "/missing", spliceai = "/missing")
status <- interp_local_status(none)
stopifnot(nrow(status) == 6L, !any(status$available))
gone <- interp_local_evidence(v, none)
stopifnot(!gone$cadd$available, !gone$known$available, !gone$fhs$sqtl_available, !gone$fhs$eqtl_available, !gone$topmed$available,
  !nrow(gone$topmed$pairs), !nrow(gone$fhs$sqtl))

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
  st <- interp_local_status(srcs)
  stopifnot(st$available[st$source == "CADD v1.7"], st$available[st$source == "TOPMed multi-tissue cis-QTL"], !st$available[st$source == "SpliceAI 1.3.1"])
} else if (native) stop("--native needs tabix and bgzip") else cat("tabix/bgzip not found: indexed-source checks skipped\n")
unlink(work, recursive = TRUE)
cat("interpretation backend tests passed\n")
