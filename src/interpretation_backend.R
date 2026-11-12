# Variant interpretation for Reg_Shiny: what one selected DNA variant does to
# transcripts and proteins (computed here from GENCODE and the reference), what
# local and public resources report about it, and a rule-based evidence summary.
# Read-only: nothing here reads a BAM, changes a saved result or refits Model 1.

# ---- Sequence helpers -------------------------------------------------------
interp_complement <- function(x) chartr("ACGTNacgtn", "TGCANtgcan", x)

interp_codon_table <- function(mitochondrial = FALSE) {
  bases <- c("T", "C", "A", "G")
  # Third base varies fastest: TTT, TTC, TTA, TTG, TCT, ... (NCBI table order).
  grid <- expand.grid(third = bases, second = bases, first = bases, stringsAsFactors = FALSE)
  codons <- paste0(grid$first, grid$second, grid$third)
  aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1L]]
  table <- stats::setNames(aa, codons)
  if (mitochondrial) {
    # Vertebrate mitochondrial code (NCBI translation table 2).
    table[c("AGA", "AGG")] <- "*"; table[["ATA"]] <- "M"; table[["TGA"]] <- "W"
  }
  table
}

interp_aa3 <- c(A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", Q = "Gln", E = "Glu", G = "Gly",
  H = "His", I = "Ile", L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro", S = "Ser", T = "Thr",
  W = "Trp", Y = "Tyr", V = "Val", U = "Sec", `*` = "Ter", X = "Xaa")

interp_is_mito <- function(chrom) toupper(sub("^chr", "", chrom, ignore.case = TRUE)) %in% c("M", "MT")

# Sequence Ontology terms in Ensembl VEP severity order (most severe first).
interp_so_terms <- function() c(
  "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost",
  "start_lost", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
  "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant",
  "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
  "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "5_prime_UTR_variant",
  "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
  "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant")

interp_so_impact <- function(term) {
  high <- c("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost")
  moderate <- c("inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant")
  low <- c("splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant")
  ifelse(term %in% high, "HIGH", ifelse(term %in% moderate, "MODERATE", ifelse(term %in% low, "LOW", "MODIFIER")))
}

interp_most_severe <- function(terms) {
  terms <- unique(unlist(strsplit(terms, "&", fixed = TRUE)))
  rank <- match(terms, interp_so_terms())
  if (!length(terms) || all(is.na(rank))) return(NA_character_)
  terms[[which.min(rank)]]
}

# ---- Variant ---------------------------------------------------------------
# One ALT allele of a VCF record, trimmed to the bases that actually change.
# start1/end1 cover the changed reference bases; an insertion has end1 = start1 - 1
# and sits between start1 - 1 and start1.
interp_variant <- function(chrom, pos1, ref, alt) {
  if (length(chrom) != 1L || length(pos1) != 1L || length(ref) != 1L || length(alt) != 1L ||
      is.na(chrom) || is.na(pos1) || is.na(ref) || is.na(alt))
    stop("Interpretation needs one chromosome, position, REF and ALT allele.")
  ref <- toupper(ref); alt <- toupper(alt)
  if (grepl(",", alt, fixed = TRUE)) stop("Interpret one ALT allele at a time.")
  if (!grepl("^[ACGTN]+$", ref) || !grepl("^[ACGTN]+$", alt))
    stop("Only sequence alleles (A, C, G, T) can be interpreted; symbolic alleles are not.")
  if (identical(ref, alt)) stop("REF and ALT are identical.")
  pos1 <- as.integer(pos1); r <- ref; a <- alt; start <- pos1
  # Trim the shared suffix, then the shared prefix (VCF anchors indels on the left).
  last <- function(x) substring(x, nchar(x))
  while (nchar(r) > 1L && nchar(a) > 1L && last(r) == last(a)) {
    r <- substring(r, 1L, nchar(r) - 1L); a <- substring(a, 1L, nchar(a) - 1L)
  }
  while (nchar(r) > 0L && nchar(a) > 0L && substring(r, 1L, 1L) == substring(a, 1L, 1L)) {
    r <- substring(r, 2L); a <- substring(a, 2L); start <- start + 1L
  }
  type <- if (nchar(r) == 1L && nchar(a) == 1L) "SNV" else if (!nchar(r)) "insertion" else
    if (!nchar(a)) "deletion" else if (nchar(r) == nchar(a)) "MNV" else "delins"
  list(chrom = chrom, pos1 = pos1, ref = ref, alt = alt, type = type, start1 = start,
    end1 = start + nchar(r) - 1L, ref_changed = r, alt_changed = a,
    id = paste(chrom, pos1, ref, alt, sep = ":"))
}

# ---- Transcript models -----------------------------------------------------
# Exons in transcript order with their transcript (cDNA) coordinates, plus the
# coding segments (CDS and stop codon) in transcript order.
interp_transcript_model <- function(transcript, exons, features) {
  strand <- transcript$strand[[1L]]
  ex <- exons[order(exons$start1), c("start1", "end1"), drop = FALSE]
  ex <- unique(ex)
  if (identical(strand, "-")) ex <- ex[rev(seq_len(nrow(ex))), , drop = FALSE]
  ex$length <- ex$end1 - ex$start1 + 1L
  ex$t_start <- cumsum(c(0L, head(ex$length, -1L))) + 1L
  ex$t_end <- ex$t_start + ex$length - 1L
  rownames(ex) <- NULL
  cds <- features[features$feature %in% c("CDS", "stop_codon"), c("start1", "end1", "frame", "feature"), drop = FALSE]
  coding <- NULL; phase <- 0L
  if (nrow(cds) && any(cds$feature == "CDS")) {
    cds <- cds[order(cds$start1, cds$end1), , drop = FALSE]
    first <- cds[cds$feature == "CDS", , drop = FALSE]
    first <- if (identical(strand, "-")) first[which.max(first$end1), ] else first[which.min(first$start1), ]
    phase <- suppressWarnings(as.integer(first$frame[[1L]])); if (is.na(phase)) phase <- 0L
    # Merge CDS with its stop codon into contiguous coding segments.
    starts <- cds$start1; ends <- cds$end1; keep_s <- starts[[1L]]; keep_e <- ends[[1L]]
    out_s <- integer(); out_e <- integer()
    for (i in seq_along(starts)[-1L]) {
      if (starts[[i]] <= keep_e + 1L) keep_e <- max(keep_e, ends[[i]]) else {
        out_s <- c(out_s, keep_s); out_e <- c(out_e, keep_e); keep_s <- starts[[i]]; keep_e <- ends[[i]]
      }
    }
    coding <- data.frame(start1 = c(out_s, keep_s), end1 = c(out_e, keep_e))
    if (identical(strand, "-")) coding <- coding[rev(seq_len(nrow(coding))), , drop = FALSE]
    coding$length <- coding$end1 - coding$start1 + 1L
    coding$c_start <- cumsum(c(0L, head(coding$length, -1L))) + 1L
    rownames(coding) <- NULL
  }
  tags <- strsplit(as.character(transcript$tags[[1L]]), ";", fixed = TRUE)[[1L]]
  stop_feature <- any(features$feature == "stop_codon")
  list(strand = strand, exons = ex, coding = coding, phase = phase, tags = tags,
    has_stop_codon = stop_feature,
    selenocysteine = features[features$feature == "Selenocysteine", c("start1", "end1"), drop = FALSE])
}

# Transcript coordinate of an exonic genomic position (NA when not exonic).
interp_tpos <- function(model, pos1) {
  ex <- model$exons; i <- which(pos1 >= ex$start1 & pos1 <= ex$end1)
  if (!length(i)) return(NA_integer_)
  i <- i[[1L]]
  if (identical(model$strand, "-")) ex$t_start[[i]] + (ex$end1[[i]] - pos1) else ex$t_start[[i]] + (pos1 - ex$start1[[i]])
}

# 1-based coding index (CDS plus stop codon) of a genomic position, NA outside.
interp_cpos <- function(model, pos1) {
  cd <- model$coding
  if (is.null(cd)) return(NA_integer_)
  i <- which(pos1 >= cd$start1 & pos1 <= cd$end1)
  if (!length(i)) return(NA_integer_)
  i <- i[[1L]]
  if (identical(model$strand, "-")) cd$c_start[[i]] + (cd$end1[[i]] - pos1) else cd$c_start[[i]] + (pos1 - cd$start1[[i]])
}

# Genomic position of a 1-based coding index (NA beyond the coding sequence).
interp_cpos_to_genome <- function(model, cpos) {
  cd <- model$coding
  vapply(cpos, function(x) {
    i <- which(x >= cd$c_start & x <= cd$c_start + cd$length - 1L)
    if (!length(i)) return(NA_integer_)
    i <- i[[1L]]; offset <- x - cd$c_start[[i]]
    as.integer(if (identical(model$strand, "-")) cd$end1[[i]] - offset else cd$start1[[i]] + offset)
  }, integer(1))
}

# HGVS-style coding coordinate of an exonic transcript position: c.N, c.-N (5' UTR)
# or c.*N (3' UTR); n.N for a transcript without CDS.
interp_hgvs_position <- function(model, tpos) {
  cd <- model$coding
  if (is.null(cd)) return(list(prefix = "n.", label = as.character(tpos)))
  five <- if (identical(model$strand, "-")) cd$end1[[1L]] else cd$start1[[1L]]
  last <- nrow(cd)
  three <- if (identical(model$strand, "-")) cd$start1[[last]] else cd$end1[[last]]
  t_first <- interp_tpos(model, five); t_last <- interp_tpos(model, three)
  if (is.na(t_first) || is.na(t_last)) return(list(prefix = "c.", label = "?"))
  label <- if (tpos < t_first) paste0("-", t_first - tpos) else if (tpos > t_last) paste0("*", tpos - t_last) else
    as.character(tpos - t_first + 1L)
  list(prefix = "c.", label = label)
}

# Where one genomic position falls in a transcript.
# Returns exon/intron number (transcript order), distances into the intron from
# the donor (5') and acceptor (3') ends, and distance to the nearest
# intron-adjacent exon edge for exonic positions.
interp_locate <- function(model, pos1) {
  ex <- model$exons; n <- nrow(ex); minus <- identical(model$strand, "-")
  hit <- which(pos1 >= ex$start1 & pos1 <= ex$end1)
  if (length(hit)) {
    i <- hit[[1L]]
    # Edge next to the upstream intron (acceptor side) and downstream intron (donor side).
    to_five <- if (minus) ex$end1[[i]] - pos1 else pos1 - ex$start1[[i]]
    to_three <- if (minus) pos1 - ex$start1[[i]] else ex$end1[[i]] - pos1
    return(list(region = "exon", number = i, total = n,
      exon_to_acceptor = if (i > 1L) to_five + 1L else NA_integer_,
      exon_to_donor = if (i < n) to_three + 1L else NA_integer_,
      donor_distance = NA_integer_, acceptor_distance = NA_integer_, donor_exon = NA_integer_))
  }
  for (i in seq_len(max(0L, n - 1L))) {
    # Intron i lies between transcript-order exons i and i + 1.
    lo <- if (minus) ex$end1[[i + 1L]] + 1L else ex$end1[[i]] + 1L
    hi <- if (minus) ex$start1[[i]] - 1L else ex$start1[[i + 1L]] - 1L
    if (pos1 >= lo && pos1 <= hi) {
      donor <- if (minus) hi - pos1 + 1L else pos1 - lo + 1L
      acceptor <- if (minus) pos1 - lo + 1L else hi - pos1 + 1L
      return(list(region = "intron", number = i, total = n - 1L, exon_to_acceptor = NA_integer_,
        exon_to_donor = NA_integer_, donor_distance = donor, acceptor_distance = acceptor, donor_exon = i))
    }
  }
  list(region = "outside", number = NA_integer_, total = n, exon_to_acceptor = NA_integer_,
    exon_to_donor = NA_integer_, donor_distance = NA_integer_, acceptor_distance = NA_integer_, donor_exon = NA_integer_)
}

# Smallest finite distance, Inf when there is none.
interp_min <- function(...) { x <- c(...); x <- x[!is.na(x)]; if (length(x)) min(x) else Inf }

# Splice terms for an intronic position (VEP definitions).
interp_intron_terms <- function(donor, acceptor) {
  terms <- character()
  if (acceptor <= 2L) terms <- c(terms, "splice_acceptor_variant")
  if (donor <= 2L) terms <- c(terms, "splice_donor_variant")
  # Donor side: +5 is its own term, +3, +4 and +6 are the donor region, and
  # +7/+8 (like -3..-8 at the acceptor) are the general splice region.
  if (donor == 5L) terms <- c(terms, "splice_donor_5th_base_variant")
  if (donor %in% c(3L, 4L, 6L)) terms <- c(terms, "splice_donor_region_variant")
  if (donor %in% 7:8 || acceptor >= 3L && acceptor <= 8L) terms <- c(terms, "splice_region_variant")
  if (acceptor >= 3L && acceptor <= 17L) terms <- c(terms, "splice_polypyrimidine_tract_variant")
  essential <- any(c("splice_acceptor_variant", "splice_donor_variant") %in% terms)
  if (!essential || length(terms) > 1L) terms <- c(terms, "intron_variant")
  unique(terms)
}

# ---- Consequences ----------------------------------------------------------
interp_empty_consequences <- function() data.frame(
  gene_name = character(), gene_id = character(), transcript_id = character(), transcript_type = character(),
  transcript_flags = character(), strand = character(), consequence = character(), most_severe = character(),
  impact = character(), location = character(), hgvs_c = character(), hgvs_p = character(),
  codon_change = character(), amino_acids = character(), protein_position = integer(),
  donor_distance = integer(), acceptor_distance = integer(), distance_to_transcript = integer(),
  notes = character(), stringsAsFactors = FALSE)

# Consequence of one variant on one transcript. fetch(chrom, start1, end1) returns
# the reference sequence of that interval (plus strand) as one upper-case string.
interp_transcript_consequence <- function(variant, transcript, model, fetch, gene_window = 5000L) {
  notes <- character(); minus <- identical(model$strand, "-")
  row <- list(consequence = "", location = "", hgvs_c = "", hgvs_p = "", codon_change = "", amino_acids = "",
    protein_position = NA_integer_, donor_distance = NA_integer_, acceptor_distance = NA_integer_,
    distance_to_transcript = 0L)
  tx_start <- transcript$start1[[1L]]; tx_end <- transcript$end1[[1L]]
  first <- variant$start1; last <- max(variant$end1, variant$start1)
  if (variant$type == "insertion") { first <- variant$start1 - 1L; last <- variant$start1 }
  if (last < tx_start || first > tx_end) {
    distance <- if (last < tx_start) tx_start - last else first - tx_end
    if (distance > gene_window) return(NULL)
    upstream <- (last < tx_start) != minus
    row$consequence <- if (upstream) "upstream_gene_variant" else "downstream_gene_variant"
    row$location <- paste0(distance, " bp ", if (upstream) "upstream" else "downstream")
    row$distance_to_transcript <- as.integer(distance)
    return(c(row, list(notes = "")))
  }
  coding_tx <- !is.null(model$coding)
  if (!coding_tx) notes <- c(notes, "no annotated CDS")
  if ("cds_start_NF" %in% model$tags) notes <- c(notes, "cds_start_NF: codon numbers are relative to the annotated partial CDS")
  if ("cds_end_NF" %in% model$tags) notes <- c(notes, "cds_end_NF")
  if (variant$type != "SNV") {
    # Indels and multi-base changes: classify by what the changed bases overlap.
    positions <- seq.int(first, last)
    where <- lapply(positions, function(p) interp_locate(model, p))
    region <- vapply(where, function(w) w$region, character(1))
    cpos <- vapply(positions, function(p) interp_cpos(model, p), integer(1))
    terms <- character()
    if (variant$type == "insertion") {
      # Between two bases: coding only when both neighbours are coding.
      in_coding <- all(!is.na(cpos))
    } else in_coding <- any(!is.na(cpos))
    intronic <- which(region == "intron")
    for (k in intronic) terms <- c(terms, interp_intron_terms(where[[k]]$donor_distance, where[[k]]$acceptor_distance))
    if (in_coding) {
      shift <- nchar(variant$alt_changed) - nchar(variant$ref_changed)
      terms <- c(terms, if (shift %% 3L != 0L) "frameshift_variant" else if (shift > 0L) "inframe_insertion" else
        if (shift < 0L) "inframe_deletion" else "protein_altering_variant")
    } else if (any(region == "exon")) {
      exonic <- positions[region == "exon"][[1L]]
      terms <- c(terms, if (!coding_tx) "non_coding_transcript_exon_variant" else {
        five <- if (minus) model$coding$end1[[1L]] else model$coding$start1[[1L]]
        if ((exonic < five) != minus) "5_prime_UTR_variant" else "3_prime_UTR_variant"
      })
    }
    edge <- vapply(where, function(w) interp_min(w$exon_to_acceptor, w$exon_to_donor) <= 3L, logical(1))
    if (any(edge & region == "exon")) terms <- c(terms, "splice_region_variant")
    if (!length(terms)) terms <- "intron_variant"
    w <- where[[1L]]
    row$consequence <- paste(unique(terms), collapse = "&")
    row$location <- paste0(w$region, " ", w$number, "/", w$total)
    row$donor_distance <- w$donor_distance; row$acceptor_distance <- w$acceptor_distance
    notes <- c(notes, "non-SNV: classified by overlap only; no protein change is computed")
    return(c(row, list(notes = paste(notes, collapse = "; "))))
  }
  pos <- variant$start1
  ref <- variant$ref_changed; alt <- variant$alt_changed
  tx_ref <- if (minus) interp_complement(ref) else ref
  tx_alt <- if (minus) interp_complement(alt) else alt
  w <- interp_locate(model, pos)
  row$location <- paste0(w$region, " ", w$number, "/", w$total)
  if (w$region == "intron") {
    terms <- interp_intron_terms(w$donor_distance, w$acceptor_distance)
    row$donor_distance <- w$donor_distance; row$acceptor_distance <- w$acceptor_distance
    ex <- model$exons
    near_donor <- w$donor_distance <= w$acceptor_distance
    # Last base of the upstream exon, or first base of the downstream exon.
    anchor <- if (near_donor) (if (minus) ex$start1[[w$number]] else ex$end1[[w$number]]) else
      (if (minus) ex$end1[[w$number + 1L]] else ex$start1[[w$number + 1L]])
    h <- interp_hgvs_position(model, interp_tpos(model, anchor))
    row$hgvs_c <- paste0(h$prefix, h$label, if (near_donor) paste0("+", w$donor_distance) else paste0("-", w$acceptor_distance),
      tx_ref, ">", tx_alt)
    row$consequence <- paste(terms, collapse = "&")
    return(c(row, list(notes = paste(notes, collapse = "; "))))
  }
  # Exonic SNV.
  terms <- character()
  if (interp_min(w$exon_to_acceptor, w$exon_to_donor) <= 3L) terms <- "splice_region_variant"
  if (!is.na(w$exon_to_donor)) row$donor_distance <- -w$exon_to_donor
  if (!is.na(w$exon_to_acceptor)) row$acceptor_distance <- -w$exon_to_acceptor
  tpos <- interp_tpos(model, pos)
  h <- interp_hgvs_position(model, tpos)
  row$hgvs_c <- paste0(h$prefix, h$label, tx_ref, ">", tx_alt)
  cpos <- interp_cpos(model, pos)
  if (!coding_tx) {
    terms <- c("non_coding_transcript_exon_variant", terms)
  } else if (is.na(cpos)) {
    five <- if (minus) model$coding$end1[[1L]] else model$coding$start1[[1L]]
    terms <- c(if ((pos < five) != minus) "5_prime_UTR_variant" else "3_prime_UTR_variant", terms)
  } else {
    index0 <- cpos - 1L - model$phase
    coding_length <- sum(model$coding$length)
    if (index0 < 0L) {
      terms <- c("coding_sequence_variant", terms)
      notes <- c(notes, "position precedes the first complete codon of a partial CDS")
    } else {
      codon0 <- index0 %/% 3L
      cidx <- model$phase + 3L * codon0 + 1:3
      if (max(cidx) > coding_length) {
        terms <- c("incomplete_terminal_codon_variant", "coding_sequence_variant", terms)
      } else {
        genome <- interp_cpos_to_genome(model, cidx)
        bases <- vapply(genome, function(g) fetch(variant$chrom, g, g), character(1))
        if (minus) bases <- interp_complement(bases)
        within <- which(genome == pos)
        if (length(within) != 1L || any(!bases %in% c("A", "C", "G", "T"))) {
          terms <- c("coding_sequence_variant", terms)
          notes <- c(notes, "reference codon could not be read")
        } else {
          if (!identical(bases[[within]], tx_ref))
            notes <- c(notes, paste0("REF allele ", ref, " differs from the reference base at this position"))
          ref_codon <- bases; ref_codon[[within]] <- tx_ref
          alt_codon <- bases; alt_codon[[within]] <- tx_alt
          table <- interp_codon_table(interp_is_mito(variant$chrom))
          ref_aa <- unname(table[[paste(ref_codon, collapse = "")]])
          alt_aa <- unname(table[[paste(alt_codon, collapse = "")]])
          sec <- model$selenocysteine
          if (nrow(sec) && any(pos >= sec$start1 & pos <= sec$end1)) {
            ref_aa <- "U"; notes <- c(notes, "annotated selenocysteine codon")
          }
          residue <- codon0 + 1L
          # The last coding codon is a stop only when the transcript annotates one.
          start_known <- !("cds_start_NF" %in% model$tags)
          term <- if (identical(ref_aa, alt_aa)) {
            if (ref_aa == "*") "stop_retained_variant" else if (residue == 1L && start_known) "start_retained_variant" else "synonymous_variant"
          } else if (residue == 1L && start_known) "start_lost" else if (alt_aa == "*") "stop_gained" else
            if (ref_aa == "*") "stop_lost" else "missense_variant"
          terms <- c(term, terms)
          show <- function(codon, k) paste(ifelse(seq_along(codon) == k, codon, tolower(codon)), collapse = "")
          row$codon_change <- paste0(show(ref_codon, within), "/", show(alt_codon, within))
          row$amino_acids <- if (identical(ref_aa, alt_aa)) ref_aa else paste0(ref_aa, "/", alt_aa)
          row$protein_position <- residue
          three <- function(a) unname(interp_aa3[[a]])
          row$hgvs_p <- if (term == "start_lost") "p.Met1?" else if (identical(ref_aa, alt_aa)) paste0("p.", three(ref_aa), residue, "=") else
            if (term == "stop_lost") paste0("p.Ter", residue, three(alt_aa), "ext*?") else
            paste0("p.", three(ref_aa), residue, three(alt_aa))
        }
      }
    }
  }
  row$consequence <- paste(unique(terms), collapse = "&")
  c(row, list(notes = paste(notes, collapse = "; ")))
}

# Consequences of one variant on every transcript of a complete annotation query
# (query_annotation() output covering the variant +/- gene_window).
interp_consequences <- function(variant, annotation, fetch, gene_window = 5000L) {
  # The mitochondrial genome is too compact for a flanking-gene window to mean anything.
  if (interp_is_mito(variant$chrom)) gene_window <- 0L
  tx <- annotation$transcripts
  out <- interp_empty_consequences()
  if (is.null(tx) || !nrow(tx)) return(out)
  cache <- new.env(parent = emptyenv())
  cached_fetch <- function(chrom, start1, end1) {
    key <- paste(chrom, start1, end1, sep = ":")
    if (is.null(cache[[key]])) cache[[key]] <- fetch(chrom, start1, end1)
    cache[[key]]
  }
  rows <- list()
  for (i in seq_len(nrow(tx))) {
    id <- tx$transcript_id[[i]]
    model <- interp_transcript_model(tx[i, , drop = FALSE],
      annotation$exons[annotation$exons$transcript_id == id, , drop = FALSE],
      annotation$features[annotation$features$transcript_id == id, , drop = FALSE])
    r <- interp_transcript_consequence(variant, tx[i, , drop = FALSE], model, cached_fetch, gene_window)
    if (is.null(r)) next
    flags <- intersect(c("MANE_Select", "MANE_Plus_Clinical", "Ensembl_canonical", "basic"), model$tags)
    severe <- interp_most_severe(r$consequence)
    rows[[length(rows) + 1L]] <- data.frame(
      gene_name = tx$gene_name[[i]], gene_id = tx$gene_id[[i]], transcript_id = id,
      transcript_type = tx$transcript_type[[i]], transcript_flags = paste(flags, collapse = ";"),
      strand = model$strand, consequence = r$consequence, most_severe = severe, impact = interp_so_impact(severe),
      location = r$location, hgvs_c = r$hgvs_c, hgvs_p = r$hgvs_p, codon_change = r$codon_change,
      amino_acids = r$amino_acids, protein_position = as.integer(r$protein_position),
      donor_distance = as.integer(r$donor_distance), acceptor_distance = as.integer(r$acceptor_distance),
      distance_to_transcript = as.integer(r$distance_to_transcript), notes = r$notes, stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(out)
  out <- do.call(rbind, rows)
  flag_rank <- ifelse(grepl("MANE_Select", out$transcript_flags), 0L,
    ifelse(grepl("Ensembl_canonical", out$transcript_flags), 1L, ifelse(grepl("basic", out$transcript_flags), 2L, 3L)))
  out <- out[order(match(out$most_severe, interp_so_terms()), flag_rank, out$gene_name, out$transcript_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# One line per gene: the most severe consequence and the reference transcript's.
interp_gene_summary <- function(consequences) {
  if (!nrow(consequences)) return(data.frame(gene_name = character(), gene_id = character(),
    most_severe = character(), impact = character(), transcripts = integer(), transcripts_most_severe = integer(),
    reference_transcript = character(), reference_consequence = character(), reference_hgvs = character(),
    stringsAsFactors = FALSE))
  genes <- unique(consequences$gene_id)
  rows <- lapply(genes, function(g) {
    x <- consequences[consequences$gene_id == g, , drop = FALSE]
    severe <- interp_most_severe(x$most_severe)
    ref <- x[order(!grepl("MANE_Select", x$transcript_flags), !grepl("Ensembl_canonical", x$transcript_flags)), , drop = FALSE][1L, ]
    data.frame(gene_name = x$gene_name[[1L]], gene_id = g, most_severe = severe, impact = interp_so_impact(severe),
      transcripts = nrow(x), transcripts_most_severe = sum(x$most_severe == severe),
      reference_transcript = paste0(ref$transcript_id, if (nzchar(ref$transcript_flags)) paste0(" (", ref$transcript_flags, ")") else ""),
      reference_consequence = ref$consequence,
      reference_hgvs = trimws(paste(ref$hgvs_c, ref$hgvs_p)), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  # Equal severity: genes whose reference transcript is MANE or canonical first.
  reference_rank <- ifelse(grepl("MANE_Select", out$reference_transcript), 0L, ifelse(grepl("Ensembl_canonical", out$reference_transcript), 1L, 2L))
  out <- out[order(match(out$most_severe, interp_so_terms()), reference_rank, out$gene_name), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---- Small helpers for nested JSON -----------------------------------------
# Nested element or default; never fails on a missing branch.
interp_get <- function(x, ..., default = NA) {
  for (key in list(...)) {
    if (is.null(x) || !is.list(x)) return(default)
    x <- if (is.numeric(key)) (if (length(x) >= key) x[[key]] else NULL) else x[[key]]
  }
  if (is.null(x) || !length(x)) default else x
}

# One data frame from a list of records; f returns a named list of scalars.
interp_rows <- function(items, f, empty) {
  rows <- lapply(items, function(item) {
    r <- f(item)
    if (is.null(r)) return(NULL)
    as.data.frame(lapply(r, function(v) if (is.null(v) || !length(v)) NA else v[[1L]]), stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  out[, names(empty), drop = FALSE]
}

interp_join <- function(x, sep = "; ", max = Inf) {
  x <- unique(as.character(unlist(x))); x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("")
  extra <- length(x) - max
  if (extra > 0) paste0(paste(head(x, max), collapse = sep), sep, "+", extra, " more") else paste(x, collapse = sep)
}

interp_p_value <- function(mantissa, exponent) {
  if (is.null(mantissa) || is.null(exponent) || is.na(mantissa) || is.na(exponent)) return(c(p = NA_real_, log10 = NA_real_))
  c(p = mantissa * 10^exponent, log10 = log10(mantissa) + exponent)
}

# ---- Public sources --------------------------------------------------------
# Every request carries only a variant identifier (chromosome, position, alleles)
# or a public gene identifier. No sample, genotype count, cohort name or user
# detail leaves the cluster.
interp_http_json <- function(url, body = NULL, timeout = 45) {
  curl <- unname(Sys.which("curl"))
  if (!nzchar(curl)) stop("curl is not available for public database lookups.")
  args <- c("-sS", "-m", as.character(timeout), "-A", "Reg_Shiny", "-H", "Accept: application/json", "-w", "\n%{http_code}")
  if (!is.null(body)) {
    payload <- tempfile("regshiny_request_", fileext = ".json"); on.exit(unlink(payload), add = TRUE)
    writeLines(body, payload, useBytes = TRUE)
    args <- c(args, "-H", "Content-Type: application/json", "--data-binary", paste0("@", payload))
  }
  p <- processx::run(curl, c(args, url), timeout = timeout + 10, error_on_status = FALSE)
  host <- sub("^https?://([^/]+).*$", "\\1", url)
  if (isTRUE(p$timeout) || p$status != 0L) stop("No response from ", host, ": ", trimws(p$stderr), call. = FALSE)
  text <- sub("\n([0-9]{3})$", "", p$stdout)
  code <- suppressWarnings(as.integer(sub("^.*\n([0-9]{3})$", "\\1", p$stdout)))
  if (is.na(code)) stop("Unreadable response from ", host, ".", call. = FALSE)
  if (code == 404L) return(NULL)
  if (code >= 400L) stop(host, " returned HTTP ", code, ".", call. = FALSE)
  jsonlite::fromJSON(text, simplifyVector = FALSE)
}

interp_bare_chrom <- function(chrom) {
  x <- sub("^chr", "", chrom, ignore.case = TRUE)
  if (toupper(x) %in% c("M", "MT")) "MT" else x
}

interp_bare_gene <- function(gene_id) sub("\\..*$", "", gene_id)

# -- Open Targets Platform: variant ------------------------------------------
interp_open_targets_url <- function() "https://api.platform.opentargets.org/api/v4/graphql"

# Three requests: the variant itself, then its GWAS and its molecular-QTL credible
# sets page by page (one combined query exceeds the API's size limit).
interp_ot_variant_query <- function() paste(
  "query V($id:String!){ variant(variantId:$id){ id rsIds hgvsId variantDescription",
  "mostSevereConsequence{label}",
  "alleleFrequencies{populationName alleleFrequency}",
  "variantEffect{method assessment score normalisedScore target{approvedSymbol}}",
  "dbXrefs{id source}",
  "evidences(size:200){ count rows{ datasourceId datatypeId score clinicalSignificances confidence",
  " diseaseFromSource disease{id name} target{approvedSymbol} literature studyId } }",
  "pharmacogenomics{ genotype genotypeAnnotationText phenotypeText pgxCategory evidenceLevel drugs{drugFromSource} target{approvedSymbol} }",
  "} }")

interp_ot_credible_query <- function(kind) paste0(
  "query V($id:String!,$index:Int!,$size:Int!){ variant(variantId:$id){ credibleSets(studyTypes:[",
  if (kind == "gwas") "gwas" else "eqtl,sqtl,pqtl,tuqtl,sceqtl", "], page:{index:$index,size:$size}){ count rows{ studyLocusId studyType isTransQtl",
  " pValueMantissa pValueExponent beta finemappingMethod locus(variantIds:[$id]){rows{posteriorProbability is95CredibleSet}} ",
  if (kind == "gwas") paste("study{ id traitFromSource nSamples nCases publicationFirstAuthor pubmedId diseases{id name} }",
    "l2GPredictions(page:{index:0,size:1}){rows{score target{approvedSymbol}}}") else
    "study{ id projectId condition nSamples target{approvedSymbol id} biosample{biosampleId biosampleName} }",
  " } } } }")

interp_ot_empty <- function() list(
  gwas = data.frame(trait = character(), diseases = character(), p_value = numeric(), log10_p = numeric(), beta = numeric(),
    posterior_probability = numeric(), in_95_credible_set = logical(), finemapping = character(), l2g_gene = character(),
    l2g_score = numeric(), n_samples = numeric(), study_id = character(), first_author = character(), pubmed_id = character(),
    stringsAsFactors = FALSE),
  qtl = data.frame(qtl_type = character(), gene = character(), gene_id = character(), biosample = character(), p_value = numeric(),
    log10_p = numeric(), beta = numeric(), posterior_probability = numeric(), in_95_credible_set = logical(), trans = logical(),
    project = character(), condition = character(), study_id = character(), stringsAsFactors = FALSE),
  clinical = data.frame(source = character(), disease = character(), clinical_significance = character(), review_status = character(),
    gene = character(), score = numeric(), literature = character(), stringsAsFactors = FALSE),
  predictors = data.frame(method = character(), assessment = character(), score = numeric(), normalised_score = numeric(),
    gene = character(), stringsAsFactors = FALSE),
  frequencies = data.frame(population = character(), allele_frequency = numeric(), stringsAsFactors = FALSE),
  pharmacogenomics = data.frame(drugs = character(), gene = character(), genotype = character(), category = character(),
    evidence_level = character(), phenotype = character(), stringsAsFactors = FALSE))

interp_parse_ot_variant <- function(response) {
  empty <- interp_ot_empty()
  if (!is.null(response$errors)) stop("Open Targets: ", interp_get(response$errors, 1L, "message", default = "query failed"), call. = FALSE)
  v <- interp_get(response, "data", "variant", default = NULL)
  if (is.null(v)) return(c(list(found = FALSE, summary = list()), empty))
  locus <- function(x, field) interp_get(x, "locus", "rows", 1L, field)
  gwas <- interp_rows(interp_get(v, "gwas", "rows", default = list()), function(x) {
    p <- interp_p_value(interp_get(x, "pValueMantissa"), interp_get(x, "pValueExponent"))
    list(trait = interp_get(x, "study", "traitFromSource", default = ""),
      diseases = interp_join(lapply(interp_get(x, "study", "diseases", default = list()), function(d) d$name)),
      p_value = p[["p"]], log10_p = p[["log10"]], beta = interp_get(x, "beta"),
      posterior_probability = locus(x, "posteriorProbability"), in_95_credible_set = locus(x, "is95CredibleSet"),
      finemapping = interp_get(x, "finemappingMethod", default = ""),
      l2g_gene = interp_get(x, "l2GPredictions", "rows", 1L, "target", "approvedSymbol", default = ""),
      l2g_score = interp_get(x, "l2GPredictions", "rows", 1L, "score"),
      n_samples = interp_get(x, "study", "nSamples"), study_id = interp_get(x, "study", "id", default = ""),
      first_author = interp_get(x, "study", "publicationFirstAuthor", default = ""),
      pubmed_id = interp_get(x, "study", "pubmedId", default = ""))
  }, empty$gwas)
  qtl <- interp_rows(interp_get(v, "qtl", "rows", default = list()), function(x) {
    p <- interp_p_value(interp_get(x, "pValueMantissa"), interp_get(x, "pValueExponent"))
    list(qtl_type = interp_get(x, "studyType", default = ""),
      gene = interp_get(x, "study", "target", "approvedSymbol", default = ""),
      gene_id = interp_get(x, "study", "target", "id", default = ""),
      biosample = interp_get(x, "study", "biosample", "biosampleName", default = ""),
      p_value = p[["p"]], log10_p = p[["log10"]], beta = interp_get(x, "beta"),
      posterior_probability = locus(x, "posteriorProbability"), in_95_credible_set = locus(x, "is95CredibleSet"),
      trans = interp_get(x, "isTransQtl"), project = interp_get(x, "study", "projectId", default = ""),
      condition = interp_get(x, "study", "condition", default = ""), study_id = interp_get(x, "study", "id", default = ""))
  }, empty$qtl)
  clinical <- interp_rows(interp_get(v, "evidences", "rows", default = list()), function(x) {
    # GWAS credible sets are already listed in their own table.
    if (identical(interp_get(x, "datasourceId"), "gwas_credible_sets")) return(NULL)
    list(source = interp_get(x, "datasourceId", default = ""),
      disease = interp_get(x, "disease", "name", default = interp_get(x, "diseaseFromSource", default = "")),
      clinical_significance = interp_join(interp_get(x, "clinicalSignificances", default = list())),
      review_status = interp_get(x, "confidence", default = ""),
      gene = interp_get(x, "target", "approvedSymbol", default = ""), score = interp_get(x, "score"),
      literature = interp_join(interp_get(x, "literature", default = list()), max = 5))
  }, empty$clinical)
  predictors <- interp_rows(interp_get(v, "variantEffect", default = list()), function(x)
    list(method = interp_get(x, "method", default = ""), assessment = interp_get(x, "assessment", default = ""),
      score = interp_get(x, "score"), normalised_score = interp_get(x, "normalisedScore"),
      gene = interp_get(x, "target", "approvedSymbol", default = "")), empty$predictors)
  frequencies <- interp_rows(interp_get(v, "alleleFrequencies", default = list()), function(x)
    list(population = interp_get(x, "populationName", default = ""), allele_frequency = interp_get(x, "alleleFrequency")),
    empty$frequencies)
  pgx <- interp_rows(interp_get(v, "pharmacogenomics", default = list()), function(x)
    list(drugs = interp_join(lapply(interp_get(x, "drugs", default = list()), function(d) d$drugFromSource)),
      gene = interp_get(x, "target", "approvedSymbol", default = ""), genotype = interp_get(x, "genotype", default = ""),
      category = interp_get(x, "pgxCategory", default = ""), evidence_level = interp_get(x, "evidenceLevel", default = ""),
      phenotype = interp_get(x, "phenotypeText", default = "")), empty$pharmacogenomics)
  if (nrow(gwas)) gwas <- gwas[order(-ifelse(is.na(gwas$posterior_probability), -1, gwas$posterior_probability), gwas$log10_p), , drop = FALSE]
  if (nrow(qtl)) qtl <- qtl[order(-ifelse(is.na(qtl$posterior_probability), -1, qtl$posterior_probability), qtl$log10_p), , drop = FALSE]
  rownames(gwas) <- NULL; rownames(qtl) <- NULL
  xref <- interp_get(v, "dbXrefs", default = list())
  list(found = TRUE,
    summary = list(id = interp_get(v, "id", default = ""), rs_ids = interp_join(interp_get(v, "rsIds", default = list())),
      description = interp_get(v, "variantDescription", default = ""),
      most_severe = interp_get(v, "mostSevereConsequence", "label", default = ""),
      clinvar_ids = interp_join(lapply(xref, function(z) if (identical(z$source, "clinvar")) z$id else NULL), max = 8),
      gwas_total = interp_get(v, "gwas", "count", default = 0), qtl_total = interp_get(v, "qtl", "count", default = 0),
      evidence_total = interp_get(v, "evidences", "count", default = 0)),
    gwas = gwas, qtl = qtl, clinical = clinical, predictors = predictors, frequencies = frequencies, pharmacogenomics = pgx)
}

interp_open_targets_variant <- function(variant, http = interp_http_json, page_size = 500L, max_pages = 4L) {
  id <- paste(interp_bare_chrom(variant$chrom), variant$pos1, variant$ref, variant$alt, sep = "_")
  ask <- function(query, variables) {
    r <- http(interp_open_targets_url(), body = jsonlite::toJSON(list(query = query, variables = variables), auto_unbox = TRUE))
    if (!is.null(r$errors)) stop("Open Targets: ", interp_get(r$errors, 1L, "message", default = "query failed"), call. = FALSE)
    r
  }
  core <- ask(interp_ot_variant_query(), list(id = id))
  if (is.null(interp_get(core, "data", "variant", default = NULL))) return(interp_parse_ot_variant(core))
  for (kind in c("gwas", "qtl")) {
    rows <- list(); count <- 0L
    for (index in seq_len(max_pages) - 1L) {
      page <- interp_get(ask(interp_ot_credible_query(kind), list(id = id, index = index, size = page_size)),
        "data", "variant", "credibleSets", default = list())
      count <- interp_get(page, "count", default = 0L)
      rows <- c(rows, interp_get(page, "rows", default = list()))
      if (length(rows) >= count || !length(interp_get(page, "rows", default = list()))) break
    }
    core$data$variant[[kind]] <- list(count = count, rows = rows)
  }
  interp_parse_ot_variant(core)
}

# -- Open Targets Platform: genes (function, pathways, GO, diseases) ----------
interp_ot_target_query <- function() paste(
  "query T($ids:[String!]!){ targets(ensemblIds:$ids){ id approvedSymbol approvedName biotype functionDescriptions",
  "proteinIds{id source} subcellularLocations{location}",
  "pathways{pathway topLevelTerm pathwayId}",
  "geneOntology{aspect evidence term{id label}}",
  "geneticConstraint{constraintType oe oeLower oeUpper score upperBin6}",
  "associatedDiseases(page:{index:0,size:15}){count rows{score disease{id name}}} } }")

interp_target_empty <- function() list(
  genes = data.frame(gene = character(), gene_id = character(), name = character(), biotype = character(), uniprot = character(),
    function_description = character(), subcellular_location = character(), lof_oe_upper = numeric(), lof_constraint_bin = numeric(),
    associated_diseases = integer(), stringsAsFactors = FALSE),
  pathways = data.frame(gene = character(), top_level = character(), pathway = character(), pathway_id = character(), stringsAsFactors = FALSE),
  go = data.frame(gene = character(), aspect = character(), term = character(), term_id = character(), evidence = character(), stringsAsFactors = FALSE),
  diseases = data.frame(gene = character(), disease = character(), disease_id = character(), score = numeric(), stringsAsFactors = FALSE))

interp_parse_ot_targets <- function(response) {
  empty <- interp_target_empty()
  if (!is.null(response$errors)) stop("Open Targets: ", interp_get(response$errors, 1L, "message", default = "query failed"), call. = FALSE)
  targets <- interp_get(response, "data", "targets", default = list())
  out <- empty
  for (t in targets) {
    symbol <- interp_get(t, "approvedSymbol", default = "")
    lof <- Filter(function(z) identical(z$constraintType, "lof"), interp_get(t, "geneticConstraint", default = list()))
    uniprot <- lapply(interp_get(t, "proteinIds", default = list()), function(z) if (identical(z$source, "uniprot_swissprot")) z$id else NULL)
    out$genes <- rbind(out$genes, data.frame(gene = symbol, gene_id = interp_get(t, "id", default = ""),
      name = interp_get(t, "approvedName", default = ""), biotype = interp_get(t, "biotype", default = ""),
      uniprot = interp_join(uniprot), function_description = interp_join(interp_get(t, "functionDescriptions", default = list()), sep = " ", max = 2),
      subcellular_location = interp_join(lapply(interp_get(t, "subcellularLocations", default = list()), function(z) z$location), max = 6),
      lof_oe_upper = as.numeric(interp_get(lof, 1L, "oeUpper")), lof_constraint_bin = as.numeric(interp_get(lof, 1L, "upperBin6")),
      associated_diseases = as.integer(interp_get(t, "associatedDiseases", "count", default = 0L)), stringsAsFactors = FALSE))
    out$pathways <- rbind(out$pathways, interp_rows(interp_get(t, "pathways", default = list()), function(z)
      list(gene = symbol, top_level = interp_get(z, "topLevelTerm", default = ""), pathway = interp_get(z, "pathway", default = ""),
        pathway_id = interp_get(z, "pathwayId", default = "")), empty$pathways))
    out$go <- rbind(out$go, interp_rows(interp_get(t, "geneOntology", default = list()), function(z)
      list(gene = symbol, aspect = c(F = "Molecular function", P = "Biological process", C = "Cellular component")[[
        match.arg(interp_get(z, "aspect", default = "P"), c("F", "P", "C"))]],
        term = interp_get(z, "term", "label", default = ""), term_id = interp_get(z, "term", "id", default = ""),
        evidence = interp_get(z, "evidence", default = "")), empty$go))
    out$diseases <- rbind(out$diseases, interp_rows(interp_get(t, "associatedDiseases", "rows", default = list()), function(z)
      list(gene = symbol, disease = interp_get(z, "disease", "name", default = ""), disease_id = interp_get(z, "disease", "id", default = ""),
        score = interp_get(z, "score")), empty$diseases))
  }
  if (nrow(out$go)) out$go <- unique(out$go)
  if (nrow(out$pathways)) out$pathways <- unique(out$pathways)
  out
}

interp_open_targets_genes <- function(gene_ids, http = interp_http_json) {
  ids <- unique(interp_bare_gene(gene_ids)); ids <- ids[grepl("^ENSG[0-9]+$", ids)]
  if (!length(ids)) return(interp_target_empty())
  body <- jsonlite::toJSON(list(query = interp_ot_target_query(), variables = list(ids = as.list(ids))), auto_unbox = TRUE)
  interp_parse_ot_targets(http(interp_open_targets_url(), body = body))
}

# -- GTEx Portal: tissue eQTL / sQTL and expression ---------------------------
interp_gtex_url <- function(path, query) paste0("https://gtexportal.org/api/v2/", path, "?",
  paste(paste0(names(query), "=", vapply(query, function(v) utils::URLencode(as.character(v), reserved = TRUE), character(1))), collapse = "&"))

interp_gtex_dataset <- function() "gtex_v10"

interp_gtex_empty <- function() list(
  qtl = data.frame(qtl_type = character(), tissue = character(), gene = character(), gene_id = character(), phenotype = character(),
    p_value = numeric(), nes = numeric(), stringsAsFactors = FALSE),
  expression = data.frame(gene = character(), tissue = character(), median_tpm = numeric(), stringsAsFactors = FALSE))

interp_gtex_all <- function(url, http) {
  # One page holds up to 250 rows; follow the paging the API reports.
  rows <- list(); page <- 0L
  repeat {
    r <- http(paste0(url, "&itemsPerPage=250&page=", page))
    if (is.null(r)) break
    rows <- c(rows, interp_get(r, "data", default = list()))
    pages <- interp_get(r, "paging_info", "numberOfPages", default = 1L)
    page <- page + 1L
    if (page >= pages || page >= 8L) break
  }
  rows
}

interp_gtex <- function(variant, gene_ids = character(), http = interp_http_json) {
  out <- interp_gtex_empty()
  out$dataset <- interp_gtex_dataset(); out$variant_in_gtex <- NA; out$genes_mapped <- character()
  # GTEx only holds SNVs and short indels on chr1-22 and chrX.
  id <- paste0(paste(variant$chrom, variant$pos1, variant$ref, variant$alt, sep = "_"), "_b38")
  found <- http(interp_gtex_url("dataset/variant", list(variantId = id, datasetId = out$dataset)))
  out$variant_in_gtex <- length(interp_get(found, "data", default = list())) > 0L
  tissue <- function(x) gsub("_", " ", interp_get(x, "tissueSiteDetailId", default = ""))
  if (isTRUE(out$variant_in_gtex)) for (type in c("eQTL", "sQTL")) {
    path <- if (type == "eQTL") "association/singleTissueEqtl" else "association/singleTissueSqtl"
    rows <- interp_gtex_all(interp_gtex_url(path, list(variantId = id, datasetId = out$dataset)), http)
    out$qtl <- rbind(out$qtl, interp_rows(rows, function(x) list(qtl_type = type, tissue = tissue(x),
      gene = interp_get(x, "geneSymbol", default = ""), gene_id = interp_get(x, "gencodeId", default = ""),
      phenotype = interp_get(x, "phenotypeId", default = ""), p_value = interp_get(x, "pValue"), nes = interp_get(x, "nes")),
      interp_gtex_empty()$qtl))
  }
  if (nrow(out$qtl)) { out$qtl <- out$qtl[order(out$qtl$p_value), , drop = FALSE]; rownames(out$qtl) <- NULL }
  ids <- unique(interp_bare_gene(gene_ids)); ids <- head(ids[grepl("^ENSG[0-9]+$", ids)], 5L)
  for (g in ids) {
    # The expression endpoint needs the gene version of the GTEx release.
    ref <- http(interp_gtex_url("reference/gene", list(geneId = g, gencodeVersion = "v39", genomeBuild = "GRCh38/hg38")))
    versioned <- interp_get(ref, "data", 1L, "gencodeId", default = NA)
    if (is.na(versioned)) next
    out$genes_mapped <- c(out$genes_mapped, versioned)
    rows <- interp_get(http(interp_gtex_url("expression/medianGeneExpression", list(gencodeId = versioned, datasetId = out$dataset))),
      "data", default = list())
    out$expression <- rbind(out$expression, interp_rows(rows, function(x) list(gene = interp_get(x, "geneSymbol", default = g),
      tissue = tissue(x), median_tpm = interp_get(x, "median")), interp_gtex_empty()$expression))
  }
  out
}

# Tissue summary of one gene's GTEx medians: rank, share of the top tissue and
# tau, the usual 0-1 tissue-specificity index (1 = one tissue only).
interp_tissue_specificity <- function(expression) {
  empty <- data.frame(gene = character(), tissues = integer(), top_tissue = character(), top_median_tpm = numeric(),
    whole_blood_tpm = numeric(), whole_blood_rank = integer(), tau = numeric(), stringsAsFactors = FALSE)
  if (is.null(expression) || !nrow(expression)) return(empty)
  rows <- lapply(split(expression, expression$gene), function(x) {
    x <- x[is.finite(x$median_tpm), , drop = FALSE]
    if (!nrow(x)) return(NULL)
    v <- log2(x$median_tpm + 1); top <- which.max(x$median_tpm)
    blood <- which(x$tissue == "Whole Blood")
    data.frame(gene = x$gene[[1L]], tissues = nrow(x), top_tissue = x$tissue[[top]], top_median_tpm = x$median_tpm[[top]],
      whole_blood_tpm = if (length(blood)) x$median_tpm[[blood[[1L]]]] else NA_real_,
      whole_blood_rank = if (length(blood)) as.integer(rank(-x$median_tpm, ties.method = "min")[[blood[[1L]]]]) else NA_integer_,
      tau = if (max(v) > 0 && nrow(x) > 1L) sum(1 - v / max(v)) / (nrow(x) - 1L) else NA_real_, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# -- MyVariant.info: ClinVar and dbSNP ----------------------------------------
interp_hgvs_g <- function(variant) {
  chrom <- paste0("chr", interp_bare_chrom(variant$chrom))
  r <- variant$ref_changed; a <- variant$alt_changed; s <- variant$start1; e <- variant$end1
  switch(variant$type,
    SNV = paste0(chrom, ":g.", s, r, ">", a),
    deletion = if (s == e) paste0(chrom, ":g.", s, "del") else paste0(chrom, ":g.", s, "_", e, "del"),
    insertion = paste0(chrom, ":g.", s - 1L, "_", s, "ins", a),
    if (s == e) paste0(chrom, ":g.", s, "delins", a) else paste0(chrom, ":g.", s, "_", e, "delins", a))
}

interp_myvariant_empty <- function() list(found = FALSE, rs_id = "", clinvar_variant_id = "", gnomad_af = NA_real_,
  revel = NA_real_, clinvar = data.frame(accession = character(), clinical_significance = character(), condition = character(),
    review_status = character(), last_evaluated = character(), stringsAsFactors = FALSE))

interp_parse_myvariant <- function(response) {
  out <- interp_myvariant_empty()
  if (is.null(response) || !is.null(response$error) || isFALSE(response$success)) return(out)
  out$found <- TRUE
  out$rs_id <- as.character(interp_get(response, "dbsnp", "rsid", default = interp_get(response, "clinvar", "rsid", default = "")))
  out$clinvar_variant_id <- as.character(interp_get(response, "clinvar", "variant_id", default = ""))
  first_number <- function(x) { x <- suppressWarnings(as.numeric(unlist(x))); x <- x[is.finite(x)]; if (length(x)) x[[1L]] else NA_real_ }
  out$gnomad_af <- first_number(interp_get(response, "gnomad_genome", "af", "af", default = interp_get(response, "gnomad_exome", "af", "af")))
  out$revel <- first_number(interp_get(response, "dbnsfp", "revel", "score"))
  rcv <- interp_get(response, "clinvar", "rcv", default = list())
  # A single record arrives as one object rather than a list of them.
  if (!is.null(rcv$accession)) rcv <- list(rcv)
  out$clinvar <- interp_rows(rcv, function(x) {
    cond <- interp_get(x, "conditions", default = list())
    names <- if (!is.null(cond$name)) cond$name else lapply(cond, function(z) z$name)
    list(accession = interp_get(x, "accession", default = ""), clinical_significance = interp_get(x, "clinical_significance", default = ""),
      condition = interp_join(names), review_status = interp_get(x, "review_status", default = ""),
      last_evaluated = interp_get(x, "last_evaluated", default = ""))
  }, out$clinvar)
  out
}

interp_myvariant <- function(variant, http = interp_http_json) {
  fields <- "dbsnp.rsid,clinvar.rsid,clinvar.variant_id,clinvar.rcv,gnomad_genome.af.af,gnomad_exome.af.af,dbnsfp.revel.score"
  url <- paste0("https://myvariant.info/v1/variant/", utils::URLencode(interp_hgvs_g(variant), reserved = TRUE),
    "?assembly=hg38&fields=", fields)
  interp_parse_myvariant(http(url))
}

# All public lookups for one variant. Each source fails on its own: an error is
# kept as that source's value and the others are still returned.
interp_public_evidence <- function(variant, gene_ids = character(), http = interp_http_json) {
  attempt <- function(expr) tryCatch(expr, error = function(e) e)
  list(queried_at = format(Sys.time(), tz = "UTC", usetz = TRUE), variant = variant$id, gene_ids = unique(interp_bare_gene(gene_ids)),
    open_targets = attempt(interp_open_targets_variant(variant, http)),
    genes = attempt(interp_open_targets_genes(gene_ids, http)),
    gtex = attempt(interp_gtex(variant, gene_ids, http)),
    myvariant = attempt(interp_myvariant(variant, http)))
}

# ---- Local sources on SCC --------------------------------------------------
# Paths can be overridden in config/runtime.Renviron; a source that is missing or
# unreadable is reported as unavailable and skipped, never guessed.
interp_local_sources <- function() {
  env <- function(name, default) { v <- Sys.getenv(name, ""); if (nzchar(v)) v else default }
  first_existing <- function(paths) { p <- paths[nzchar(paths) & file.exists(paths)]; if (length(p)) p[[1L]] else "" }
  lab <- "/restricted/projectnb/mtdna-alcohol"
  sqtl <- strsplit(env("REGSHINY_FHS_SQTL", paste(
    file.path(lab, "topmed-qtl/sqtl/sqtl-batch12-combined-result-1e-4-annot-filtered-meta.csv"),
    file.path(lab, "topmed-qtl/sqtl/sqtl-batch3-combined-result-1e-4-annot-filtered-meta.csv"), sep = ";")), ";", fixed = TRUE)[[1L]]
  list(
    tabix = first_existing(c(Sys.getenv("REGSHINY_TABIX", ""), unname(Sys.which("tabix")),
      rev(Sys.glob("/share/pkg.8/htslib/*/install/bin/tabix")))),
    look = unname(Sys.which("look")),
    cadd_snv = env("REGSHINY_CADD_SNV", file.path(lab, "Jinjie_Jay/codex_shared_resources/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz")),
    cadd_indel = env("REGSHINY_CADD_INDEL", file.path(lab, "Jinjie_Jay/codex_shared_resources/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz")),
    vep_cache = env("REGSHINY_VEP_CACHE", "/share/pkg.8/vep/112.0-singularity/install/data/homo_sapiens/112_GRCh38"),
    topmed_qtl = env("REGSHINY_TOPMED_QTL", file.path(lab, "projectnb/RNAseq/TOPMed_cis-sQTL")),
    fhs_sqtl = stats::setNames(sqtl, c("discovery", "replication")[seq_along(sqtl)]),
    fhs_eqtl = env("REGSHINY_FHS_EQTL", file.path(lab, "Jinjie_Jay/Reg_Shiny_output/resources/topmed_qtl/fhs_whole_blood_eqtl_gene_combined_1e-7.sorted.csv")),
    gwas_catalog = env("REGSHINY_GWAS_CATALOG", file.path(lab, "Jinjie_Jay/Reg_Shiny_output/resources/disease/gwas_catalog_e114_r2025-07-21.b38.tsv.gz")),
    clinvar_vcf = env("REGSHINY_CLINVAR_VCF", file.path(lab, "Jinjie_Jay/Reg_Shiny_output/resources/disease/clinvar_2025-07-21.GRCh38.vcf.gz")),
    reactome = env("REGSHINY_REACTOME", file.path(lab, "Jinjie_Jay/Reg_Shiny_output/resources/pathway/reactome_ensembl_all_levels_2026-09-22.human.tsv")),
    msigdb = env("REGSHINY_MSIGDB", "/restricted/projectnb/cbmhive/MSigDB/2026.1.Hs"),
    tissue_model = c(eQTL = env("REGSHINY_TISSUE_MODEL_EQTL", file.path(lab, "Jinjie_Jay/Reg_Shiny_output/resources/topmed_qtl/tissue_model/tissue_sharing_eqtl.json")),
      sQTL = env("REGSHINY_TISSUE_MODEL_SQTL", file.path(lab, "Jinjie_Jay/Reg_Shiny_output/resources/topmed_qtl/tissue_model/tissue_sharing_sqtl.json"))),
    spliceai = env("REGSHINY_SPLICEAI", "/share/pkg.8/spliceai/1.3.1/install/bin/spliceai"))
}

interp_readable <- function(path) length(path) == 1L && !is.na(path) && nzchar(path) && file.exists(path) && file.access(path, 4) == 0L

interp_local_status <- function(sources = interp_local_sources()) {
  tabix <- interp_readable(sources$tabix)
  row <- function(source, content, ok, path) data.frame(source = source, content = content,
    available = ok, path = paste(path, collapse = "; "), stringsAsFactors = FALSE)
  rbind(
    row("CADD v1.7", "Deleteriousness score of every possible SNV; gnomAD indels", tabix && interp_readable(sources$cadd_snv), sources$cadd_snv),
    row("Ensembl VEP 112 variation cache", "rsID, ClinVar significance (2023-10), 1000 Genomes and gnomAD allele frequencies",
      tabix && dir.exists(sources$vep_cache), sources$vep_cache),
    row("TOPMed multi-tissue cis-QTL", "cis-eQTL and cis-sQTL in lung, monocytes, nasal epithelium, PBMC, T cells and whole blood",
      tabix && dir.exists(file.path(sources$topmed_qtl, "cis-eqtl", "nominal")), sources$topmed_qtl),
    row("Cross-context sharing model", "Model trained on TOPMed lead cis-QTLs: effect and its certainty in each context (scripts/train_tissue_model.R)",
      length(sources$tissue_model) > 0L && any(vapply(sources$tissue_model, interp_readable, logical(1))), sources$tissue_model),
    row("FHS whole-blood isoform-ratio sQTL", "Discovery and replication sets, p < 1e-4",
      nzchar(sources$look) && all(vapply(sources$fhs_sqtl, interp_readable, logical(1))), sources$fhs_sqtl),
    row("FHS whole-blood gene eQTL", "Combined set, p < 1e-7 (position-sorted copy)",
      nzchar(sources$look) && interp_readable(sources$fhs_eqtl), sources$fhs_eqtl),
    row("GWAS Catalog", "Published trait associations reported at the variant's position (position-indexed copy)",
      tabix && interp_readable(sources$gwas_catalog) && interp_readable(paste0(sources$gwas_catalog, ".tbi")), sources$gwas_catalog),
    row("ClinVar VCF", "Germline classification, review status and conditions of this allele",
      tabix && interp_readable(sources$clinvar_vcf) && interp_readable(paste0(sources$clinvar_vcf, ".tbi")), sources$clinvar_vcf),
    row("Reactome gene-to-pathway mapping", "Pathways of the genes the variant overlaps or is a QTL for (gene-sorted copy)",
      interp_readable(sources$look) && interp_readable(sources$reactome), sources$reactome),
    row("MSigDB gene sets", "Hallmark, canonical pathways (KEGG, WikiPathways, PID, BioCarta) and GO biological process",
      length(interp_msigdb_files(sources)) > 0L, sources$msigdb),
    row("SpliceAI 1.3.1", "Predicted gain or loss of splice sites within 500 bp (run on request)", interp_readable(sources$spliceai), sources$spliceai))
}

interp_run_lines <- function(command, args, timeout = 120, env = NULL) {
  p <- processx::run(command, args, timeout = timeout, error_on_status = FALSE, env = env)
  if (isTRUE(p$timeout)) stop(basename(command), " timed out.", call. = FALSE)
  # tabix and look return 1 for "nothing found"; real failures write to stderr.
  if (!p$status %in% c(0L, 1L) || (p$status == 1L && grepl("rror|ail", p$stderr) && !nzchar(p$stdout)))
    stop(basename(command), " failed: ", trimws(p$stderr), call. = FALSE)
  lines <- strsplit(p$stdout, "\n", fixed = TRUE)[[1L]]
  lines[nzchar(lines)]
}

interp_tabix <- function(sources, file, chrom, start1, end1 = start1)
  interp_run_lines(sources$tabix, c(file, paste0(chrom, ":", start1, "-", end1)))

# CADD v1.7 score of this exact allele (files use bare chromosome names).
interp_cadd <- function(variant, sources = interp_local_sources()) {
  out <- list(available = FALSE, found = FALSE, raw = NA_real_, phred = NA_real_, source = "CADD v1.7 (GRCh38)")
  file <- if (variant$type == "SNV") sources$cadd_snv else sources$cadd_indel
  if (!interp_readable(sources$tabix) || !interp_readable(file) || interp_is_mito(variant$chrom)) return(out)
  out$available <- TRUE
  # CADD lists an SNV at its own position and an indel in VCF form.
  snv <- variant$type == "SNV"
  pos <- if (snv) variant$start1 else variant$pos1
  ref <- if (snv) variant$ref_changed else variant$ref; alt <- if (snv) variant$alt_changed else variant$alt
  for (line in interp_tabix(sources, file, interp_bare_chrom(variant$chrom), pos)) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(f) >= 6L && f[[2L]] == as.character(pos) && f[[3L]] == ref && f[[4L]] == alt) {
      out$found <- TRUE; out$raw <- as.numeric(f[[5L]]); out$phred <- as.numeric(f[[6L]]); break
    }
  }
  out
}

# dbSNP record of this allele from the Ensembl VEP 112 cache: rsID, ClinVar
# significance and population allele frequencies. Works without network.
interp_known_variant <- function(variant, sources = interp_local_sources()) {
  out <- list(available = FALSE, found = FALSE, rs_id = "", clinvar_significance = "", pubmed_n = 0L,
    frequencies = data.frame(population = character(), allele = character(), allele_frequency = numeric(), stringsAsFactors = FALSE),
    source = "Ensembl VEP 112 cache (dbSNP 156, ClinVar 2023-10, 1000 Genomes phase 3, gnomAD v3.1.2 genomes)")
  chrom <- interp_bare_chrom(variant$chrom)
  file <- file.path(sources$vep_cache, chrom, "all_vars.gz")
  if (!interp_readable(sources$tabix) || !interp_readable(file)) return(out)
  out$available <- TRUE
  info <- file.path(sources$vep_cache, "info.txt")
  cols <- c("chr", "variation_name", "failed", "somatic", "start", "end", "allele_string", "strand", "clin_sig",
    "phenotype_or_disease", "clin_sig_allele", "pubmed", "var_synonyms", "AF", "AFR", "AMR", "EAS", "EUR", "SAS")
  if (file.exists(info)) {
    line <- grep("^variation_cols\t", readLines(info, warn = FALSE), value = TRUE)
    if (length(line)) cols <- strsplit(sub("^variation_cols\t", "", line[[1L]]), ",", fixed = TRUE)[[1L]]
  }
  # The cache writes an indel with its changed bases only ("-" for none) at the first changed base.
  ref <- if (nzchar(variant$ref_changed)) variant$ref_changed else "-"
  alt <- if (nzchar(variant$alt_changed)) variant$alt_changed else "-"
  start <- if (variant$type == "insertion") variant$start1 - 1L else variant$start1
  for (line in interp_tabix(sources, file, chrom, min(start, variant$start1), max(variant$end1, variant$start1))) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(f) < length(cols)) f <- c(f, rep(".", length(cols) - length(f)))
    x <- stats::setNames(as.list(f[seq_along(cols)]), cols)
    alleles <- strsplit(x$allele_string, "/", fixed = TRUE)[[1L]]
    if (!startsWith(x$variation_name, "rs") || !identical(alleles[[1L]], ref) || !alt %in% alleles[-1L]) next
    out$found <- TRUE; out$rs_id <- x$variation_name
    # Only the significance entries recorded for this ALT allele.
    sig <- strsplit(if (x$clin_sig_allele == ".") "" else x$clin_sig_allele, ";", fixed = TRUE)[[1L]]
    sig <- sub("^[^:]*:", "", sig[startsWith(sig, paste0(alt, ":"))])
    out$clinvar_significance <- interp_join(unlist(strsplit(sig, "/", fixed = TRUE)), sep = ", ")
    if (x$pubmed != ".") out$pubmed_n <- length(strsplit(x$pubmed, ",", fixed = TRUE)[[1L]])
    freq <- list()
    for (name in setdiff(cols[-(1:13)], "")) {
      value <- x[[name]]; if (is.null(value) || value %in% c(".", "")) next
      for (part in strsplit(value, ",", fixed = TRUE)[[1L]]) {
        kv <- strsplit(part, ":", fixed = TRUE)[[1L]]
        allele <- if (length(kv) == 2L) kv[[1L]] else alt; number <- suppressWarnings(as.numeric(kv[[length(kv)]]))
        if (identical(allele, alt) && is.finite(number))
          freq[[length(freq) + 1L]] <- data.frame(population = name, allele = allele, allele_frequency = number, stringsAsFactors = FALSE)
      }
    }
    if (length(freq)) out$frequencies <- do.call(rbind, freq)
    break
  }
  out
}

# Position range that covers the variant in VCF form and its changed bases.
interp_variant_span <- function(variant) c(min(variant$pos1, variant$start1), max(variant$pos1 + nchar(variant$ref) - 1L, variant$end1))

# GWAS Catalog associations reported at this position (position-indexed copy,
# see scripts/build_gwas_index.sh). The catalog keys a row by rsID and position,
# not by allele: risk_allele_is says whether the reported risk allele is this
# variant's ALT, its REF, or neither / not reported.
interp_gwas_catalog_empty <- function() data.frame(rs_id = character(), risk_allele = character(), risk_allele_is = character(),
  trait = character(), mapped_trait = character(), p_value = numeric(), effect = numeric(), ci_text = character(),
  risk_allele_frequency = character(), mapped_gene = character(), pubmed = character(), first_author = character(),
  date = character(), study_accession = character(), initial_sample = character(), stringsAsFactors = FALSE)

interp_gwas_catalog <- function(variant, sources = interp_local_sources()) {
  out <- list(available = FALSE, found = FALSE, associations = interp_gwas_catalog_empty(), release = "",
    source = "GWAS Catalog associations (local position-indexed copy)")
  file <- sources$gwas_catalog
  if (!interp_readable(sources$tabix) || !interp_readable(file) || !interp_readable(paste0(file, ".tbi"))) return(out)
  out$available <- TRUE
  out$release <- sub("^.*_(e[0-9]+_r[0-9-]+).*$", "\\1", basename(file))
  if (identical(out$release, basename(file))) out$release <- ""
  span <- interp_variant_span(variant)
  lines <- interp_tabix(sources, file, variant$chrom, span[[1L]], span[[2L]])
  lines <- lines[!startsWith(lines, "#")]
  if (!length(lines)) return(out)
  f <- strsplit(lines, "\t", fixed = TRUE)
  field <- function(i) vapply(f, function(x) if (length(x) >= i) x[[i]] else "", character(1))
  allele <- toupper(sub("^.*-", "", field(4L))); allele[!grepl("-", field(4L), fixed = TRUE)] <- ""
  is <- ifelse(allele == variant$alt | allele == variant$alt_changed, "ALT",
    ifelse(allele == variant$ref | allele == variant$ref_changed, "REF", ifelse(allele %in% c("", "?"), "not reported", "another allele")))
  x <- data.frame(rs_id = field(3L), risk_allele = field(4L), risk_allele_is = is, trait = field(5L), mapped_trait = field(6L),
    p_value = suppressWarnings(as.numeric(field(7L))), effect = suppressWarnings(as.numeric(field(9L))), ci_text = field(10L),
    risk_allele_frequency = field(11L), mapped_gene = field(12L), pubmed = field(13L), first_author = field(14L), date = field(15L),
    study_accession = field(16L), initial_sample = field(17L), stringsAsFactors = FALSE)
  # A p-value too small for a double is kept as 0 and ordered by its stored -log10.
  mlog <- suppressWarnings(as.numeric(field(8L)))
  x <- x[order(-ifelse(is.finite(mlog), mlog, -Inf), x$trait), , drop = FALSE]; rownames(x) <- NULL
  out$associations <- x; out$found <- TRUE
  out
}

# ClinVar record of this exact allele from a tabix-indexed ClinVar VCF.
interp_clinvar_empty <- function() data.frame(variation_id = character(), clinical_significance = character(),
  review_status = character(), conditions = character(), gene = character(), molecular_consequence = character(),
  rs_id = character(), stringsAsFactors = FALSE)

interp_clinvar_local <- function(variant, sources = interp_local_sources()) {
  out <- list(available = FALSE, found = FALSE, records = interp_clinvar_empty(), other_alleles = 0L, file_date = "",
    source = "ClinVar VCF (local)")
  file <- sources$clinvar_vcf
  if (!interp_readable(sources$tabix) || !interp_readable(file) || !interp_readable(paste0(file, ".tbi"))) return(out)
  out$available <- TRUE
  header <- tryCatch(interp_run_lines(sources$tabix, c("-H", file)), error = function(e) character())
  date <- sub("^##fileDate=", "", grep("^##fileDate=", header, value = TRUE))
  if (length(date)) { out$file_date <- date[[1L]]; out$source <- paste0("ClinVar VCF ", date[[1L]], " (local)") }
  span <- interp_variant_span(variant)
  lines <- interp_tabix(sources, file, interp_bare_chrom(variant$chrom), span[[1L]], span[[2L]])
  rows <- list()
  for (line in lines) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(f) < 8L) next
    if (!identical(f[[2L]], as.character(variant$pos1)) || !identical(f[[4L]], variant$ref) || !identical(f[[5L]], variant$alt)) {
      if (!identical(f[[5L]], ".")) out$other_alleles <- out$other_alleles + 1L
      next
    }
    info <- strsplit(f[[8L]], ";", fixed = TRUE)[[1L]]
    kv <- stats::setNames(sub("^[^=]*=", "", info), sub("=.*$", "", info))
    get <- function(key) if (key %in% names(kv)) utils::URLdecode(gsub("_", " ", kv[[key]])) else ""
    rows[[length(rows) + 1L]] <- data.frame(variation_id = f[[3L]], clinical_significance = get("CLNSIG"),
      review_status = get("CLNREVSTAT"), conditions = gsub("|", "; ", get("CLNDN"), fixed = TRUE),
      gene = sub(":.*$", "", get("GENEINFO")), molecular_consequence = gsub("SO:[0-9]+\\|", "", get("MC")),
      rs_id = if (nzchar(get("RS"))) paste0("rs", get("RS")) else "", stringsAsFactors = FALSE)
  }
  if (length(rows)) { out$records <- do.call(rbind, rows); out$found <- TRUE }
  out
}

# One variant's rows from a CSV sorted on its first column (the SNP id).
interp_look_rows <- function(sources, file, variant) {
  header <- strsplit(readLines(file, n = 1L, warn = FALSE), ",", fixed = TRUE)[[1L]]
  prefix <- paste0(variant$chrom, ":", variant$pos1, ":", variant$ref, ":", variant$alt, ":")
  lines <- interp_run_lines(sources$look, c(prefix, file), env = c("current", LC_ALL = "C"))
  if (!length(lines)) return(NULL)
  x <- utils::read.csv(text = paste(lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE,
    col.names = header, check.names = FALSE)
  x[startsWith(x[[1L]], prefix), , drop = FALSE]
}

interp_fhs_empty <- function() data.frame(set = character(), feature = character(), gene = character(), cis = logical(),
  beta_alt = numeric(), t_value = numeric(), r_squared = numeric(), log10_p = numeric(), alt_frequency = numeric(), n = integer(),
  stringsAsFactors = FALSE)

interp_fhs_rows <- function(x, set) {
  if (is.null(x) || !nrow(x)) return(interp_fhs_empty())
  eaf <- if ("EAF" %in% names(x)) x$EAF else if ("MAF" %in% names(x)) x$MAF else NA_real_
  data.frame(set = set, feature = x$Tx, gene = x$Tx_Gene, cis = x$Is_Cis == 1, beta_alt = x$Fx, t_value = x$T,
    r_squared = x$RSq, log10_p = x$log10P, alt_frequency = eaf, n = as.integer(x$N), stringsAsFactors = FALSE)
}

# FHS whole-blood QTL summary statistics: isoform-ratio sQTL (discovery and
# replication) and gene-level eQTL. Effects are per ALT allele.
interp_fhs_qtl <- function(variant, sources = interp_local_sources()) {
  out <- list(sqtl_available = FALSE, eqtl_available = FALSE, sqtl = interp_fhs_empty(), eqtl = interp_fhs_empty(),
    sqtl_note = "Isoform ratio (transcript / gene) in FHS whole blood; discovery N about 2,622, replication about 1,093; rows with p < 1e-4.",
    eqtl_note = "Gene expression in FHS whole blood, N about 2,622; rows with p < 1e-7.")
  if (!nzchar(sources$look)) return(out)
  if (all(vapply(sources$fhs_sqtl, interp_readable, logical(1)))) {
    out$sqtl_available <- TRUE
    parts <- lapply(names(sources$fhs_sqtl), function(set) interp_fhs_rows(interp_look_rows(sources, sources$fhs_sqtl[[set]], variant), set))
    out$sqtl <- do.call(rbind, parts)
    if (nrow(out$sqtl)) out$sqtl <- out$sqtl[order(!out$sqtl$cis, out$sqtl$log10_p), , drop = FALSE]
    rownames(out$sqtl) <- NULL
  }
  if (interp_readable(sources$fhs_eqtl)) {
    out$eqtl_available <- TRUE
    out$eqtl <- interp_fhs_rows(interp_look_rows(sources, sources$fhs_eqtl, variant), "combined")
    if (nrow(out$eqtl)) out$eqtl <- out$eqtl[order(!out$eqtl$cis, out$eqtl$log10_p), , drop = FALSE]
    rownames(out$eqtl) <- NULL
  }
  out
}

# Isoform-ratio sQTL rows that replicate: same transcript in both sets with the
# same direction of effect.
interp_fhs_replicated <- function(sqtl) {
  d <- sqtl[sqtl$set == "discovery", , drop = FALSE]; r <- sqtl[sqtl$set == "replication", , drop = FALSE]
  m <- merge(d, r, by = "feature", suffixes = c("_discovery", "_replication"))
  m <- m[sign(m$beta_alt_discovery) == sign(m$beta_alt_replication), , drop = FALSE]
  if (!nrow(m)) return(data.frame(feature = character(), gene = character(), cis = logical(), beta_discovery = numeric(),
    log10_p_discovery = numeric(), beta_replication = numeric(), log10_p_replication = numeric(), stringsAsFactors = FALSE))
  out <- data.frame(feature = m$feature, gene = m$gene_discovery, cis = m$cis_discovery, beta_discovery = m$beta_alt_discovery,
    log10_p_discovery = m$log10_p_discovery, beta_replication = m$beta_alt_replication, log10_p_replication = m$log10_p_replication,
    stringsAsFactors = FALSE)
  out <- out[order(out$log10_p_discovery), , drop = FALSE]; rownames(out) <- NULL; out
}

# -- TOPMed multi-tissue cis-QTL ---------------------------------------------
interp_topmed_contexts <- function() c(Lung = "Lung", Monocyte = "Monocytes", Nasal_epithelial = "Nasal epithelium",
  PBMC = "PBMC", T_cell = "T cells", Whole_blood = "Whole blood")

interp_topmed_empty <- function() data.frame(qtl_type = character(), context = character(), phenotype = character(),
  gene_id = character(), gene = character(), alt_frequency = numeric(), n_estimated = integer(), slope_alt = numeric(),
  slope_se = numeric(), p_nominal = numeric(), p_threshold = numeric(), significant = logical(), pip = numeric(),
  credible_set = character(), tss_distance = integer(), stringsAsFactors = FALSE)

# Every tested pair of this variant in each context. `significant` uses the
# phenotype's own permutation threshold (tensorQTL, FDR 5% across phenotypes);
# `pip` is the SuSiE posterior inclusion probability where the variant is in a
# credible set. Effects are per ALT allele on normalised phenotypes.
interp_topmed_qtl <- function(variant, sources = interp_local_sources(), keep_p = 1e-3, keep_top = 3L) {
  out <- list(available = FALSE, pairs = interp_topmed_empty(), contexts = data.frame(qtl_type = character(), context = character(),
    tested = logical(), pairs_tested = integer(), significant_pairs = integer(), min_p = numeric(), stringsAsFactors = FALSE))
  root <- sources$topmed_qtl
  if (!interp_readable(sources$tabix) || !dir.exists(file.path(root, "cis-eqtl", "nominal"))) return(out)
  out$available <- TRUE
  id <- paste(variant$chrom, variant$pos1, variant$ref, variant$alt, sep = "_")
  layout <- list(eQTL = list(dir = "cis-eqtl", cs = ".cs.txt"), sQTL = list(dir = "cis-sqtl2", cs = ".cs-per-intron.txt"))
  pairs <- list(); summary <- list()
  for (type in names(layout)) for (context in names(interp_topmed_contexts())) {
    base <- file.path(root, layout[[type]]$dir)
    # Whole blood also has a run down to MAF 0.1%; use it when the 1% run has no row.
    sets <- if (context == "Whole_blood") c("maf001", "maf0001") else "maf001"
    rows <- NULL; used <- sets[[1L]]
    for (set in sets) {
      file <- file.path(base, "nominal", paste0(context, ".", set, ".tsv.gz"))
      if (!interp_readable(file) || !interp_readable(paste0(file, ".tbi"))) next
      lines <- interp_tabix(sources, file, variant$chrom, variant$pos1)
      lines <- lines[vapply(strsplit(lines, "\t", fixed = TRUE), function(f) length(f) >= 12L && f[[5L]] == id, logical(1))]
      if (length(lines)) { rows <- do.call(rbind, strsplit(lines, "\t", fixed = TRUE)); used <- set; break }
    }
    label <- unname(interp_topmed_contexts()[[context]])
    if (is.null(rows)) {
      summary[[length(summary) + 1L]] <- data.frame(qtl_type = type, context = label, tested = FALSE, pairs_tested = 0L,
        significant_pairs = 0L, min_p = NA_real_, stringsAsFactors = FALSE)
      next
    }
    x <- data.frame(qtl_type = type, context = label, phenotype = rows[, 4L], tss_distance = as.integer(rows[, 6L]),
      alt_frequency = as.numeric(rows[, 7L]), ma_count = as.numeric(rows[, 9L]), p_nominal = as.numeric(rows[, 10L]),
      slope_alt = as.numeric(rows[, 11L]), slope_se = as.numeric(rows[, 12L]), stringsAsFactors = FALSE)
    x$gene_id <- sub("^.*:", "", x$phenotype)
    minor <- pmin(x$alt_frequency, 1 - x$alt_frequency)
    x$n_estimated <- as.integer(round(x$ma_count / (2 * minor)))
    x$p_threshold <- NA_real_
    hits <- Sys.glob(file.path(base, "top-hits", paste0(context, ".*", used, ".permutations.txt.gz")))
    if (!length(hits)) hits <- Sys.glob(file.path(base, "top-hits", paste0(context, ".", used, ".permutations.txt.gz")))
    if (length(hits)) {
      top <- data.table::fread(hits[[1L]], select = c("phenotype_id", "qval", "pval_nominal_threshold"), showProgress = FALSE)
      top <- top[top$qval <= 0.05, ]
      # eQTL thresholds are per gene; sQTL thresholds are per gene (grouped introns).
      key <- if (type == "eQTL") top$phenotype_id else sub("^.*:", "", top$phenotype_id)
      x$p_threshold <- top$pval_nominal_threshold[match(x$gene_id, key)]
    }
    x$significant <- !is.na(x$p_threshold) & x$p_nominal <= x$p_threshold
    x$pip <- NA_real_; x$credible_set <- ""
    cs <- file.path(base, "susie", paste0(context, ".", used, layout[[type]]$cs))
    if (interp_readable(cs)) {
      found <- interp_run_lines(unname(Sys.which("grep")), c("-F", paste0("\t", id, "\t"), cs))
      for (line in found) {
        f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
        k <- which(x$phenotype == f[[1L]])
        if (length(k) && length(f) >= 5L) { x$pip[k] <- as.numeric(f[[3L]]); x$credible_set[k] <- f[[5L]] }
      }
    }
    summary[[length(summary) + 1L]] <- data.frame(qtl_type = type, context = label, tested = TRUE, pairs_tested = nrow(x),
      significant_pairs = sum(x$significant), min_p = min(x$p_nominal, na.rm = TRUE), stringsAsFactors = FALSE)
    x <- x[order(x$p_nominal), , drop = FALSE]
    x$top <- seq_len(nrow(x)) <= keep_top
    x$gene <- ""
    pairs[[length(pairs) + 1L]] <- x
  }
  out$contexts <- do.call(rbind, summary)
  if (length(pairs)) {
    all <- do.call(rbind, pairs)
    # Keep every context's row of a phenotype that is significant somewhere, so
    # the comparison across contexts sees the non-significant ones too.
    key <- paste(all$qtl_type, interp_phenotype_key(all$qtl_type, all$phenotype))
    keep <- all$significant | !is.na(all$pip) | all$p_nominal <= keep_p | all$top | key %in% key[all$significant]
    out$pairs <- all[keep, names(interp_topmed_empty()), drop = FALSE]; rownames(out$pairs) <- NULL
  }
  out
}

# Same phenotype across contexts: intron coordinates plus gene for sQTL (cluster
# numbers differ between runs), gene for eQTL.
interp_phenotype_key <- function(qtl_type, phenotype)
  ifelse(qtl_type == "sQTL", sub("^([^:]+:[0-9]+:[0-9]+):.*:(ENSG[^:]*)$", "\\1:\\2", phenotype), phenotype)

# Is the effect shared across contexts or specific to some? Fixed-effect
# heterogeneity (Cochran's Q, I-squared) of the per-allele slopes, for every
# phenotype tested in at least two contexts. Descriptive: the contexts share
# donors, so Q is not an exact test.
interp_context_heterogeneity <- function(pairs) {
  empty <- data.frame(qtl_type = character(), phenotype = character(), gene_id = character(), gene = character(),
    contexts_tested = integer(), contexts_significant = character(), strongest_context = character(), min_p = numeric(),
    pooled_slope = numeric(), q_p_value = numeric(), i_squared = numeric(), pattern = character(), stringsAsFactors = FALSE)
  if (is.null(pairs) || !nrow(pairs)) return(empty)
  pairs$key <- interp_phenotype_key(pairs$qtl_type, pairs$phenotype)
  rows <- lapply(split(pairs, paste(pairs$qtl_type, pairs$key)), function(x) {
    x <- x[is.finite(x$slope_alt) & is.finite(x$slope_se) & x$slope_se > 0, , drop = FALSE]
    if (nrow(x) < 1L || !any(x$significant)) return(NULL)
    w <- 1 / x$slope_se^2; pooled <- sum(w * x$slope_alt) / sum(w)
    q <- sum(w * (x$slope_alt - pooled)^2); df <- nrow(x) - 1L
    sig <- x$context[x$significant]
    same_sign <- all(sign(x$slope_alt[x$significant]) == sign(x$slope_alt[x$significant][[1L]]))
    q_p <- if (df > 0L) stats::pchisq(q, df, lower.tail = FALSE) else NA_real_
    pattern <- if (nrow(x) == 1L) "tested in one context only" else if (!same_sign) "opposite directions between contexts" else
      if (length(sig) == nrow(x)) "shared: significant in every tested context" else
      if (q_p >= 0.05) paste0("consistent effect size; significant in ", length(sig), " of ", nrow(x), " (the others may lack power)") else
      if (length(sig) == 1L) paste0("context-specific: significant only in ", sig, ", effect differs elsewhere") else
      paste0("heterogeneous: strongest in ", x$context[[which.max(abs(x$slope_alt))]])
    data.frame(qtl_type = x$qtl_type[[1L]], phenotype = x$key[[1L]], gene_id = x$gene_id[[1L]], gene = x$gene[[1L]],
      contexts_tested = nrow(x), contexts_significant = paste(sig, collapse = ", "),
      strongest_context = x$context[[which.min(x$p_nominal)]], min_p = min(x$p_nominal), pooled_slope = pooled,
      q_p_value = q_p,
      i_squared = if (df > 0L && q > 0) max(0, (q - df) / q) else NA_real_, pattern = pattern, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows); out <- out[order(out$min_p), , drop = FALSE]; rownames(out) <- NULL; out
}

# -- SpliceAI (run on request) -----------------------------------------------
interp_spliceai_empty <- function() data.frame(gene = character(), acceptor_gain = numeric(), acceptor_loss = numeric(),
  donor_gain = numeric(), donor_loss = numeric(), acceptor_gain_pos1 = integer(), acceptor_loss_pos1 = integer(),
  donor_gain_pos1 = integer(), donor_loss_pos1 = integer(), max_delta = numeric(), stringsAsFactors = FALSE)

interp_parse_spliceai <- function(vcf_lines, variant) {
  out <- interp_spliceai_empty()
  data <- vcf_lines[!startsWith(vcf_lines, "#")]
  if (!length(data)) return(out)
  info <- strsplit(data[[1L]], "\t", fixed = TRUE)[[1L]]
  if (length(info) < 8L || !grepl("SpliceAI=", info[[8L]], fixed = TRUE)) return(out)
  value <- sub("^.*SpliceAI=([^;]*).*$", "\\1", info[[8L]])
  rows <- lapply(strsplit(value, ",", fixed = TRUE)[[1L]], function(entry) {
    f <- strsplit(entry, "|", fixed = TRUE)[[1L]]
    if (length(f) < 10L) return(NULL)
    score <- suppressWarnings(as.numeric(f[3:6])); offset <- suppressWarnings(as.integer(f[7:10]))
    data.frame(gene = f[[2L]], acceptor_gain = score[[1L]], acceptor_loss = score[[2L]], donor_gain = score[[3L]],
      donor_loss = score[[4L]], acceptor_gain_pos1 = variant$pos1 + offset[[1L]], acceptor_loss_pos1 = variant$pos1 + offset[[2L]],
      donor_gain_pos1 = variant$pos1 + offset[[3L]], donor_loss_pos1 = variant$pos1 + offset[[4L]],
      max_delta = suppressWarnings(max(score, na.rm = TRUE)), stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows)) do.call(rbind, rows) else out
}

# Runs the SCC SpliceAI module on a one-line VCF. About a minute on CPU, most of
# it loading the models. Scores cover the genes in SpliceAI's own annotation.
interp_spliceai <- function(variant, fasta, sources = interp_local_sources(), distance = 500L, timeout = 900) {
  if (!interp_readable(sources$spliceai)) stop("SpliceAI is not available on this system.")
  if (!interp_readable(fasta)) stop("SpliceAI needs the reference FASTA of the connected dataset.")
  if (nchar(variant$ref) > 1L && nchar(variant$alt) > 1L) stop("SpliceAI scores SNVs, insertions and deletions only.")
  work <- tempfile("regshiny_spliceai_"); dir.create(work, mode = "0700"); on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "in.vcf"); output <- file.path(work, "out.vcf")
  writeLines(c("##fileformat=VCFv4.2", paste0("##contig=<ID=", variant$chrom, ">"),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    paste(variant$chrom, variant$pos1, ".", variant$ref, variant$alt, ".", ".", ".", sep = "\t")), input)
  p <- processx::run(sources$spliceai, c("-I", input, "-O", output, "-R", fasta, "-A", "grch38", "-D", as.character(distance)),
    timeout = timeout, error_on_status = FALSE)
  if (isTRUE(p$timeout)) stop("SpliceAI timed out.")
  if (p$status != 0L || !file.exists(output)) stop("SpliceAI failed: ", trimws(tail(strsplit(p$stderr, "\n")[[1L]], 1L)))
  list(scores = interp_parse_spliceai(readLines(output, warn = FALSE), variant), distance = distance,
    model = "SpliceAI 1.3.1, GRCh38, masked = no", ran_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

# One row per predicted change of at least `min_delta`, with the coordinate the
# app uses for a splice site (the intronic base next to it) and whether the
# selected GENCODE release annotates a site there. SpliceAI reports the exonic
# boundary base; strand is the gene's.
interp_spliceai_sites <- function(scores, strand, introns = NULL, min_delta = 0.1) {
  empty <- data.frame(gene = character(), effect = character(), delta = numeric(), exon_boundary_pos1 = integer(),
    site_type = character(), site_position = integer(), annotated_site = character(), stringsAsFactors = FALSE)
  if (is.null(scores) || !nrow(scores)) return(empty)
  kinds <- list(acceptor_gain = c("acceptor", "gain"), acceptor_loss = c("acceptor", "loss"),
    donor_gain = c("donor", "gain"), donor_loss = c("donor", "loss"))
  rows <- list()
  for (i in seq_len(nrow(scores))) for (k in names(kinds)) {
    delta <- scores[[k]][[i]]
    if (is.na(delta) || delta < min_delta) next
    pos <- scores[[paste0(k, "_pos1")]][[i]]; type <- kinds[[k]][[1L]]; s <- strand[[i]]
    # Donor: the intron follows the exon end; acceptor: it precedes the exon start.
    step <- if (is.na(s)) NA_integer_ else if ((type == "donor") == (s == "+")) 1L else -1L
    site <- pos + step
    annotated <- if (is.null(introns) || is.na(site)) "unknown" else {
      known <- if (type == "donor") introns$donor1[introns$strand == s] else introns$acceptor1[introns$strand == s]
      if (site %in% known) "yes" else "no"
    }
    rows[[length(rows) + 1L]] <- data.frame(gene = scores$gene[[i]], effect = paste(type, kinds[[k]][[2L]]), delta = delta,
      exon_boundary_pos1 = as.integer(pos), site_type = type, site_position = as.integer(site), annotated_site = annotated,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows); out <- out[order(-out$delta), , drop = FALSE]; rownames(out) <- NULL; out
}

# ---- Genes linked to the variant, and their pathways -------------------------
# Genes the variant can act through: the ones it overlaps, and the ones whose
# splicing or expression is associated with it in the local QTL sources.
interp_linked_genes <- function(consequences, local = NULL) {
  rows <- list()
  add <- function(gene_id, gene, relation) if (length(gene_id)) rows[[length(rows) + 1L]] <<- data.frame(
    gene_id = interp_bare_gene(gene_id), gene = gene, relation = relation, stringsAsFactors = FALSE)
  genic <- consequences[!consequences$most_severe %in% c("upstream_gene_variant", "downstream_gene_variant"), , drop = FALSE]
  add(genic$gene_id, genic$gene_name, "overlaps")
  topmed <- if (interp_ok(local) && interp_ok(local$topmed) && isTRUE(local$topmed$available)) local$topmed$pairs else NULL
  if (!is.null(topmed) && nrow(topmed)) {
    hit <- topmed[topmed$significant %in% TRUE, , drop = FALSE]
    add(hit$gene_id[hit$qtl_type == "sQTL"], hit$gene[hit$qtl_type == "sQTL"], "sQTL target (TOPMed)")
    add(hit$gene_id[hit$qtl_type == "eQTL"], hit$gene[hit$qtl_type == "eQTL"], "eQTL target (TOPMed)")
  }
  if (!length(rows)) return(data.frame(gene_id = character(), gene = character(), relation = character(), stringsAsFactors = FALSE))
  x <- unique(do.call(rbind, rows))
  # One row per gene; a name from any row, relations joined in a fixed order.
  out <- do.call(rbind, lapply(split(x, x$gene_id), function(g) data.frame(gene_id = g$gene_id[[1L]],
    gene = { n <- g$gene[nzchar(g$gene) & !startsWith(g$gene, "ENSG")]; if (length(n)) n[[1L]] else "" },
    relation = paste(unique(g$relation), collapse = "; "), stringsAsFactors = FALSE)))
  out <- out[order(!grepl("overlaps", out$relation), out$gene, out$gene_id), , drop = FALSE]; rownames(out) <- NULL; out
}

# Reactome's top-level pathways. The mapping file lists every level, so a gene's
# rows with one of these names are its top-level pathways.
interp_reactome_top_levels <- function() c("Autophagy", "Cell Cycle", "Cell-Cell communication", "Cellular responses to stimuli",
  "Chromatin organization", "Circadian clock", "Developmental Biology", "Digestion and absorption", "Disease", "DNA Repair",
  "DNA Replication", "Drug ADME", "Extracellular matrix organization", "Gene expression (Transcription)", "Hemostasis",
  "Immune System", "Metabolism", "Metabolism of proteins", "Metabolism of RNA", "Muscle contraction", "Neuronal System",
  "Organelle biogenesis and maintenance", "Programmed Cell Death", "Protein localization", "Reproduction", "Sensory Perception",
  "Signal Transduction", "Transport of small molecules", "Vesicle-mediated transport")

interp_msigdb_files <- function(sources) {
  if (length(sources$msigdb) != 1L || !nzchar(sources$msigdb) || !dir.exists(sources$msigdb)) return(character())
  pick <- function(pattern) { f <- Sys.glob(file.path(sources$msigdb, pattern)); if (length(f)) f[[1L]] else NA_character_ }
  files <- c(Hallmark = pick("h.all.v*.Hs.symbols.gmt"), `Canonical pathways` = pick("c2.cp.v*.Hs.symbols.gmt"),
    `GO biological process` = pick("c5.go.bp.v*.Hs.symbols.gmt"))
  files[!is.na(files) & vapply(files, interp_readable, logical(1))]
}

interp_pathways_empty <- function() list(
  reactome = data.frame(gene = character(), relation = character(), pathway = character(), pathway_id = character(),
    top_level = logical(), evidence = character(), stringsAsFactors = FALSE),
  gene_sets = data.frame(gene = character(), relation = character(), collection = character(), gene_set = character(),
    genes_in_set = integer(), stringsAsFactors = FALSE))

# Pathways and gene sets of the linked genes, from local files only. They
# describe the genes; they are not evidence that the variant changes a pathway.
interp_local_pathways <- function(genes, sources = interp_local_sources()) {
  out <- c(list(available = FALSE, reactome_available = FALSE, msigdb_available = FALSE, genes = genes), interp_pathways_empty(),
    list(source = character()))
  if (is.null(genes) || !nrow(genes)) return(out)
  label <- ifelse(nzchar(genes$gene), genes$gene, genes$gene_id)
  if (interp_readable(sources$look) && interp_readable(sources$reactome)) {
    out$reactome_available <- TRUE; out$source <- c(out$source, paste0("Reactome (", basename(sources$reactome), ")"))
    rows <- list()
    for (i in seq_len(nrow(genes))) {
      lines <- interp_run_lines(sources$look, c(paste0(genes$gene_id[[i]], "\t"), sources$reactome), env = c("current", LC_ALL = "C"))
      if (!length(lines)) next
      f <- strsplit(lines, "\t", fixed = TRUE)
      x <- data.frame(gene = label[[i]], relation = genes$relation[[i]], pathway = vapply(f, `[[`, character(1), 4L),
        pathway_id = vapply(f, `[[`, character(1), 2L), evidence = vapply(f, `[[`, character(1), 5L), stringsAsFactors = FALSE)
      # The same pathway can be listed once per evidence code.
      x <- do.call(rbind, lapply(split(x, x$pathway_id), function(g) { g$evidence <- paste(sort(unique(g$evidence)), collapse = ", "); g[1L, ] }))
      x$top_level <- x$pathway %in% interp_reactome_top_levels()
      rows[[length(rows) + 1L]] <- x[order(!x$top_level, x$pathway), c("gene", "relation", "pathway", "pathway_id", "top_level", "evidence")]
    }
    if (length(rows)) { out$reactome <- do.call(rbind, rows); rownames(out$reactome) <- NULL }
  }
  files <- interp_msigdb_files(sources)
  named <- nzchar(genes$gene)
  if (length(files)) {
    out$msigdb_available <- TRUE
    out$source <- c(out$source, paste0("MSigDB ", sub("^.*\\.v([0-9.]+)\\.Hs\\..*$", "\\1", basename(files[[1L]]))))
    rows <- list()
    if (any(named)) for (collection in names(files)) {
      sets <- strsplit(readLines(files[[collection]], warn = FALSE), "\t", fixed = TRUE)
      # Reactome sets are left to the Reactome mapping above, which is keyed by gene ID.
      sets <- sets[!startsWith(vapply(sets, `[[`, character(1), 1L), "REACTOME_")]
      for (i in which(named)) {
        hit <- vapply(sets, function(set) genes$gene[[i]] %in% set[-(1:2)], logical(1))
        if (any(hit)) rows[[length(rows) + 1L]] <- data.frame(gene = genes$gene[[i]], relation = genes$relation[[i]], collection = collection,
          gene_set = vapply(sets[hit], `[[`, character(1), 1L), genes_in_set = lengths(sets[hit]) - 2L, stringsAsFactors = FALSE)
      }
    }
    if (length(rows)) {
      x <- do.call(rbind, rows)
      out$gene_sets <- x[order(match(x$gene, genes$gene), match(x$collection, names(files)), x$genes_in_set), , drop = FALSE]
      rownames(out$gene_sets) <- NULL
    }
  }
  out$available <- out$reactome_available || out$msigdb_available
  out
}

# All local lookups; each source fails on its own.
interp_local_evidence <- function(variant, sources = interp_local_sources()) {
  attempt <- function(expr) tryCatch(expr, error = function(e) e)
  list(cadd = attempt(interp_cadd(variant, sources)), known = attempt(interp_known_variant(variant, sources)),
    fhs = attempt(interp_fhs_qtl(variant, sources)), topmed = attempt(interp_topmed_qtl(variant, sources)),
    gwas_catalog = attempt(interp_gwas_catalog(variant, sources)), clinvar = attempt(interp_clinvar_local(variant, sources)),
    status = interp_local_status(sources))
}

# ---- Evidence summary ------------------------------------------------------
# Fixed rules that turn the tables above into one line per question. Levels:
#   Strong      direct or fine-mapped evidence that this variant is responsible
#   Moderate    an association or a prediction, but not both
#   Limited     weak or indirect support only
#   None found  the sources were consulted and report nothing
#   Not assessed  the sources needed were unavailable or not queried
# The rules are listed in docs/METRICS.md. They rank evidence; they are not a
# clinical classification.
interp_levels <- function() c("Strong", "Moderate", "Limited", "None found", "Not assessed")

interp_ok <- function(x) !is.null(x) && !inherits(x, "error")

interp_fmt_p <- function(p) vapply(p, function(x) if (is.na(x)) "NA" else if (x < 1e-300) "< 1e-300" else
  formatC(x, format = "e", digits = 1), character(1))

# Cohort RNA evidence from the open result's Model 1 tables (view-time values).
interp_rna_evidence <- function(bin_ratio = NULL, junctions = NULL) {
  out <- list(available = FALSE, sites = character(), junctions = character(), supported = FALSE, details = character())
  ok_rows <- function(x) is.data.frame(x) && nrow(x) && all(c("status", "p_permutation", "beta") %in% names(x))
  if (ok_rows(bin_ratio)) {
    out$available <- TRUE
    x <- bin_ratio[bin_ratio$status == "OK" & is.finite(bin_ratio$p_permutation), , drop = FALSE]
    if (nrow(x)) {
      flat <- if ("bin_check" %in% names(x)) grepl("NO_DEPTH_DROP_AT_SITE", x$bin_check, fixed = TRUE) else rep(FALSE, nrow(x))
      hit <- x$p_permutation <= 0.05 / nrow(x) & !flat
      if (any(hit)) {
        h <- x[hit, , drop = FALSE]; h <- h[order(h$p_permutation, -abs(h$beta)), , drop = FALSE]
        out$sites <- h$site_id
        out$details <- c(out$details, sprintf("bin ratio at %s: %+.3f per ALT allele (permutation p %s, n = %d)",
          h$site_id, h$beta, interp_fmt_p(h$p_permutation), h$n))
      }
      if (any(flat)) out$details <- c(out$details, paste0(sum(flat), " site(s) set aside: no depth drop at the site itself"))
    }
  }
  if (ok_rows(junctions)) {
    out$available <- TRUE
    x <- junctions[junctions$status == "OK" & is.finite(junctions$p_permutation), , drop = FALSE]
    if (nrow(x)) {
      hit <- x$p_permutation <= 0.05 / nrow(x)
      if (any(hit)) {
        h <- x[hit, , drop = FALSE]; h <- h[order(h$p_permutation, -abs(h$beta)), , drop = FALSE]
        id <- paste0(h$chrom, ":", h$intron_start1, "-", h$intron_end1)
        out$junctions <- id
        out$details <- c(out$details, sprintf("junction %s: %+.3f reads per 100 per ALT allele (permutation p %s, n = %d)",
          id, h$beta, interp_fmt_p(h$p_permutation), h$n))
      }
    }
  }
  out$supported <- length(out$sites) > 0L || length(out$junctions) > 0L
  out
}

interp_summary_row <- function(question, level, finding, evidence = character(), sources = character())
  data.frame(question = question, level = level, finding = finding, evidence = paste(evidence[nzchar(evidence)], collapse = " | "),
    sources = paste(unique(sources[nzchar(sources)]), collapse = "; "), stringsAsFactors = FALSE)

interp_evidence_summary <- function(consequences, local = NULL, public = NULL, spliceai = NULL, rna = NULL) {
  rows <- list()
  genic <- consequences[!consequences$most_severe %in% c("upstream_gene_variant", "downstream_gene_variant"), , drop = FALSE]
  ot <- if (interp_ok(public) && interp_ok(public$open_targets) && isTRUE(public$open_targets$found)) public$open_targets else NULL
  ot_queried <- interp_ok(public) && interp_ok(public$open_targets)
  gtex <- if (interp_ok(public) && interp_ok(public$gtex)) public$gtex else NULL
  myv <- if (interp_ok(public) && interp_ok(public$myvariant) && isTRUE(public$myvariant$found)) public$myvariant else NULL
  cadd <- if (interp_ok(local) && interp_ok(local$cadd) && isTRUE(local$cadd$found)) local$cadd$phred else NA_real_
  known <- if (interp_ok(local) && interp_ok(local$known) && isTRUE(local$known$found)) local$known else NULL
  topmed <- if (interp_ok(local) && interp_ok(local$topmed) && isTRUE(local$topmed$available)) local$topmed else NULL
  fhs <- if (interp_ok(local) && interp_ok(local$fhs)) local$fhs else NULL

  # 1. Protein sequence.
  coding_terms <- c("stop_gained", "frameshift_variant", "stop_lost", "start_lost", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_acceptor_variant", "splice_donor_variant")
  src <- "GENCODE + reference (computed here)"
  if (!nrow(genic)) {
    rows[[1L]] <- interp_summary_row("Protein sequence", "None found", "The variant lies outside every annotated transcript.", sources = src)
  } else {
    top <- genic[1L, ]
    reference <- genic[grepl("MANE_Select|Ensembl_canonical", genic$transcript_flags), , drop = FALSE]
    ref_top <- if (nrow(reference)) reference[1L, ] else top
    n_severe <- sum(genic$most_severe == top$most_severe & genic$gene_id == top$gene_id)
    # Same denominator as the per-gene table: every transcript of the gene within reach.
    n_gene <- sum(consequences$gene_id == top$gene_id)
    where <- sprintf("%s in %s, %d of %d transcripts", top$most_severe, top$gene_name, n_severe, n_gene)
    hgvs <- trimws(paste(ref_top$hgvs_c, ref_top$hgvs_p))
    evidence <- c(where, if (nzchar(hgvs)) paste0("reference transcript ", ref_top$transcript_id, ": ", ref_top$consequence, " ", hgvs) else "",
      if (is.finite(cadd)) sprintf("CADD PHRED %.1f", cadd) else "")
    predictors <- if (!is.null(ot)) ot$predictors else NULL
    am <- if (!is.null(predictors)) predictors[predictors$method == "AlphaMissense", , drop = FALSE] else NULL
    if (!is.null(am) && nrow(am)) evidence <- c(evidence, sprintf("AlphaMissense %s (%.3f)", am$assessment[[1L]], am$score[[1L]]))
    revel <- if (!is.null(myv)) myv$revel else NA_real_
    if (is.finite(revel)) evidence <- c(evidence, sprintf("REVEL %.2f", revel))
    # An essential splice site of an intron in a UTR or a non-coding transcript
    # (c.-N+1, c.*N+1, n.N+1) does not by itself change the protein sequence.
    outside_cds <- function(x) x$most_severe %in% c("splice_donor_variant", "splice_acceptor_variant") & grepl("^(c\\.[-*]|n\\.)", x$hgvs_c)
    high <- genic[genic$impact == "HIGH" & genic$gene_id == top$gene_id, , drop = FALSE]
    if (top$impact == "HIGH" && all(outside_cds(high))) {
      level <- "Limited"
      finding <- paste0("Essential splice site of an intron outside the coding sequence (", nrow(high), " transcript(s) of ", top$gene_name,
        "); the protein sequence is not changed directly. See Splicing.")
    } else if (top$impact == "HIGH") {
      on_ref <- ref_top$impact == "HIGH" && !outside_cds(ref_top)
      level <- if (on_ref) "Strong" else "Moderate"
      finding <- paste0("Predicted loss of the normal product: ", high$most_severe[!outside_cds(high)][[1L]], " in ", top$gene_name,
        if (on_ref) " on the reference transcript." else " on non-reference transcripts only.")
    } else if (top$most_severe %in% coding_terms) {
      damaging <- (!is.null(am) && nrow(am) && identical(am$assessment[[1L]], "likely_pathogenic")) ||
        (is.finite(revel) && revel >= 0.75) || (is.finite(cadd) && cadd >= 25)
      tolerated <- (!is.null(am) && nrow(am) && identical(am$assessment[[1L]], "likely_benign")) && (!is.finite(cadd) || cadd < 20)
      level <- if (tolerated) "Limited" else "Moderate"
      finding <- paste0("Changes the protein of ", top$gene_name, " (", trimws(paste(top$most_severe, top$hgvs_p)), ")",
        if (damaging) "; predictors call it damaging." else if (tolerated) "; predictors call it tolerated." else ".")
    } else {
      level <- "None found"
      finding <- paste0("No change to the protein sequence is predicted (", top$most_severe, " in ", top$gene_name, ").")
    }
    rows[[1L]] <- interp_summary_row("Protein sequence", level, finding, evidence,
      c(src, if (is.finite(cadd)) "CADD v1.7", if (!is.null(am) && nrow(am)) "Open Targets", if (is.finite(revel)) "MyVariant.info"))
  }

  # 2. Protein abundance (pQTL credible sets).
  if (!ot_queried) {
    rows[[2L]] <- interp_summary_row("Protein abundance (pQTL)", "Not assessed", "Public databases were not queried; no pQTL results are stored locally.")
  } else {
    pq <- if (!is.null(ot)) ot$qtl[ot$qtl$qtl_type == "pqtl", , drop = FALSE] else NULL
    if (is.null(pq) || !nrow(pq)) {
      rows[[2L]] <- interp_summary_row("Protein abundance (pQTL)", "None found", "Not in any pQTL credible set in Open Targets.", sources = "Open Targets")
    } else {
      pip <- ifelse(is.na(pq$posterior_probability), 0, pq$posterior_probability)
      cis <- !(pq$trans %in% TRUE)
      level <- if (any(pip >= 0.9 & cis)) "Strong" else if (any(pip >= 0.1)) "Moderate" else "Limited"
      # Cis sets first: they concern the protein encoded next to the variant.
      order_best <- order(!cis, -pip)
      best <- pq[order_best, , drop = FALSE]; best_pip <- pip[order_best]; best_cis <- cis[order_best]
      rows[[2L]] <- interp_summary_row("Protein abundance (pQTL)", level,
        sprintf("In %d pQTL credible set(s), %d in cis; leading set: %s in %s (%s), posterior probability %.2f.", nrow(pq), sum(cis),
          best$gene[[1L]], best$biosample[[1L]], if (best_cis[[1L]]) "cis" else "trans", best_pip[[1L]]),
        sprintf("%s in %s: PIP %.2f, beta %+.2f, p %s%s", head(best$gene, 6L), head(best$biosample, 6L), head(best_pip, 6L),
          head(best$beta, 6L), interp_fmt_p(head(best$p_value, 6L)), ifelse(head(best_cis, 6L), "", " (trans)")),
        "Open Targets")
    }
  }

  # 3. Splicing.
  evidence <- character(); sources <- character()
  splice_terms <- unique(unlist(strsplit(genic$consequence, "&", fixed = TRUE)))
  essential <- any(c("splice_donor_variant", "splice_acceptor_variant") %in% splice_terms)
  region <- any(c("splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant") %in% splice_terms)
  if (essential || region) {
    hit <- genic[grepl("splice_", genic$consequence), , drop = FALSE]
    evidence <- c(evidence, sprintf("annotation: %s in %d transcript(s) of %s", interp_join(setdiff(splice_terms[grepl("^splice_", splice_terms)], ""), ", "),
      nrow(hit), interp_join(hit$gene_name, ", ")))
    sources <- c(sources, src)
  }
  sai <- if (interp_ok(spliceai) && is.data.frame(spliceai$scores) && nrow(spliceai$scores)) spliceai$scores else NULL
  sai_max <- if (!is.null(sai)) max(sai$max_delta, na.rm = TRUE) else NA_real_
  if (!is.null(sai)) {
    s <- sai[which.max(sai$max_delta), ]
    kinds <- c(acceptor_gain = "acceptor gain", acceptor_loss = "acceptor loss", donor_gain = "donor gain", donor_loss = "donor loss")
    parts <- vapply(names(kinds), function(k) if (s[[k]] >= 0.1) sprintf("%s %.2f at exon base %d", kinds[[k]], s[[k]], s[[paste0(k, "_pos1")]]) else "", character(1))
    evidence <- c(evidence, paste0("SpliceAI (", s$gene, "): ", if (any(nzchar(parts))) paste(parts[nzchar(parts)], collapse = ", ") else
      sprintf("largest change %.2f", sai_max)))
    sources <- c(sources, "SpliceAI")
  }
  observed <- interp_ok(rna) && isTRUE(rna$supported)
  if (observed) { evidence <- c(evidence, paste0("this cohort: ", head(rna$details, 3L))); sources <- c(sources, "Model 1 of this result") }
  s_pairs <- if (!is.null(topmed)) topmed$pairs[topmed$pairs$qtl_type == "sQTL" & topmed$pairs$significant, , drop = FALSE] else NULL
  topmed_sqtl <- !is.null(s_pairs) && nrow(s_pairs) > 0L
  finemapped <- topmed_sqtl && any(s_pairs$pip >= 0.5, na.rm = TRUE)
  if (topmed_sqtl) {
    best <- s_pairs[which.min(s_pairs$p_nominal), ]
    evidence <- c(evidence, sprintf("TOPMed sQTL significant in %d of %d contexts (%s); strongest %s, slope %+.2f, p %s%s",
      length(unique(s_pairs$context)), sum(topmed$contexts$qtl_type == "sQTL" & topmed$contexts$tested), interp_join(s_pairs$context, ", "),
      best$context, best$slope_alt, interp_fmt_p(best$p_nominal),
      if (any(is.finite(s_pairs$pip))) sprintf("; highest SuSiE PIP %.2f", max(s_pairs$pip, na.rm = TRUE)) else ""))
    sources <- c(sources, "TOPMed cis-sQTL")
  }
  replicated <- if (!is.null(fhs) && nrow(fhs$sqtl)) interp_fhs_replicated(fhs$sqtl) else NULL
  replicated <- if (!is.null(replicated)) replicated[replicated$cis, , drop = FALSE] else NULL
  fhs_sqtl <- !is.null(replicated) && nrow(replicated) > 0L
  if (fhs_sqtl) {
    evidence <- c(evidence, sprintf("FHS whole-blood isoform ratio: %d transcript(s) of %s replicate; strongest %s, log10 p %.1f (discovery) and %.1f (replication)",
      nrow(replicated), interp_join(replicated$gene, ", "), replicated$feature[[1L]], replicated$log10_p_discovery[[1L]], replicated$log10_p_replication[[1L]]))
    sources <- c(sources, "FHS sQTL")
  }
  g_sqtl <- if (!is.null(gtex)) gtex$qtl[gtex$qtl$qtl_type == "sQTL", , drop = FALSE] else NULL
  gtex_sqtl <- !is.null(g_sqtl) && nrow(g_sqtl) > 0L
  if (gtex_sqtl) {
    evidence <- c(evidence, sprintf("GTEx sQTL in %d tissue(s) for %s; smallest p %s (%s)", length(unique(g_sqtl$tissue)),
      interp_join(g_sqtl$gene, ", ", max = 4), interp_fmt_p(min(g_sqtl$p_value)), g_sqtl$tissue[[which.min(g_sqtl$p_value)]]))
    sources <- c(sources, "GTEx")
  }
  ot_s <- if (!is.null(ot)) ot$qtl[ot$qtl$qtl_type %in% c("sqtl", "tuqtl") & ot$qtl$posterior_probability >= 0.5 &
    !is.na(ot$qtl$posterior_probability), , drop = FALSE] else NULL
  ot_sqtl <- !is.null(ot_s) && nrow(ot_s) > 0L
  if (ot_sqtl) {
    evidence <- c(evidence, sprintf("Open Targets: %d splicing or transcript-usage QTL credible set(s) with PIP >= 0.5 (%s in %s, PIP %.2f)",
      nrow(ot_s), ot_s$gene[[1L]], ot_s$biosample[[1L]], ot_s$posterior_probability[[1L]]))
    sources <- c(sources, "Open Targets")
  }
  associated <- observed || topmed_sqtl || fhs_sqtl || gtex_sqtl
  predicted <- essential || (is.finite(sai_max) && sai_max >= 0.5)
  pinned <- predicted || finemapped || ot_sqtl
  weak <- region || (is.finite(sai_max) && sai_max >= 0.2)
  assessed <- nrow(genic) > 0L || !is.null(topmed) || !is.null(sai) || interp_ok(rna)
  level <- if (associated && pinned) "Strong" else if (associated || predicted || ot_sqtl) "Moderate" else if (weak) "Limited" else
    if (assessed) "None found" else "Not assessed"
  finding <- switch(level,
    Strong = "Splicing differs by genotype, and prediction or fine-mapping points to this variant itself.",
    Moderate = if (associated) "Splicing differs by genotype; this variant is not singled out from its neighbours in LD." else
      "Predicted to alter splicing; no RNA association was found in the sources consulted.",
    Limited = "Near a splice site or weakly predicted; no RNA association found.",
    `None found` = "No splicing effect predicted or observed.", "No splicing source was available.")
  if (is.null(sai)) evidence <- c(evidence, "SpliceAI not run")
  rows[[3L]] <- interp_summary_row("Splicing", level, finding, evidence, sources)

  # 4. Expression.
  evidence <- character(); sources <- character()
  e_pairs <- if (!is.null(topmed)) topmed$pairs[topmed$pairs$qtl_type == "eQTL" & topmed$pairs$significant, , drop = FALSE] else NULL
  topmed_eqtl <- !is.null(e_pairs) && nrow(e_pairs) > 0L
  name <- function(x) ifelse(nzchar(x$gene), x$gene, x$gene_id)
  if (topmed_eqtl) {
    best <- e_pairs[which.min(e_pairs$p_nominal), ]
    evidence <- c(evidence, sprintf("TOPMed eQTL for %s in %s; strongest %s in %s, slope %+.2f, p %s", interp_join(name(e_pairs), ", ", max = 5),
      interp_join(e_pairs$context, ", "), name(best), best$context, best$slope_alt, interp_fmt_p(best$p_nominal)))
    sources <- c(sources, "TOPMed cis-eQTL")
  }
  e_fine <- topmed_eqtl && any(e_pairs$pip >= 0.5, na.rm = TRUE)
  fhs_eqtl <- !is.null(fhs) && nrow(fhs$eqtl) > 0L
  if (fhs_eqtl) {
    evidence <- c(evidence, sprintf("FHS whole-blood eQTL for %s (log10 p %.1f)", interp_join(fhs$eqtl$gene, ", ", max = 5), min(fhs$eqtl$log10_p)))
    sources <- c(sources, "FHS eQTL")
  }
  g_eqtl <- if (!is.null(gtex)) gtex$qtl[gtex$qtl$qtl_type == "eQTL", , drop = FALSE] else NULL
  gtex_eqtl <- !is.null(g_eqtl) && nrow(g_eqtl) > 0L
  if (gtex_eqtl) {
    evidence <- c(evidence, sprintf("GTEx eQTL in %d tissue(s) for %s; smallest p %s (%s)", length(unique(g_eqtl$tissue)),
      interp_join(g_eqtl$gene, ", ", max = 4), interp_fmt_p(min(g_eqtl$p_value)), g_eqtl$tissue[[which.min(g_eqtl$p_value)]]))
    sources <- c(sources, "GTEx")
  }
  ot_e <- if (!is.null(ot)) ot$qtl[ot$qtl$qtl_type %in% c("eqtl", "sceqtl"), , drop = FALSE] else NULL
  ot_eqtl <- !is.null(ot_e) && nrow(ot_e) > 0L
  ot_fine <- ot_eqtl && any(ot_e$posterior_probability >= 0.5, na.rm = TRUE)
  if (ot_eqtl) {
    evidence <- c(evidence, sprintf("Open Targets: %d eQTL credible set(s), highest PIP %.2f (%s, %s)", nrow(ot_e),
      max(ot_e$posterior_probability, na.rm = TRUE), ot_e$gene[[1L]], ot_e$biosample[[1L]]))
    sources <- c(sources, "Open Targets")
  }
  any_e <- topmed_eqtl || fhs_eqtl || gtex_eqtl || ot_eqtl
  level <- if (e_fine || ot_fine) "Strong" else if (any_e) "Moderate" else if (!is.null(topmed) || !is.null(gtex) || ot_queried) "None found" else "Not assessed"
  finding <- switch(level, Strong = "Associated with expression, and fine-mapping gives this variant a posterior probability of at least 0.5.",
    Moderate = "Associated with expression; not fine-mapped to this variant.", `None found` = "No eQTL reported for this variant.",
    "No expression-QTL source was available.")
  rows[[4L]] <- interp_summary_row("Gene expression (eQTL)", level, finding, evidence, sources)

  # 5. Disease and traits.
  evidence <- character(); sources <- character(); level <- "None found"
  rank_up <- function(current, new) interp_levels()[[min(match(current, interp_levels()), match(new, interp_levels()))]]
  clinvar_local <- if (interp_ok(local) && interp_ok(local$clinvar) && isTRUE(local$clinvar$available)) local$clinvar else NULL
  catalog <- if (interp_ok(local) && interp_ok(local$gwas_catalog) && isTRUE(local$gwas_catalog$available)) local$gwas_catalog else NULL
  significance <- character(); detail <- ""
  if (!is.null(myv) && nrow(myv$clinvar)) { significance <- tolower(myv$clinvar$clinical_significance); sources <- c(sources, "ClinVar via MyVariant.info") }
  if (!length(significance) && !is.null(clinvar_local) && isTRUE(clinvar_local$found)) {
    r <- clinvar_local$records
    significance <- tolower(trimws(unlist(strsplit(r$clinical_significance, "[|/,]"))))
    conditions <- setdiff(unique(trimws(unlist(strsplit(r$conditions, ";", fixed = TRUE)))), c("", "not provided", "not specified"))
    detail <- paste0(if (length(conditions)) paste0(" for ", interp_join(conditions, ", ", max = 4)) else "",
      if (any(nzchar(r$review_status))) paste0(" (", interp_join(unique(r$review_status), "; "), ")") else "")
    sources <- c(sources, clinvar_local$source)
  }
  if (!length(significance) && !is.null(known) && nzchar(known$clinvar_significance)) {
    significance <- tolower(gsub("_", " ", strsplit(known$clinvar_significance, ", ", fixed = TRUE)[[1L]]))
    sources <- c(sources, "ClinVar 2023-10 (VEP cache)")
  }
  significance <- unique(significance[nzchar(significance)])
  if (length(significance)) {
    evidence <- c(evidence, paste0("ClinVar: ", interp_join(significance, ", ", max = 6), detail))
    pathogenic <- any(grepl("^(likely )?pathogenic", significance))
    benign <- any(grepl("benign", significance))
    level <- rank_up(level, if (pathogenic && !benign) "Strong" else if (pathogenic || any(grepl("risk|drug|association|protective|affects|conflicting", significance))) "Moderate" else "Limited")
  }
  if (!is.null(ot) && nrow(ot$gwas)) {
    g <- ot$gwas; pip <- ifelse(is.na(g$posterior_probability), 0, g$posterior_probability)
    top <- g[order(-pip, g$log10_p), , drop = FALSE]; top_pip <- sort(pip, decreasing = TRUE)
    traits <- ifelse(nzchar(top$diseases), top$diseases, top$trait)
    keep <- !duplicated(traits)
    evidence <- c(evidence, sprintf("GWAS credible sets: %d (%d with PIP >= 0.5, %d with PIP >= 0.1); top: %s", ot$summary$gwas_total,
      sum(pip >= 0.5), sum(pip >= 0.1), paste(sprintf("%s (PIP %.2f)", head(traits[keep], 5L), head(top_pip[keep], 5L)), collapse = ", ")))
    sources <- c(sources, "Open Targets")
    level <- rank_up(level, if (any(pip >= 0.5)) "Strong" else if (any(pip >= 0.1)) "Moderate" else "Limited")
  }
  if (!is.null(catalog) && isTRUE(catalog$found)) {
    a <- catalog$associations
    genome_wide <- is.finite(a$p_value) & a$p_value <= 5e-8
    traits <- unique(ifelse(nzchar(a$mapped_trait), a$mapped_trait, a$trait))
    evidence <- c(evidence, sprintf("GWAS Catalog: %d association(s) reported at this position, %d at p <= 5e-8, %d distinct trait(s); strongest: %s",
      nrow(a), sum(genome_wide), length(traits), interp_join(traits, "; ", max = 5)))
    sources <- c(sources, trimws(paste("GWAS Catalog", catalog$release)))
    # A published association at the variant itself; alone it does not single the variant out from its LD neighbours.
    level <- rank_up(level, if (any(genome_wide)) "Moderate" else "Limited")
  }
  if (!is.null(ot) && nrow(ot$pharmacogenomics)) {
    evidence <- c(evidence, paste0("Pharmacogenomics: ", interp_join(ot$pharmacogenomics$drugs, ", ", max = 6)))
    level <- rank_up(level, "Moderate"); sources <- c(sources, "Open Targets")
  }
  consulted <- c(if (!is.null(clinvar_local)) clinvar_local$source else if (!is.null(known)) "ClinVar 2023-10 (VEP cache)",
    if (!is.null(catalog)) trimws(paste("GWAS Catalog", catalog$release)), if (!is.null(myv)) "ClinVar via MyVariant.info",
    if (ot_queried) "Open Targets")
  # ClinVar and published associations are both needed before saying "none".
  local_complete <- !is.null(catalog) && (!is.null(clinvar_local) || !is.null(myv))
  if (level == "None found" && !ot_queried && !local_complete) level <- "Not assessed"
  if (!ot_queried) evidence <- c(evidence, "GWAS fine-mapping (credible sets) not consulted: public query not run")
  finding <- switch(level, Strong = "Reported as disease-causing, or fine-mapped as the likely causal variant of a GWAS signal.",
    Moderate = "Linked to a disease, trait or drug response, without being singled out as causal.",
    Limited = "Present in ClinVar or in GWAS records with weak support for this variant.",
    `None found` = paste0("No clinical record and no published GWAS association for this variant in ", interp_join(consulted, ", "), "."),
    if (length(consulted)) paste0("Only ", interp_join(consulted, ", "), " could be consulted and reports nothing; run the public query for GWAS and curated evidence.") else
      "Run the public query to consult ClinVar, GWAS credible sets and curated evidence.")
  rows[[5L]] <- interp_summary_row("Disease and traits", level, finding, evidence, if (length(evidence) > !ot_queried) sources else consulted)

  # 6. Pathway (context for the genes the variant acts through; not evidence about the variant).
  genes <- if (interp_ok(public) && interp_ok(public$genes) && nrow(public$genes$genes)) public$genes else NULL
  lp <- if (interp_ok(local) && interp_ok(local$pathways) && isTRUE(local$pathways$available) && nrow(local$pathways$genes)) local$pathways else NULL
  if (is.null(genes) && is.null(lp)) {
    rows[[6L]] <- interp_summary_row("Pathway", "Not assessed", if (nrow(genic)) paste0("No pathway source was available for ",
      interp_join(genic$gene_name, ", ", max = 4), "; run the public query.") else "No gene overlaps the variant or is a QTL target of it.")
  } else {
    evidence <- character(); sources <- character(); finding <- character()
    if (!is.null(lp)) {
      g <- utils::head(lp$genes, 6L); label <- ifelse(nzchar(g$gene), g$gene, g$gene_id)
      pretty <- function(set) tolower(gsub("_", " ", sub("^[A-Z]+_", "", set)))
      bare <- character()
      for (i in seq_len(nrow(g))) {
        r <- lp$reactome[lp$reactome$gene == label[[i]], , drop = FALSE]
        m <- lp$gene_sets[lp$gene_sets$gene == label[[i]], , drop = FALSE]
        top <- r$pathway[r$top_level]
        # Without a Reactome entry, name the gene's most specific (smallest) gene sets.
        named <- if (length(top)) interp_join(top, ", ", max = 4) else if (nrow(m))
          interp_join(pretty(m$gene_set[order(m$genes_in_set)]), ", ", max = 3) else "no pathway annotation"
        if (!nrow(r) && !nrow(m)) { bare <- c(bare, label[[i]]); next }
        finding <- c(finding, paste0(label[[i]], " (", g$relation[[i]], "): ", named))
        evidence <- c(evidence, sprintf("%s: %d Reactome pathway(s), %d gene set(s)", label[[i]], nrow(r), nrow(m)))
      }
      if (length(bare)) evidence <- c(evidence, paste0("no pathway annotation: ", interp_join(bare, ", ")))
      if (!length(finding)) finding <- "None of the linked genes has a pathway annotation in the local sources."
      if (nrow(lp$genes) > 6L) evidence <- c(evidence, paste0(nrow(lp$genes) - 6L, " more linked gene(s) in the table"))
      sources <- c(sources, lp$source)
    }
    if (!is.null(genes)) {
      if (is.null(lp)) finding <- vapply(genes$genes$gene, function(gene) {
        top <- unique(genes$pathways$top_level[genes$pathways$gene == gene])
        paste0(gene, ": ", if (length(top)) interp_join(top, ", ", max = 5) else "no Reactome pathway")
      }, character(1))
      evidence <- c(evidence, sprintf("Open Targets, %s: %d Reactome pathway(s), %d GO term(s), %d associated disease(s)", genes$genes$gene,
        vapply(genes$genes$gene, function(gene) sum(genes$pathways$gene == gene), integer(1)),
        vapply(genes$genes$gene, function(gene) sum(genes$go$gene == gene), integer(1)), genes$genes$associated_diseases))
      sources <- c(sources, "Reactome and GO via Open Targets")
    }
    rows[[6L]] <- interp_summary_row("Pathway", "Context", paste(finding, collapse = " · "), evidence, sources)
  }

  # 7. Tissue.
  evidence <- character(); sources <- character()
  if (!is.null(topmed) && nrow(topmed$pairs)) {
    h <- interp_context_heterogeneity(topmed$pairs)
    if (nrow(h)) {
      h <- head(h, 4L)
      what <- ifelse(h$qtl_type == "sQTL", sub("^[^:]+:([0-9]+):([0-9]+):.*$", " intron \\1-\\2", h$phenotype), "")
      # The trained cross-context model's reading of the same gene or intron, when a model is available.
      tissue <- if (interp_ok(local) && interp_ok(local$tissue) && isTRUE(local$tissue$available)) local$tissue$phenotypes else NULL
      modelled <- if (!is.null(tissue) && nrow(tissue)) tissue$summary[match(paste(h$qtl_type, h$phenotype), paste(tissue$qtl_type, tissue$phenotype))] else rep(NA_character_, nrow(h))
      evidence <- c(evidence, sprintf("TOPMed %s %s%s: %s%s", h$qtl_type, ifelse(nzchar(h$gene), h$gene, h$gene_id), what, h$pattern,
        ifelse(is.na(modelled), "", paste0("; model: ", modelled))))
      sources <- c(sources, "TOPMed cis-QTL (6 contexts)", if (any(!is.na(modelled))) "cross-context model")
    }
  }
  if (!is.null(gtex) && nrow(gtex$qtl)) {
    evidence <- c(evidence, sprintf("GTEx: QTL in %d of 50 tissues; whole blood %s", length(unique(gtex$qtl$tissue)),
      if ("Whole Blood" %in% gtex$qtl$tissue) "included" else "not among them"))
    sources <- c(sources, "GTEx")
  }
  spec <- if (!is.null(gtex)) interp_tissue_specificity(gtex$expression) else NULL
  if (!is.null(spec) && nrow(spec)) {
    evidence <- c(evidence, sprintf("%s expression: highest in %s (%.0f TPM); whole blood %.1f TPM, rank %d of %d; tau %.2f", spec$gene, spec$top_tissue,
      spec$top_median_tpm, spec$whole_blood_tpm, spec$whole_blood_rank, spec$tissues, spec$tau))
    sources <- c(sources, "GTEx")
  }
  rows[[7L]] <- if (length(evidence)) interp_summary_row("Tissue", "Context",
    "Where the effect and the gene are seen; this cohort's RNA is whole blood.", evidence, sources) else
    interp_summary_row("Tissue", "Not assessed", "No tissue source was available or the variant is not a QTL in any of them.")
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# Fill gene names into TOPMed pairs from a gene_id -> gene_name table.
interp_name_genes <- function(pairs, names) {
  if (is.null(pairs) || !nrow(pairs) || is.null(names) || !nrow(names)) return(pairs)
  pairs$gene <- names$gene_name[match(interp_bare_gene(pairs$gene_id), interp_bare_gene(names$gene_id))]
  pairs$gene[is.na(pairs$gene)] <- ""
  pairs
}
# ---- Cross-context sharing model ---------------------------------------------
# An empirical-Bayes model of how a cis-QTL effect is shared between contexts,
# in the manner of multivariate adaptive shrinkage (Urbut et al. 2019). For one
# variant-phenotype pair with estimates x_hat and standard errors s in the
# contexts where it was tested,
#   x_hat ~ N(x, S V S),   x ~ pi_0 * delta_0 + sum_kl pi_kl * N(0, w_l^2 * U_k)
# S = diag(s). V is the correlation of the estimates between contexts when there
# is no effect (shared donors); U_k are patterns of sharing (the same in all
# contexts, one context only, patterns learned from the strongest signals);
# w_l is a grid of effect sizes. A model works on one of two scales:
#   "slope": x is the per-allele slope, so contexts share effect sizes;
#   "z":     x is slope / s and S is the identity, so contexts share z-scores.
# scripts/train_tissue_model.R learns V, the learned U_k and pi from the TOPMed
# cis-QTL files and keeps the scale that predicts held-out pairs better; the
# functions here are used both there and at view time. A context that was not
# tested is left out of the likelihood.

interp_mash_grid <- function(scale = "z") if (identical(scale, "slope")) 2^(seq(-10, 3) / 2) else 2^(seq(-2, 12) / 2)

# Fixed sharing patterns for R contexts.
interp_mash_canonical <- function(contexts) {
  r <- length(contexts)
  out <- list(identity = diag(r), equal_effects = matrix(1, r, r))
  for (rho in c(.25, .5, .75)) { u <- matrix(rho, r, r); diag(u) <- 1; out[[paste0("simple_het_", rho)]] <- u }
  for (i in seq_len(r)) { u <- matrix(0, r, r); u[i, i] <- 1; out[[paste0("only_", contexts[[i]])]] <- u }
  lapply(out, function(u) { dimnames(u) <- list(contexts, contexts); u })
}

# Class of a pattern, for reporting which kind carries the posterior weight.
interp_mash_class <- function(name) ifelse(name == "null", "no effect", ifelse(startsWith(name, "only_"), "one context",
  ifelse(name %in% c("equal_effects") | startsWith(name, "simple_het_"), "shared", ifelse(name == "identity", "independent", "learned"))))

interp_mash_model <- function(contexts, V, U, pi = NULL, grid = interp_mash_grid(scale), meta = list(), scale = "z") {
  r <- length(contexts)
  stopifnot(scale %in% c("z", "slope"), is.matrix(V), all(dim(V) == r), length(U) > 0L,
    all(vapply(U, function(u) is.matrix(u) && all(dim(u) == r), logical(1))))
  n <- length(U) * length(grid) + 1L
  if (is.null(pi)) pi <- rep(1 / n, n)
  stopifnot(length(pi) == n, all(pi >= 0), abs(sum(pi) - 1) < 1e-6)
  structure(list(contexts = contexts, scale = scale, V = V, U = U, grid = grid, pi = pi, meta = meta), class = "regshiny_tissue_model")
}

# Component index: 1 is the null, then every pattern k at every scale l.
interp_mash_components <- function(model) data.frame(index = seq_along(model$pi),
  pattern = c("null", rep(names(model$U), each = length(model$grid))), scale = c(0, rep(model$grid, times = length(model$U))),
  stringsAsFactors = FALSE)

interp_mash_logdmvnorm <- function(X, S) {
  # Rows of X ~ N(0, S); returns one log density per row.
  ch <- chol(S)
  y <- backsolve(ch, t(X), transpose = TRUE)
  -0.5 * colSums(y^2) - sum(log(diag(ch))) - 0.5 * ncol(X) * log(2 * pi)
}

# Sampling covariance of the estimates in the tested contexts `o`.
interp_mash_noise <- function(model, o, se = NULL) {
  v <- model$V[o, o, drop = FALSE]
  if (identical(model$scale, "slope")) v * outer(se[o], se[o]) else v
}

# n x components matrix of log-likelihoods. X is n x R with NA for untested
# contexts; SE (same shape) is needed by a slope-scale model.
interp_mash_loglik <- function(X, model, SE = NULL) {
  X <- as.matrix(X); comps <- interp_mash_components(model); n <- nrow(comps)
  out <- matrix(NA_real_, nrow(X), n)
  prior <- lapply(seq_len(n)[-1L], function(j) comps$scale[[j]]^2 * model$U[[comps$pattern[[j]]]])
  if (identical(model$scale, "slope")) {
    stopifnot(!is.null(SE), all(dim(SE) == dim(X)))
    for (i in seq_len(nrow(X))) {
      o <- which(is.finite(X[i, ]) & is.finite(SE[i, ]))
      if (!length(o)) { out[i, ] <- 0; next }
      x <- matrix(X[i, o], 1L); noise <- interp_mash_noise(model, o, SE[i, ])
      out[i, 1L] <- interp_mash_logdmvnorm(x, noise)
      for (j in seq_len(n)[-1L]) out[i, j] <- interp_mash_logdmvnorm(x, noise + prior[[j - 1L]][o, o, drop = FALSE])
    }
    return(out)
  }
  # z scale: the covariance depends only on which contexts were tested.
  pattern <- apply(is.finite(X), 1L, function(x) paste(which(x), collapse = ","))
  for (p in unique(pattern)) {
    rows <- which(pattern == p); o <- as.integer(strsplit(p, ",", fixed = TRUE)[[1L]])
    if (!length(o)) { out[rows, ] <- 0; next }
    x <- X[rows, o, drop = FALSE]; noise <- interp_mash_noise(model, o)
    out[rows, 1L] <- interp_mash_logdmvnorm(x, noise)
    for (j in seq_len(n)[-1L]) out[rows, j] <- interp_mash_logdmvnorm(x, noise + prior[[j - 1L]][o, o, drop = FALSE])
  }
  out
}

# Mixture weights by EM on a likelihood matrix. null_weight > 1 is a penalty
# that favours "no effect" when the data do not decide, as in ashr.
interp_mash_fit_pi <- function(loglik, null_weight = 10, iterations = 500L, tolerance = 1e-7) {
  lik <- exp(loglik - apply(loglik, 1L, max))
  k <- ncol(lik); pi <- rep(1 / k, k); prior <- c(null_weight, rep(1, k - 1L)) - 1
  objective <- -Inf
  for (i in seq_len(iterations)) {
    w <- sweep(lik, 2L, pi, `*`); total <- rowSums(w)
    new <- sum(log(total)) + sum(prior * log(pmax(pi, 1e-300)))
    pi <- colSums(w / total) + prior; pi <- pi / sum(pi)
    if (abs(new - objective) < tolerance * abs(new)) break
    objective <- new
  }
  pi
}

# Posterior of the true effect in every context of the model for one pair, on
# the model's scale. One row per context: for a tested context the re-estimated
# effect, for an untested one what the tested contexts imply. lfsr is the local
# false sign rate: the posterior probability that the effect is absent or has
# the other sign.
interp_mash_posterior <- function(x_hat, model, se = NULL) {
  r <- length(model$contexts); comps <- interp_mash_components(model)
  slope <- identical(model$scale, "slope")
  o <- which(is.finite(x_hat) & (if (slope) is.finite(se) & se > 0 else TRUE))
  out <- data.frame(context = model$contexts, tested = seq_len(r) %in% o, estimate = as.numeric(x_hat), posterior_mean = NA_real_,
    posterior_sd = NA_real_, p_positive = NA_real_, p_negative = NA_real_, lfsr = NA_real_, stringsAsFactors = FALSE)
  if (!length(o)) return(list(contexts = out, weights = data.frame(class = character(), weight = numeric(), stringsAsFactors = FALSE),
    top_single = NULL, log_bf = NA_real_))
  x <- x_hat[o]; noise <- interp_mash_noise(model, o, se)
  n <- nrow(comps); ll <- numeric(n); mu <- matrix(0, n, r); va <- matrix(0, n, r)
  ll[[1L]] <- interp_mash_logdmvnorm(matrix(x, 1L), noise)
  for (j in seq_len(n)[-1L]) {
    P <- comps$scale[[j]]^2 * model$U[[comps$pattern[[j]]]]
    S <- noise + P[o, o, drop = FALSE]
    ll[[j]] <- interp_mash_logdmvnorm(matrix(x, 1L), S)
    A <- P[, o, drop = FALSE] %*% solve(S)
    mu[j, ] <- A %*% x
    va[j, ] <- pmax(0, diag(P) - rowSums(A * P[, o, drop = FALSE]))
  }
  lw <- log(pmax(model$pi, 1e-300)) + ll; w <- exp(lw - max(lw)); w <- w / sum(w)
  mean <- colSums(w * mu); second <- colSums(w * (va + mu^2))
  sd <- sqrt(va); positive <- ifelse(sd > 0, stats::pnorm(mu / ifelse(sd > 0, sd, 1)), 0)
  negative <- ifelse(sd > 0, stats::pnorm(-mu / ifelse(sd > 0, sd, 1)), 0)
  out$posterior_mean <- mean; out$posterior_sd <- sqrt(pmax(0, second - mean^2))
  out$p_positive <- colSums(w * positive); out$p_negative <- colSums(w * negative)
  out$lfsr <- pmin(1, pmax(0, 1 - pmax(out$p_positive, out$p_negative)))
  class <- interp_mash_class(comps$pattern)
  weights <- stats::aggregate(weight ~ class, data.frame(class = class, weight = w, stringsAsFactors = FALSE), sum)
  # The single one-context pattern with the most weight, if any.
  one <- tapply(w[class == "one context"], comps$pattern[class == "one context"], sum)
  list(contexts = out, weights = weights[order(-weights$weight), , drop = FALSE],
    top_single = if (length(one)) list(context = sub("^only_", "", names(one)[[which.max(one)]]), weight = max(one)) else NULL,
    # Bayes factor of "some effect" against "no effect in any tested context".
    log_bf = max(lw[-1L]) + log(sum(exp(lw[-1L] - max(lw[-1L])))) - log(max(sum(model$pi[-1L]), 1e-300)) - ll[[1L]])
}

# What the other contexts predict for the observed estimate of context `leave`
# (which must be tested). Used to check the model on pairs it was not fitted on.
interp_mash_predict_left_out <- function(x_hat, model, leave, se = NULL) {
  slope <- identical(model$scale, "slope")
  tested <- which(is.finite(x_hat) & (if (slope) is.finite(se) & se > 0 else TRUE))
  o <- setdiff(tested, leave); comps <- interp_mash_components(model)
  noise <- interp_mash_noise(model, seq_along(model$contexts), if (slope) ifelse(is.finite(se), se, 1) else NULL)
  if (!length(o)) return(c(mean = 0, sd = sqrt(noise[leave, leave])))
  x <- x_hat[o]; n <- nrow(comps); ll <- numeric(n); mu <- numeric(n); va <- numeric(n)
  for (j in seq_len(n)) {
    total <- noise + if (j == 1L) 0 else comps$scale[[j]]^2 * model$U[[comps$pattern[[j]]]]
    S <- total[o, o, drop = FALSE]
    ll[[j]] <- interp_mash_logdmvnorm(matrix(x, 1L), S)
    a <- total[leave, o, drop = FALSE] %*% solve(S)
    mu[[j]] <- a %*% x; va[[j]] <- max(0, total[leave, leave] - sum(a * total[leave, o]))
  }
  lw <- log(pmax(model$pi, 1e-300)) + ll; w <- exp(lw - max(lw)); w <- w / sum(w)
  m <- sum(w * mu)
  c(mean = m, sd = sqrt(max(0, sum(w * (va + mu^2)) - m^2)))
}

# Correlation of z between contexts among pairs with no sign of an effect
# (|z| < 2 in both), from pairs chosen without looking at their p-values.
interp_mash_null_correlation <- function(Z, threshold = 2, min_pairs = 50L) {
  r <- ncol(Z); V <- diag(r)
  for (i in seq_len(r - 1L)) for (j in (i + 1L):r) {
    keep <- is.finite(Z[, i]) & is.finite(Z[, j]) & abs(Z[, i]) < threshold & abs(Z[, j]) < threshold
    if (sum(keep) >= min_pairs) V[i, j] <- V[j, i] <- stats::cor(Z[keep, i], Z[keep, j])
  }
  interp_mash_nearest_pd(V, keep_diagonal = TRUE)
}

interp_mash_nearest_pd <- function(M, floor = 1e-6, keep_diagonal = FALSE) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  out <- e$vectors %*% diag(pmax(e$values, floor), nrow(M)) %*% t(e$vectors)
  if (keep_diagonal) { d <- sqrt(diag(out)); out <- out / outer(d, d) }
  dimnames(out) <- dimnames(M); out
}

# Sharing patterns learned from the strongest signals: the second moment of
# their estimates with the sampling noise removed, and its leading directions.
# SE is given for slope-scale estimates and left out for z-scores.
interp_mash_learned <- function(X, V, contexts, SE = NULL, n_directions = 3L) {
  r <- ncol(X); M <- matrix(0, r, r)
  for (i in seq_len(r)) for (j in i:r) {
    keep <- is.finite(X[, i]) & is.finite(X[, j])
    noise <- if (is.null(SE)) V[i, j] else V[i, j] * mean(SE[keep, i] * SE[keep, j])
    M[i, j] <- M[j, i] <- if (sum(keep) > 1L) mean(X[keep, i] * X[keep, j]) - noise else 0
  }
  e <- eigen(M, symmetric = TRUE); e$values <- pmax(e$values, 0)
  scale <- function(u) { u <- u / max(diag(u)); dimnames(u) <- list(contexts, contexts); u }
  out <- list(learned_all = scale(e$vectors %*% diag(e$values, r) %*% t(e$vectors)))
  k <- min(n_directions, sum(e$values > 0))
  for (i in seq_len(k)) out[[paste0("learned_direction_", i)]] <- scale(tcrossprod(e$vectors[, i]))
  if (k > 1L) out$learned_low_rank <- scale(e$vectors[, seq_len(k), drop = FALSE] %*% diag(e$values[seq_len(k)], k) %*% t(e$vectors[, seq_len(k), drop = FALSE]))
  out
}

interp_tissue_model_write <- function(model, path) {
  body <- list(format = "regshiny_tissue_model_v1", scale = model$scale, contexts = model$contexts, V = model$V, patterns = names(model$U),
    U = unname(model$U), grid = model$grid, pi = model$pi, meta = model$meta)
  tmp <- paste0(path, ".tmp")
  writeLines(jsonlite::toJSON(body, digits = NA, auto_unbox = TRUE, pretty = TRUE, na = "null"), tmp)
  invisible(file.rename(tmp, path))
}

interp_tissue_model_read <- function(path) {
  x <- jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyMatrix = TRUE, simplifyDataFrame = TRUE)
  if (!identical(x$format, "regshiny_tissue_model_v1")) stop("Not a tissue sharing model file: ", path)
  r <- length(x$contexts)
  U <- if (is.array(x$U) && length(dim(x$U)) == 3L) lapply(seq_len(dim(x$U)[[1L]]), function(i) x$U[i, , ]) else lapply(x$U, function(u) matrix(unlist(u), r, r, byrow = TRUE))
  names(U) <- x$patterns
  interp_mash_model(x$contexts, matrix(unlist(x$V), r, r), U, as.numeric(x$pi), as.numeric(x$grid), x$meta,
    scale = if (is.null(x$scale)) "z" else x$scale)
}

# -- The model applied to one variant ---------------------------------------------
interp_tissue_empty <- function() list(
  contexts = data.frame(qtl_type = character(), phenotype = character(), gene = character(), context = character(), tested = logical(),
    slope_alt = numeric(), slope_se = numeric(), z = numeric(), posterior_slope = numeric(), posterior_sd = numeric(), posterior_z = numeric(),
    lfsr = numeric(), call = character(), stringsAsFactors = FALSE),
  phenotypes = data.frame(qtl_type = character(), phenotype = character(), gene = character(), contexts_tested = integer(),
    effect_in = character(), no_clear_effect_in = character(), predicted_in = character(), similar_size_in = character(),
    leading_pattern = character(), pattern_weight = numeric(), log10_bayes_factor = numeric(), summary = character(), stringsAsFactors = FALSE))

interp_tissue_models <- function(sources = interp_local_sources()) {
  files <- sources$tissue_model
  if (!length(files)) return(list())
  files <- files[vapply(files, interp_readable, logical(1))]
  out <- lapply(files, function(f) tryCatch(interp_tissue_model_read(f), error = function(e) e))
  out[!vapply(out, inherits, logical(1), "error")]
}

# Posterior effect of this variant in every context of the model, for each gene
# (eQTL) or intron (sQTL) that is significant in at least one context. `lfsr`
# below `cut` is called an effect. For a context where the pair was not tested
# the model gives what the tested contexts imply (a slope from a slope-scale model, z units from a z-scale one).
interp_tissue_sharing <- function(pairs, sources = interp_local_sources(), models = interp_tissue_models(sources), cut = 0.05) {
  out <- c(list(available = length(models) > 0L, models = lapply(models, function(m) c(m$meta, list(scale = m$scale))), cut = cut), interp_tissue_empty())
  if (!length(models) || is.null(pairs) || !nrow(pairs)) return(out)
  pairs$key <- interp_phenotype_key(pairs$qtl_type, pairs$phenotype)
  ctx_rows <- list(); ph_rows <- list()
  groups <- split(pairs, paste(pairs$qtl_type, pairs$key))
  groups <- groups[order(vapply(groups, function(x) min(x$p_nominal, na.rm = TRUE), numeric(1)))]
  for (x in groups) {
    model <- models[[x$qtl_type[[1L]]]]
    x <- x[is.finite(x$slope_alt) & is.finite(x$slope_se) & x$slope_se > 0 & !duplicated(x$context), , drop = FALSE]
    if (is.null(model) || !nrow(x) || !any(x$significant %in% TRUE)) next
    k <- match(model$contexts, x$context)
    if (all(is.na(k))) next
    z <- x$slope_alt[k] / x$slope_se[k]
    # A slope-scale model re-estimates the slope itself, also where the pair was not tested;
    # a z-scale model re-estimates z, which converts to a slope only where a standard error exists.
    on_slope <- identical(model$scale, "slope")
    post <- if (on_slope) interp_mash_posterior(x$slope_alt[k], model, x$slope_se[k]) else interp_mash_posterior(z, model)
    p <- post$contexts
    effect <- is.finite(p$lfsr) & p$lfsr < cut
    gene <- if (nzchar(x$gene[[1L]])) x$gene[[1L]] else interp_bare_gene(x$gene_id[[1L]])
    slope <- if (on_slope) p$posterior_mean else p$posterior_mean * x$slope_se[k]
    slope_sd <- if (on_slope) p$posterior_sd else p$posterior_sd * x$slope_se[k]
    ctx_rows[[length(ctx_rows) + 1L]] <- data.frame(qtl_type = x$qtl_type[[1L]], phenotype = x$key[[1L]], gene = gene, context = p$context,
      tested = p$tested, slope_alt = x$slope_alt[k], slope_se = x$slope_se[k], z = z, posterior_slope = slope, posterior_sd = slope_sd,
      posterior_z = if (on_slope) slope / x$slope_se[k] else p$posterior_mean, lfsr = p$lfsr,
      call = ifelse(p$tested, ifelse(effect, "effect", "no clear effect"), ifelse(effect, "predicted effect", "not predicted")), stringsAsFactors = FALSE)
    # Similar size: tested contexts with an effect whose posterior slope has the sign of the largest and at least half its size.
    lead <- if (any(effect & p$tested)) max(abs(slope[effect & p$tested])) else NA_real_
    lead_sign <- if (is.finite(lead)) sign(slope[effect & p$tested][[which.max(abs(slope[effect & p$tested]))]]) else 0
    similar <- effect & p$tested & is.finite(slope) & sign(slope) == lead_sign & abs(slope) >= lead / 2
    w <- post$weights
    leading <- if (nrow(w)) w$class[[1L]] else ""
    if (identical(leading, "one context") && !is.null(post$top_single)) leading <- paste0("one context (", post$top_single$context, ")")
    n_tested <- sum(p$tested); n_effect <- sum(effect & p$tested)
    summary <- if (!n_effect) "no tested context has a confident effect" else if (n_effect == n_tested)
      paste0("effect in all ", n_tested, " tested contexts", if (sum(similar) < n_effect) paste0(", similar size in ", sum(similar)) else ", similar size") else
      paste0("effect in ", n_effect, " of ", n_tested, " tested contexts (", paste(p$context[effect & p$tested], collapse = ", "), ")")
    if (any(effect & !p$tested)) summary <- paste0(summary, "; predicted in ", paste(p$context[effect & !p$tested], collapse = ", "))
    ph_rows[[length(ph_rows) + 1L]] <- data.frame(qtl_type = x$qtl_type[[1L]], phenotype = x$key[[1L]], gene = gene, contexts_tested = n_tested,
      effect_in = paste(p$context[effect & p$tested], collapse = ", "), no_clear_effect_in = paste(p$context[!effect & p$tested], collapse = ", "),
      predicted_in = paste(p$context[effect & !p$tested], collapse = ", "), similar_size_in = paste(p$context[similar], collapse = ", "),
      leading_pattern = leading, pattern_weight = if (nrow(w)) w$weight[[1L]] else NA_real_, log10_bayes_factor = post$log_bf / log(10),
      summary = summary, stringsAsFactors = FALSE)
  }
  if (length(ph_rows)) { out$contexts <- do.call(rbind, ctx_rows); out$phenotypes <- do.call(rbind, ph_rows); rownames(out$contexts) <- rownames(out$phenotypes) <- NULL }
  out
}
