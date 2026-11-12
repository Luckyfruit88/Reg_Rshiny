# Rscript tests/test_backend.R [--unit]
# backend.R: CIGAR/coordinate parsing, BAM contig choice and RNA base counting.
# --unit skips the native samtools / RegTools section and the optional real-BAM
# header check (REGTOOLS_TEST_BAM).
unit_only <- "--unit" %in% commandArgs(trailingOnly = TRUE)
source("src/backend.R", encoding = "UTF-8")
cat("Core CIGAR / coordinate tests\n")
a <- parse_cigar("20M100N30M", 181L)
stopifnot(identical(as.integer(a$introns[1, ]), c(200L, 300L)))
stopifnot(identical(as.integer(a$blocks[1, ]), c(180L, 200L)))
stopifnot(nrow(parse_cigar("20M100D30M", 181L)$introns) == 0L)
stopifnot(nrow(parse_cigar("20M100N30M100N20M", 181L)$introns) == 2L)
b <- parse_cigar("5S10=2X8M3I100N20M2D5M", 181L)
stopifnot(identical(as.integer(b$introns[1, ]), c(200L, 300L)))
# Same-fragment paired ends both infer the same transcript strand.
stopifnot(all(infer_regtools_strand(c(99L, 147L), c(NA, NA), "RF") == "-"))
stopifnot(all(infer_regtools_strand(c(99L, 147L), c(NA, NA), "FR") == "+"))
stopifnot(infer_regtools_strand(0L, NA_character_, "XS") == "?")

s <- tempfile(fileext = ".sam"); make_demo_sam(s)
lines <- readLines(s); records <- lines[!startsWith(lines, "@")]
cfg <- default_demo_config()
# Filter flags and MAPQ as samtools would before parse_sam_records().
fields <- strsplit(records, "\t", fixed = TRUE)
keep <- vapply(fields, function(x) as.integer(x[[5]]) >= 20L &&
                 bitwAnd(as.integer(x[[2]]), 2820L) == 0L, logical(1))
p <- parse_sam_records(records[keep], cfg)
stopifnot(nrow(p$reads) == 4L, p$span_only_excluded == 1L, nrow(p$events) == 3L)
unlink(s)
# Junction scope: the exonic base next to either boundary inside the interval.
# Interval 201-450; boundaries are 0-based intron start (= last exonic base,
# 1-based) and 0-based exclusive end (flanking exonic base = end + 1).
stopifnot(identical(junction_in_scope(c(270L, 430L, 120L, 20L, 180L, 200L, 450L, 451L, 100L),
                                      c(370L, 730L, 220L, 100L, 480L, 449L, 900L, 900L, 199L), 201L, 450L),
                    c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)),
          is.character(junction_scope_rule()), nzchar(junction_scope_rule()))
cat("Core tests: PASS\n")


x <- parse_bam_contigs(c("@HD\tVN:1.6", "@SQ\tSN:chr2\tLN:242193529", "@SQ\tSN:chr1\tLN:248956422", "@SQ\tSN:chrX\tLN:156040895"))
stopifnot(identical(x$chrom, c("chr2","chr1","chrX")), choose_bam_contig(x) == "chr1",
          choose_bam_contig(x, "chrX") == "chrX", choose_bam_contig(x, "chrDemo") == "chr1")
y <- parse_bam_contigs(c("@SQ\tSN:2\tLN:200", "@SQ\tSN:1\tLN:100"))
stopifnot(choose_bam_contig(y) == "1")
z <- parse_bam_contigs("@SQ\tSN:scaffoldA\tLN:123")
stopifnot(choose_bam_contig(z) == "scaffoldA")
bad <- list(character(), "@SQ\tSN:chr1", "@SQ\tSN:chr1\tLN:no",
            c("@SQ\tSN:chr1\tLN:100", "@SQ\tSN:chr1\tLN:100"))
stopifnot(all(vapply(bad, function(h) inherits(try(parse_bam_contigs(h), silent=TRUE), "try-error"), logical(1))))
bam <- Sys.getenv("REGTOOLS_TEST_BAM")
if (!unit_only && nzchar(bam)) {
  actual <- read_bam_contigs(bam)
  stopifnot(all(c("chr1","chr22","chrX","chrY") %in% actual$chrom),
            !"chrDemo" %in% actual$chrom, choose_bam_contig(actual, "chrDemo") == "chr1")
  cat("Real BAM header choices: PASS; contigs=", nrow(actual), "\n", sep="")
}
cat("Contig selection tests: PASS\n")

rec <- function(cigar, seq, qual = paste(rep("I", nchar(seq)), collapse = ""), flag = 0L) {
  paste(c("synthetic", flag, "chrDemo", 100, 60, cigar, "*", 0, 0, seq, qual), collapse = "\t")
}
# Soft clips and insertions consume query, N and D consume reference only.
x <- count_rna_bases(c(rec("2S3M2I2M", "TTAAAGGCT"),
                       rec("3M2D2M", "AAACC"),
                       rec("3M2N2M", "AAACC"),
                       rec("5M", "AAAGT", flag = 16L),
                       rec("5M", "AAAAT", "III!I"),
                       rec("5M", "AAAAA", "*")), 103L, 20L)
stopifnot(x$count[x$base == "C"] == 1L, x$count[x$base == "G"] == 1L,
          sum(x$count) == 2L, attr(x, "excluded")[["deleted_D"]] == 1L,
          attr(x, "excluded")[["skipped_N"]] == 1L,
          attr(x, "excluded")[["low_baseq"]] == 1L,
          attr(x, "excluded")[["missing_quality"]] == 1L)
stopifnot(sum(count_rna_bases(character(), 103L)$count) == 0L)
stopifnot(sum(count_rna_bases(rec("5M", "AAAAA"), 110L)$count) == 0L)
cat("RNA base CIGAR, reverse-strand, BQ and missingness tests: PASS\n")

if (unit_only) {
  cat("Native integration tests: skipped (--unit)\n")
} else {
  cat("Native samtools / RegTools integration tests\n")
  r <- analyze_bam(cfg)
  m <- summarize_target(r, 201L, 300L, 5L)
  stopifnot(r$native_audit_passed, m$denominator_reads == 4L,
            m$junction_reads == 3L, m$junction_exact_count == 2L,
            m$junction_near_count == 1L, m$junction_per_100_reads == 50,
            abs(m$mean_read_depth - 200/350) < 1e-8,
            sum(r$depth$depth) == 200L, r$span_only_excluded == 1L)
  stopifnot(summarize_target(r, 201L, 300L, 0L)$junction_near_count == 0L)
  stopifnot(summarize_target(r, 201L, 300L, 1L)$junction_near_count == 0L)
  stopifnot(summarize_target(r, 201L, 300L, 2L)$junction_near_count == 1L)
  stopifnot(inherits(try(summarize_target(r, 101L, 300L, 5L), silent = TRUE), "try-error"))
  # Tool versions and help are asked once per process and executable, and a
  # repeated analysis returns the same versions, help text and command log.
  facts <- length(ls(.tool_facts)); stopifnot(facts == 3L)
  again <- analyze_bam(cfg)
  stopifnot(length(ls(.tool_facts)) == facts, identical(again$versions, r$versions),
            identical(gsub("regtools_shiny_[A-Za-z0-9]+", "W", again$log), gsub("regtools_shiny_[A-Za-z0-9]+", "W", r$log)),
            any(again$log == "samtools '--version'"), isTRUE(warm_tool_facts()), length(ls(.tool_facts)) == facts)
  # No retained reads: no-junction is zero, but a zero-denominator rate is NA.
  cfg$mapq <- 61L
  z <- analyze_bam(cfg); mz <- summarize_target(z, 201L, 300L, 5L)
  stopifnot(mz$denominator_reads == 0L, mz$junction_exact_count == 0L,
            is.na(mz$junction_per_100_reads), sum(z$depth$depth) == 0L)
  # A junction that leaves the interval is reported in full when the exonic
  # base next to one boundary is inside; one with both flanks outside is not,
  # even when its read is retained for another block.
  scope_dir <- tempfile("scope_"); dir.create(scope_dir)
  rec <- function(name, pos, cigar, len) paste(c(name, 0, "chrDemo", pos, 60, cigar, "*", 0, 0,
    paste(rep("A", len), collapse = ""), paste(rep("I", len), collapse = ""), "NH:i:1", "XS:A:+"), collapse = "\t")
  writeLines(c("@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chrDemo\tLN:1000",
    rec("out_then_in", 1, "20M80N20M100N30M", 70L),   # introns 21-100 (out of scope), 121-220 (right flank inside)
    rec("leaves", 401, "30M300N20M", 50L),            # intron 431-730, left flank inside
    rec("leaves2", 411, "20M300N20M", 40L),           # same junction, second read
    rec("span_only", 151, "30M300N30M", 60L),         # intron 181-480, no block inside
    rec("inside", 251, "20M100N30M", 50L)),           # intron 271-370
    file.path(scope_dir, "scope.sam"))
  scope_bam <- file.path(scope_dir, "scope.bam")
  for (args in list(c("sort", "-o", scope_bam, file.path(scope_dir, "scope.sam")), c("index", scope_bam)))
    processx::run("samtools", args)
  scope_cfg <- default_demo_config(); scope_cfg$demo <- FALSE; scope_cfg$bam <- scope_bam
  scope_cfg$start1 <- 201L; scope_cfg$end1 <- 450L
  sr <- analyze_bam(scope_cfg); sj <- sr$junctions[order(sr$junctions$intron_start1), ]
  stopifnot(sr$native_audit_passed, nrow(sr$reads) == 4L, sr$span_only_excluded == 1L,
            identical(sj$intron_start1, c(121L, 271L, 431L)), identical(sj$intron_end1, c(220L, 370L, 730L)),
            identical(sj$score, c(1L, 1L, 2L)), nrow(sr$events) == 4L)
  unlink(scope_dir, recursive = TRUE)
  fake_dir <- tempfile("fact_"); dir.create(fake_dir); fake <- file.path(fake_dir, "tool"); calls <- file.path(fake_dir, "calls")
  write_fake <- function(answer, status) { writeLines(c("#!/bin/sh", paste0("echo x >> ", shQuote(calls)),
    paste("echo", answer), paste("exit", status)), fake); Sys.chmod(fake, "0700") }
  n_calls <- function() if (file.exists(calls)) length(readLines(calls)) else 0L
  write_fake("version-one", 0L)
  stopifnot(trimws(tool_fact(fake, "--version")$stdout) == "version-one", trimws(tool_fact(fake, "--version")$stdout) == "version-one",
            n_calls() == 1L, trimws(tool_fact(fake, "-h")$stdout) == "version-one", n_calls() == 2L)
  # A replaced executable is asked again; a failing answer is never kept.
  write_fake("version-two-longer", 0L)
  stopifnot(trimws(tool_fact(fake, "--version")$stdout) == "version-two-longer", n_calls() == 3L)
  write_fake("broken-tool-answer-x", 1L)
  stopifnot(tool_fact(fake, "--version")$status == 1L, tool_fact(fake, "--version")$status == 1L, n_calls() == 5L)
  unlink(fake_dir, recursive = TRUE)
  # Per-sample tool calls must not scan the node's process table (cleanup_tree).
  stopifnot(!any(grepl("cleanup_tree = TRUE", deparse(analyze_bam), fixed = TRUE)))
  cat("Native integration tests: PASS\n")
}
