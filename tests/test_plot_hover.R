# Hover read-outs: marks are recorded only inside hover_capture(), as fractions
# of the device, per panel, with labels that carry values but no sample IDs.
options(warn = 2)
suppressPackageStartupMessages(library(shiny))
source("src/plot_hover.R")
source("src/annotation_module.R")
source("src/splice_evidence_module.R")
source("src/interpretation_module.R")
source("src/variant_module.R")

grDevices::png(tempfile(fileext = ".png"), width = 800, height = 600)
on_device <- function(code) code

# Outside a capture nothing is recorded and plotting is unchanged.
graphics::plot(1:3, 1:3)
stopifnot(!hover_recording(), is.null(hover_points(1, 1, "a")), is.null(.hover_state$groups))

# Points: fractions of the device, inside 0-1, increasing with the data.
marks <- hover_capture({ graphics::plot(1:3, 1:3); hover_points(1:3, 1:3, paste("p", 1:3)); hover_points(NA, 1, "dropped") })
stopifnot(length(marks) == 1L, marks[[1L]]$kind == "xy", identical(unclass(marks[[1L]]$label), paste("p", 1:3)),
          all(marks[[1L]]$x > 0 & marks[[1L]]$x < 1), !is.unsorted(marks[[1L]]$x), !is.unsorted(marks[[1L]]$y),
          !hover_recording())
# Round trip: the stored fraction is the device position of the data value.
stopifnot(isTRUE(all.equal(graphics::grconvertX(unclass(marks[[1L]]$x), "ndc", "user"), c(1, 2, 3))))

# Panels of one figure get their own positions: the second panel lies to the right.
marks <- hover_capture({
  old <- graphics::par(mfrow = c(1L, 2L))
  graphics::plot(0:1, 0:1); hover_points(.5, .5, "left"); hover_columns(c(0, 1), c("a", "b"))
  graphics::plot(0:1, 0:1); hover_points(.5, .5, "right")
  graphics::par(old)
})
stopifnot(length(marks) == 3L, marks[[1L]]$x < .5, marks[[3L]]$x > .5, marks[[2L]]$kind == "x",
          length(marks[[2L]]$box) == 4L, marks[[2L]]$box[[2L]] < .5, marks[[2L]]$box[[1L]] < marks[[2L]]$box[[2L]])

# Long tracks are thinned; every kept label still names its own position.
marks <- hover_capture({ graphics::plot(c(1, 10000), 0:1); hover_columns(1:10000, as.character(1:10000), max_n = 100L) })
kept <- as.integer(unclass(marks[[1L]]$label))
stopifnot(length(kept) <= 100L, kept[[1L]] == 1L, kept[[length(kept)]] == 10000L,
          isTRUE(all.equal(graphics::grconvertX(unclass(marks[[1L]]$x), "ndc", "user"), kept)))

# Rectangles are normalised to lower-left / upper-right.
marks <- hover_capture({ graphics::plot(0:1, 0:1); hover_rects(.8, .9, .2, .1, "cell") })
stopifnot(marks[[1L]]$x0 < marks[[1L]]$x1, marks[[1L]]$y0 < marks[[1L]]$y1)
stopifnot(identical(hover_number(c(1234.567, NA, 0, 26.9876, 15.08)), c("1,235", "unavailable", "0", "26.99", "15.08")))

# Depth plot: one read-out per position listing every genotype; NA reads as unavailable.
fixture <- list(variant = data.frame(chrom = "c", pos1 = 12L), bin_ratio = list(disabled = TRUE),
  depth = data.frame(genotype = rep(c("0/0", "0/1"), each = 3L), pos1 = rep(11:13, 2L), mean_depth = c(4, 6, NA, 3, 0, 9)))
marks <- hover_capture(plot_variant_depth(fixture, list(start1 = 11L, end1 = 13L), c("0/0" = "red", "0/1" = "blue")))
stopifnot(length(marks) == 1L, marks[[1L]]$kind == "x",
          identical(unclass(marks[[1L]]$label), c("c:11\n0/0: 4\n0/1: 3", "c:12\n0/0: 6\n0/1: 0", "c:13\n0/0: unavailable\n0/1: 9")))

# Bin ratio: sample read-outs carry the ratio only, never an identifier; clipped
# points report their true value; each panel also has its three dosage means.
rows <- data.frame(site_id = paste0("donor:+:chr1:", 1:2), n = 30L, n_dosage0 = 10L, n_dosage1 = 10L, n_dosage2 = 10L,
  mean_ratio_dosage0 = .9, mean_ratio_dosage1 = .7, mean_ratio_dosage2 = .5, beta = -.2, p_value = 1e-5,
  p_permutation = .01, caution = "", bin_check = "OK", status = "OK", stringsAsFactors = FALSE)
samples <- data.frame(sample = paste0("SECRET", 1:60), site_id = rep(rows$site_id, each = 30L), dosage = rep(rep(0:2, each = 10L), 2L),
  ratio = rep(c(.9, .7, .5), each = 10L) + c(-50, rep(0, 59L)), in_model = TRUE, stringsAsFactors = FALSE)
marks <- hover_capture(plot_bin_ratio_by_dosage(samples, rows))
labels <- unlist(lapply(marks, function(g) unclass(g$label)))
stopifnot(length(marks) == 4L, sum(lengths(lapply(marks, `[[`, "label"))) == 66L, !any(grepl("SECRET", labels, fixed = TRUE)),
          any(grepl("ratio -49.1 (drawn at the axis edge)", labels, fixed = TRUE)),
          any(labels == "donor:+:chr1:2\ndosage 1 · mean ratio 0.7 · n = 10"),
          marks[[2L]]$x[[1L]] < .5, marks[[4L]]$x[[1L]] > .5)

# Junction rate: one read-out per dosage mean.
rate <- data.frame(intron_start1 = 100L, intron_end1 = 200L, strands = "+", status = "OK", log10_p_value = -5, n = 30L,
  n_dosage0 = 10L, n_dosage1 = 10L, n_dosage2 = 10L, mean_rate_dosage0 = 1, mean_rate_dosage1 = 2, mean_rate_dosage2 = 3,
  beta = 1, p_value = 1e-5, p_permutation = .01, caution = "", stringsAsFactors = FALSE)
marks <- hover_capture(plot_junction_rate_by_dosage(rate))
stopifnot(length(marks) == 1L, identical(unclass(marks[[1L]]$label)[[3L]],
  "100-200 (+)\ndosage 2 · mean 3 exact reads per 100 denominator reads · n = 10"))

# Transcript track: a transcript band and its exons; the exon read-out gives coordinates and length.
annotation <- list(metadata = list(release = "48"),
  transcripts = data.frame(gene_name = "G", transcript_id = "T1", strand = "+", start1 = 100L, end1 = 900L,
    transcript_type = "protein_coding", tags = "basic", stringsAsFactors = FALSE),
  exons = data.frame(transcript_id = "T1", start1 = c(100L, 700L), end1 = c(199L, 900L)), features = NULL)
marks <- hover_capture(tryCatch(annotation_plot_models(annotation, c(50, 1000)), error = function(e) stop(e)))
labels <- unlist(lapply(marks, function(g) unclass(g$label)))
stopifnot(all(vapply(marks, `[[`, "", "kind") == "rect"), any(grepl("transcript 100-900", labels, fixed = TRUE)),
          any(grepl("exon 100-199 · 100 bp", labels, fixed = TRUE)), any(grepl("exon 700-900 · 201 bp", labels, fixed = TRUE)))

# The message a plot sends: the namespaced output id and JSON arrays even for one mark.
sent <- NULL
session <- list(ns = shiny::NS("mod"), sendCustomMessage = function(type, message) sent <<- list(type = type, message = message))
hover_send(session, "depth_plot", { graphics::plot(0:1, 0:1); hover_points(.5, .5, "only") })
json <- jsonlite::toJSON(sent$message, auto_unbox = TRUE, digits = NA)
stopifnot(sent$type == "regshiny-hover", sent$message$id == "mod-depth_plot", grepl('"label":["only"]', json, fixed = TRUE),
          grepl('"kind":"xy"', json, fixed = TRUE))
# A plot without marks clears the earlier ones.
hover_send(session, "depth_plot", graphics::plot.new())
stopifnot(length(sent$message$groups) == 0L)
grDevices::dev.off()

# The page loads the script, and the script handles the three kinds of mark.
app <- paste(readLines("app.R"), collapse = "\n"); js <- paste(readLines("www/plot_hover.js"), collapse = "\n")
stopifnot(grepl('plot_hover.js?v=', app, fixed = TRUE), grepl('source("src/plot_hover.R", local = TRUE)', app, fixed = TRUE),
          grepl("regshiny-hover", js, fixed = TRUE), all(vapply(c("'xy'", "'rect'", "'x'"), grepl, logical(1), js, fixed = TRUE)),
          grepl(".regshiny-tip", paste(readLines("www/frontend.css"), collapse = "\n"), fixed = TRUE))
# Every plot output registers its marks and is redrawn, not replayed, on resize.
for (f in list.files("src", "_module\\.R$", full.names = TRUE)) {
  code <- paste(readLines(f), collapse = "\n")
  n_plots <- lengths(regmatches(code, gregexpr("renderPlot(", code, fixed = TRUE)))
  stopifnot(n_plots == lengths(regmatches(code, gregexpr("hover_send(session", code, fixed = TRUE))),
            n_plots == lengths(regmatches(code, gregexpr("execOnResize = TRUE", code, fixed = TRUE))))
}
cat("Plot hover tests passed\n")
