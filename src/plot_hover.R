# Hover read-outs for the base-graphics plots. A plot function calls
# hover_points(), hover_columns() or hover_rects() right after drawing a mark;
# outside hover_capture() these do nothing, so the functions still serve
# downloads and tests unchanged. Positions are stored as fractions of the whole
# device (0-1, origin bottom left), converted while the panel is current, so
# multi-panel plots work and the browser needs no knowledge of the axes.
# www/plot_hover.js shows the label of the mark under the pointer.
# Display only: nothing here feeds a statistic, an export or a saved result.
.hover_state <- new.env(parent = emptyenv())
.hover_state$groups <- NULL

hover_recording <- function() is.list(.hover_state$groups)

hover_number <- function(x, digits = 4L) {
  out <- rep("unavailable", length(x)); ok <- is.finite(x)
  out[ok] <- format(signif(x[ok], digits), trim = TRUE, big.mark = ",")
  out
}
hover_position <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)

hover_add <- function(group) {
  .hover_state$groups[[length(.hover_state$groups) + 1L]] <- group
  invisible(NULL)
}
hover_panel_box <- function() {
  usr <- graphics::par("usr")
  c(sort(graphics::grconvertX(usr[1:2], "user", "ndc")), sort(graphics::grconvertY(usr[3:4], "user", "ndc")))
}

# Marks hit when the pointer is within a few pixels: points, arc samples.
hover_points <- function(x, y, label) {
  if (!hover_recording() || !length(x)) return(invisible(NULL))
  label <- rep_len(as.character(label), length(x)); keep <- is.finite(x) & is.finite(y)
  if (!any(keep)) return(invisible(NULL))
  hover_add(list(kind = "xy", x = I(graphics::grconvertX(x[keep], "user", "ndc")),
    y = I(graphics::grconvertY(y[keep], "user", "ndc")), label = I(label[keep])))
}

# One label per x position, hit anywhere inside the current panel at the
# nearest x: curves and per-base tracks. More positions than max_n are thinned
# evenly (each label names its own position, so a thinned read-out is still exact).
hover_columns <- function(x, label, max_n = 3000L) {
  if (!hover_recording() || !length(x)) return(invisible(NULL))
  label <- as.character(label)
  if (length(x) > max_n) { i <- unique(round(seq(1, length(x), length.out = max_n))); x <- x[i]; label <- label[i] }
  hover_add(list(kind = "x", box = I(hover_panel_box()), x = I(graphics::grconvertX(x, "user", "ndc")), label = I(label)))
}

# Marks hit when the pointer is inside them: cells, bars, exons. The smallest
# rectangle under the pointer wins.
hover_rects <- function(x0, y0, x1, y1, label) {
  if (!hover_recording() || !length(x0)) return(invisible(NULL))
  n <- max(length(x0), length(y0), length(x1), length(y1))
  x0 <- rep_len(x0, n); x1 <- rep_len(x1, n); y0 <- rep_len(y0, n); y1 <- rep_len(y1, n)
  label <- rep_len(as.character(label), n); keep <- is.finite(x0) & is.finite(x1) & is.finite(y0) & is.finite(y1)
  if (!any(keep)) return(invisible(NULL))
  ax <- graphics::grconvertX(x0[keep], "user", "ndc"); bx <- graphics::grconvertX(x1[keep], "user", "ndc")
  ay <- graphics::grconvertY(y0[keep], "user", "ndc"); by <- graphics::grconvertY(y1[keep], "user", "ndc")
  hover_add(list(kind = "rect", x0 = I(pmin(ax, bx)), x1 = I(pmax(ax, bx)), y0 = I(pmin(ay, by)), y1 = I(pmax(ay, by)),
    label = I(label[keep])))
}

# Points along a junction arc (parabola of the given height between two ends).
hover_arc <- function(start, end, height, label, n = 15L) {
  t <- seq(0, 1, length.out = n)
  hover_points(start + t * (end - start), 4 * height * t * (1 - t), label)
}

# Runs plotting code and returns the marks it registered.
hover_capture <- function(expr) {
  previous <- .hover_state$groups
  .hover_state$groups <- list()
  on.exit(.hover_state$groups <- previous)
  force(expr)
  .hover_state$groups
}

# Inside renderPlot(..., execOnResize = TRUE): draws the plot and sends its
# marks to the browser. A replayed plot would keep marks of the old size, hence
# execOnResize. Write the plotting call at the call site (not through a
# wrapper function), so it is evaluated again on every render.
hover_send <- function(session, id, expr) {
  groups <- hover_capture(expr)
  session$sendCustomMessage("regshiny-hover", list(id = session$ns(id), groups = groups))
  invisible(NULL)
}
