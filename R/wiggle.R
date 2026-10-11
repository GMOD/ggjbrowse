#' Read a region of a BigWig as JBrowse draws it
#'
#' Runs JBrowse's own quantitative display in V8: the scores it fetches at
#' this zoom, the marks its encoder lays out from them, the color each paints
#' in and the score axis it autoscales to. Coordinates are 0-based,
#' half-open, on the reference.
#'
#' @param file a BigWig, a path or a URL
#' @param chrom,start,end the region, 0-based half-open
#' @param display the display's settings, as a track config's `displays`
#'   entry writes them for a `LinearWiggleDisplay`: a list with any of `mark`
#'   (`"bar"`, `"line"`, `"point"`, or `"span"` for density), `interpolate`
#'   (`"step"` or `"linear"`), `aggregate` (`"mean"`, `"min"`, `"max"`),
#'   `extent` (`"min-max"` or `"none"`), `origin`, `color` (a CSS color, or a
#'   list of the color object's slots) and `scales` (`list(y = list(type =,
#'   domainMin =, domainMax =, domainQuantile =, zero =, rules =))`)
#' @param width the width in pixels the display draws the region across,
#'   which picks the BigWig zoom level that answers and so each bin's width
#' @return a `jbrowse_wiggle` list:
#'   * `bins`: each bin's scores as the file serves them at this zoom:
#'     `source`, `start`, `end`, `score`, and `min` and `max` where the file
#'     summarizes the bin
#'   * `marks`: what the display draws from them, in paint order: `layer`,
#'     `source`, `row`, `start`, `end`, `score`, `low` (where the mark reaches
#'     down to: the origin, or a whisker band's min), `band`, `color` and
#'     `opacity`; under a min-max extent the max and min bands come back as
#'     layers of their own, tinted as the display tints them
#'   * `sources`: `name`, `label`, `group`, `color`
#'   * `rules`: the reference lines `scales$y$rules` draws
#'
#'   `domain` is the score axis the display autoscales to, `scale_type` its
#'   type, and `rendering` what the display draws (`"xyplot"`, `"line"`,
#'   `"linecenter"`, `"scatter"`, `"density"`).
#' @export
#' @examples
#' bw <- system.file("extdata", "volvox_microarray_posneg.bw",
#'   package = "ggjbrowse")
#' w <- jbrowse_wiggle(bw, "ctgA", 0, 50000)
#' head(w$marks)
jbrowse_wiggle <- function(file, chrom, start, end, display = NULL,
                           width = 800) {
  args <- list(
    adapterConfig = list(type = "BigWigAdapter",
      bigWigLocation = location(file)),
    region = list(refName = chrom, start = start, end = end),
    bpPerPx = (end - start) / width)
  if (!is.null(display)) {
    if (!is.null(display$color) && is.list(display$color)) {
      display$color <- color_object(display$color)
    }
    args$display <- display
  }
  t <- call_engine("wiggle", args)
  structure(
    list(bins = as_frame(t$bins), marks = as_frame(t$marks),
      sources = as_frame(t$sources), rules = as_frame(t$rules),
      domain = unlist(t$domain), scale_type = t$scaleType,
      symlog_constant = t$symlogConstant, rendering = t$rendering,
      origin = t$origin, row_layout = t$rowLayout,
      gap_limit_bp = t$gapLimitBp,
      region = list(chrom = chrom, start = start, end = end)),
    class = "jbrowse_wiggle")
}

#' @export
print.jbrowse_wiggle <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_wiggle> %s:%.0f-%.0f, %d bins from %d source(s), drawn as %s on %s axis %g to %g\n",
    x$region$chrom, x$region$start, x$region$end, nrow(x$bins),
    nrow(x$sources), x$rendering, x$scale_type, x$domain[1], x$domain[2]))
  invisible(x)
}

# The display's symlog axis, d3's sign(x) log(1 + |x| / c)
symlog_transform <- function(constant) {
  scales::trans_new("symlog",
    function(x) sign(x) * log1p(abs(x) / constant),
    function(y) sign(y) * constant * expm1(abs(y)))
}

y_scale <- function(x) {
  switch(x$scale_type,
    log = ggplot2::scale_y_log10(),
    symlog = ggplot2::scale_y_continuous(
      transform = symlog_transform(x$symlog_constant)),
    ggplot2::scale_y_continuous())
}

#' Draw a BigWig as JBrowse's quantitative track
#'
#' The marks the display draws, in its colors and paint order: bars from the
#' origin, a stepped or interpolated line with its whisker band, points, or
#' density, on the score axis the display autoscales to, with its reference
#' lines.
#'
#' @param object a [jbrowse_wiggle()] result
#' @param ... unused
#' @return a ggplot, to extend or to stack over an alignments track with
#'   patchwork
#' @exportS3Method ggplot2::autoplot
#' @examples
#' bw <- system.file("extdata", "volvox_microarray_posneg.bw",
#'   package = "ggjbrowse")
#' autoplot(jbrowse_wiggle(bw, "ctgA", 0, 50000))
autoplot.jbrowse_wiggle <- function(object, ...) {
  m <- object$marks
  xlim <- c(object$region$start, object$region$end)
  density <- object$rendering == "density"
  band <- m[m$band, , drop = FALSE]
  marks <- m[!m$band, , drop = FALSE]
  rules <- object$rules
  rules$color[rules$color == ""] <- "grey50"
  drawn <- switch(object$rendering,
    density = ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end,
      ymin = 0, ymax = 1, fill = color), marks),
    scatter = ggplot2::geom_point(ggplot2::aes(x = (start + end) / 2,
      y = score, colour = color), marks, size = 0.6),
    line = ggplot2::geom_segment(ggplot2::aes(x = start, xend = end,
      y = score, yend = score, colour = color), marks),
    linecenter = ggplot2::geom_path(ggplot2::aes(x = (start + end) / 2,
      y = score, colour = color, group = run),
      line_runs(marks, object$gap_limit_bp)),
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = low,
      ymax = score, fill = color), marks))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = low,
      ymax = score, fill = color, alpha = opacity), band) +
    drawn
  if (!density) {
    p <- p +
      ggplot2::geom_hline(ggplot2::aes(yintercept = value, colour = color),
        rules, linetype = "dashed") +
      y_scale(object)
  }
  p <- p +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::coord_cartesian(xlim = xlim,
      ylim = if (density) c(0, 1) else object$domain, expand = FALSE) +
    ggplot2::labs(x = sprintf("%s (bp)", object$region$chrom), y = NULL) +
    ggplot2::theme_minimal()
  if (object$row_layout) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(source))
  }
  if (density) p + no_y else p
}

# An interpolated line's points split into runs where the display breaks it:
# across a hole wider than the gap limit
line_runs <- function(marks, limit) {
  mid <- (marks$start + marks$end) / 2
  gap <- c(FALSE, diff(mid) > (limit %||% Inf) |
    utils::head(marks$source, -1) != utils::tail(marks$source, -1))
  marks$run <- cumsum(gap)
  marks
}
