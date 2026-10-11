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
  marks <- as_frame(t$marks)
  marks$color <- jbrowse_color(marks$color)
  rules <- as_frame(t$rules)
  rules$color <- jbrowse_color(ifelse(rules$color == "", "grey50",
    rules$color))
  structure(
    list(bins = as_frame(t$bins), marks = marks,
      sources = as_frame(t$sources), rules = rules, keys = list(),
      domain = unlist(t$domain), scale_type = t$scaleType,
      symlog_constant = t$symlogConstant, rendering = t$rendering,
      origin = t$origin, row_layout = t$rowLayout,
      gap_limit_bp = t$gapLimitBp,
      region = list(chrom = chrom, start = start, end = end)),
    primary = "marks",
    class = c("jbrowse_wiggle", "jbrowse_result"))
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

#' A quantitative track's marks as ggplot layers
#'
#' `geom_jbrowse_wiggle()` draws the marks as the display's plot type draws
#' them: bars from the origin, a stepped or interpolated line with its
#' whisker band, points, or density, each in the display's color. Map
#' `score`, `source` or any other column of the marks table to draw it your
#' own way. `geom_jbrowse_rules()` draws the reference lines `scales$y$rules`
#' declares.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_wiggle()] result, or `NULL` for the plot's
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' bw <- system.file("extdata", "volvox_microarray_posneg.bw",
#'   package = "ggjbrowse")
#' w <- jbrowse_wiggle(bw, "ctgA", 0, 50000)
#' library(ggplot2)
#' ggplot(w) + geom_jbrowse_wiggle(aes(fill = abs(score))) +
#'   scale_fill_viridis_c()
geom_jbrowse_wiggle <- function(mapping = NULL, data = NULL, ...) {
  of <- function(rendering, pick) function(x) {
    d <- pick(x)
    if (x$rendering %in% rendering) d else d[0, , drop = FALSE]
  }
  marks <- function(x) x$marks[!x$marks$band, , drop = FALSE]
  band <- function(x) x$marks[x$marks$band, , drop = FALSE]
  list(
    jbrowse_layer(ggplot2::GeomRect, data, band,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = low, ymax = score,
        fill = color, alpha = I(opacity)), list(...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomRect, data, of("xyplot", marks),
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = low, ymax = score,
        fill = color), list(...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomRect, data, of("density", marks),
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 1,
        fill = color), list(...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomPoint, data, of("scatter", marks),
      for_geom(mapping, ggplot2::GeomPoint),
      ggplot2::aes(x = (start + end) / 2, y = score, colour = color),
      list(size = 0.6, ...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomSegment, data, of("line", marks),
      for_geom(mapping, ggplot2::GeomSegment),
      ggplot2::aes(x = start, xend = end, y = score, yend = score,
        colour = color), list(...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomPath, data,
      of("linecenter", function(x) line_runs(marks(x), x$gap_limit_bp)),
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = (start + end) / 2, y = score, colour = color,
        group = run), list(...), keyed = FALSE))
}

#' @rdname geom_jbrowse_wiggle
#' @export
geom_jbrowse_rules <- function(mapping = NULL, data = NULL, ...) {
  jbrowse_layer(ggplot2::GeomHline, data, function(x) {
    if (x$rendering == "density") x$rules[0, ] else x$rules
  }, mapping, ggplot2::aes(yintercept = value, colour = color),
    list(linetype = "dashed", ...), keyed = FALSE)
}

#' Draw a BigWig as JBrowse's quantitative track
#'
#' [geom_jbrowse_wiggle()] and [geom_jbrowse_rules()] on the score axis the
#' display autoscales to, a row per source where the display gives each one.
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
  density <- object$rendering == "density"
  p <- ggplot2::ggplot(object) +
    geom_jbrowse_wiggle() +
    geom_jbrowse_rules() +
    (if (density) NULL else y_scale(object)) +
    ggplot2::coord_cartesian(
      xlim = c(object$region$start, object$region$end),
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
