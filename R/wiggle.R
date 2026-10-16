#' Read a region of a BigWig as JBrowse draws it
#'
#' Runs JBrowse's own quantitative display in V8: the scores it fetches at
#' this zoom, the marks its encoder lays out from them, the color each paints
#' in and the score axis it autoscales to. Coordinates are 0-based,
#' half-open, on the reference.
#'
#' @param file a BigWig, a path or a URL; or a JBrowse track's `adapter`
#'   config, as a list, such as a `MultiWiggleAdapter`'s
#' @param region the regions to read, read as one view and drawn end to end
#'   on one axis as a JBrowse view lays them out: loci as JBrowse's location
#'   box writes them (`"ctgA:1,001-9,000"`, 1-based and inclusive), contig
#'   names alone for whole contigs, a data frame of `chrom`, `start` and `end`
#'   (0-based, half-open, as BED writes them), or a GRanges. With `NULL`, the
#'   default, the reader reads the whole file: each contig it lists but those
#'   its index shows empty, to the length the file or `assembly` states. Several
#'   regions share one score axis and color ramp
#' @param display the display's settings, as a track config's `displays`
#'   entry writes them for a `LinearWiggleDisplay`: a list with any of `mark`
#'   (`"bar"`, `"line"`, `"point"`, or `"span"` for density), `interpolate`
#'   (`"step"` or `"linear"`), `aggregate` (`"mean"`, `"min"`, `"max"`),
#'   `extent` (`"min-max"` or `"none"`), `origin`, `color` (a CSS color, or a
#'   list of the color object's slots) and `scales` (`list(y = list(type =,
#'   domainMin =, domainMax =, domainQuantile =, zero =, rules =))`)
#' @param width the width in pixels the display draws the region across,
#'   which picks the BigWig zoom level that answers and so each bin's width
#' @param assembly a [jbrowse_assembly()], an assembly as a JBrowse config's
#'   `assemblies` entry writes it, or its FASTA or 2bit alone: `region` may
#'   name a contig by any alias it lists, and the reader names the contig in
#'   its tables as `region` does
#' @return a `jbrowse_wiggle` list, beside the slots [jbrowse_result] lists:
#'   * `bins`: each bin's scores as the file serves them at this zoom:
#'     `source`, `start`, `end`, `score`, and `min` and `max` where the file
#'     summarizes the bin
#'   * `marks`: what the display draws from them, in paint order: `layer`,
#'     `source`, `row`, `start`, `end`, `score`, `low` (where the mark reaches
#'     down to: the origin, or a whisker band's min), `band`, `color` and
#'     `opacity`; under a min-max extent the max and min bands come back as
#'     layers of their own, tinted as the display tints them
#'   * `sources`: `name`, `label`, `group`, `color`; `bins` and `marks` name
#'     their source by `name`, and assigning `sources` drops the bins and
#'     marks of the sources it leaves out
#'   * `rules`: the reference lines `scales$y$rules` draws
#'
#'   `y_scale` is the score axis the display autoscales to, `ramps` the
#'   score's color bar where one ramp paints every row (density, or a
#'   gradient), and `rendering` what the display draws (`"xyplot"`, `"line"`,
#'   `"linecenter"`, `"scatter"`, `"density"`).
#' @export
#' @examples
#' bw <- system.file("extdata", "gapdh_phylop100way.bw",
#'   package = "ggjbrowse")
#' w <- jbrowse_wiggle(bw, "chr12:6,534,201-6,538,700")
#' head(w$marks)
jbrowse_wiggle <- function(file, region = NULL, display = NULL,
                           width = DEFAULT_WIDTH, assembly = NULL) {
  adapter <- adapter_of(file, function()
    list(type = "BigWigAdapter", bigWigLocation = location(file)))
  assembly <- assembly_config(assembly)
  regions <- reader_regions(adapter, region, assembly)
  args <- list(adapterConfig = adapter, regions = engine_regions(regions),
    bpPerPx = sum(regions$end - regions$start) / width)
  args$assembly <- assembly
  if (isTRUE(attr(regions, "whole"))) args$rowLimit <- WHOLE_FILE_ROWS
  args$display <- display_config(display)
  t <- call_engine("wiggle", args)
  marks <- as_frame(t$marks)
  bins <- as_frame(t$bins)
  rules <- as_frame(t$rules)
  in_one_call(structure(
    list(bins = bins, marks = marks, sources = as_frame(t$sources),
      rules = rules, display = t$config, keys = list(),
      ramps = ramp_frames(t$ramps), y_scale = y_scale(t), fields = list(),
      rendering = t$rendering, origin = t$origin, row_layout = t$rowLayout,
      height = t$height, bp_per_px = args$bpPerPx),
    primary = "marks",
    class = c("jbrowse_wiggle", "jbrowse_result")),
    regions, t, width)
}

#' A quantitative track's sources as the display's clustering reads them
#'
#' The matrix the multi-row quantitative display's "Cluster rows" clusters: a
#' row per source of `x$sources`, in its order, named by its `name`, and a
#' column per pixel of the `width` the track was read at, holding the mean of
#' the display's `aggregate` over the bins covering that pixel, 0 where none
#' does. `hclust(dist(as.matrix(w)), "average")` is the R script the
#' display's clustering dialog writes.
#' @param x a [jbrowse_wiggle()] result
#' @param ... unused
#' @return a numeric matrix, its columns named by the locus each pixel starts
#'   at
#' @export
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' cells <- c("gm12878", "k562", "hepg2")
#' sources <- list(type = "MultiWiggleAdapter", subadapters = lapply(cells,
#'   function(cell) list(type = "BigWigAdapter", name = toupper(cell),
#'     bigWigLocation = list(uri = e(sprintf("hbb_dnase_%s.bw", cell))))))
#' w <- jbrowse_wiggle(sources, "chr11:5,270,001-5,320,000",
#'   display = list(mark = "span"))
#' dim(as.matrix(w))
#' w$tree <- hclust(dist(as.matrix(w)), "average")
as.matrix.jbrowse_wiggle <- function(x, ...) {
  r <- x$regions
  per_px <- 1 / x$bp_per_px
  columns <- pmax(0, floor((r$end - r$start) * per_px))
  offset <- cumsum(c(0, columns))[seq_along(columns)]
  b <- x$bins[x$bins$source %in% x$sources$name, , drop = FALSE]
  value <- switch(x$display$aggregate, min = b$min, max = b$max,
    b$score) %||% b$score
  i <- b$region
  from <- pmax(0, trunc((b$start - r$start[i]) * per_px))
  to <- pmin(columns[i], pmax(trunc((b$end - r$start[i]) * per_px), from + 1))
  keep <- b$end > r$start[i] & to > from & !is.na(value)
  n <- (to - from)[keep]
  row <- rep(match(b$source[keep], x$sources$name), n)
  column <- rep(offset[i][keep] + from[keep], n) + sequence(n)
  cell <- (column - 1) * nrow(x$sources) + row
  sums <- rowsum(rep(value[keep], n), cell)
  counts <- rowsum(rep(1, length(cell)), cell)
  m <- matrix(0, nrow(x$sources), sum(columns),
    dimnames = list(x$sources$name, unlist(lapply(seq_len(nrow(r)),
      function(k) sprintf("%s:%.0f", r$chrom[k],
        r$start[k] + floor((seq_len(columns[k]) - 1) * x$bp_per_px) + 1)))))
  m[as.numeric(rownames(sums))] <- sums / counts
  m
}

#' @export
print.jbrowse_wiggle <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_wiggle> %s, %d bins from %d source(s), drawn as %s on %s axis %g to %g\n",
    print_loci(x), nrow(x$bins),
    nrow(x$sources), x$rendering, x$y_scale$type, x$y_scale$domain[1],
    x$y_scale$domain[2]))
  invisible(x)
}

# The display's symlog axis, d3's sign(x) log(1 + |x| / c)
symlog_transform <- function(constant) {
  scales::trans_new("symlog",
    function(x) sign(x) * log1p(abs(x) / constant),
    function(y) sign(y) * constant * expm1(abs(y)))
}

#' The display's score axis as a ggplot scale
#'
#' The y scale a quantitative display plots its scores on: linear, log, or
#' symlog around the display's constant.
#' @param object a [jbrowse_wiggle()] result
#' @param ... passed to the ggplot scale, such as `name`
#' @return a ggplot scale
#' @export
#' @examples
#' bw <- system.file("extdata", "gapdh_phylop100way.bw",
#'   package = "ggjbrowse")
#' w <- jbrowse_wiggle(bw, "chr12:6,534,201-6,538,700",
#'   display = list(scales = list(y = list(type = "symlog"))))
#' library(ggplot2)
#' ggplot(w) + geom_jbrowse_wiggle() + scale_y_jbrowse(w)
scale_y_jbrowse <- function(object, ...) {
  s <- object$y_scale
  switch(s$type,
    log = ggplot2::scale_y_log10(...),
    symlog = scale_with(ggplot2::scale_y_continuous,
      list(transform = symlog_transform(s$symlog_constant)), ...),
    ggplot2::scale_y_continuous(...))
}

# Marks as the display stacks its rows: under a row layout each source's
# scores fill its own slice of the axis, first source on top, with no gap
# between rows, on the display's own scale; `base` is the slice's floor
# The display draws a mark narrower than a pixel a pixel wide, so a sparse
# track of 1 bp sites stays visible zoomed out.
pixel_wide <- function(x, keep) {
  m <- x$marks[keep, , drop = FALSE]
  if (!is.null(x$bp_per_px)) m$end <- pmax(m$end, m$start + x$bp_per_px)
  m
}

in_rows <- function(x, d) {
  n <- nrow(x$sources)
  d$base <- if (isTRUE(x$row_layout)) n - 1 - d$row else rep(0, nrow(d))
  if (isTRUE(x$row_layout)) {
    t <- switch(x$y_scale$type, log = scales::transform_log10(),
      symlog = symlog_transform(x$y_scale$symlog_constant),
      scales::transform_identity())
    lim <- t$transform(x$y_scale$domain)
    at <- function(v) {
      d$base + pmin(1, pmax(0, (t$transform(v) - lim[1]) / diff(lim)))
    }
    d$low <- at(d$low)
    d$score <- at(d$score)
  }
  d
}

# A line segment that changes color across the origin draws each side in its
# own color, as the display's shader does: a point at the crossing takes the
# next point's color, since a path segment takes its start point's.
cut_at_origin <- function(d, origin) {
  n <- nrow(d)
  if (n < 2) return(d)
  origin <- origin %||% 0
  i <- which(d$run[-1] == d$run[-n] & d$color[-1] != d$color[-n] &
    (d$score[-n] - origin) * (d$score[-1] - origin) < 0)
  if (!length(i)) return(d)
  cut <- d[i + 1, , drop = FALSE]
  cut$x <- d$x[i] + (origin - d$score[i]) / (d$score[i + 1] - d$score[i]) *
    (d$x[i + 1] - d$x[i])
  cut$score <- origin
  rbind(d, cut)[order(c(seq_len(n), i + 0.5)), , drop = FALSE]
}

#' A quantitative track's marks as ggplot layers
#'
#' `geom_jbrowse_wiggle()` draws the marks as the display's plot type draws
#' them: bars from the origin, a stepped or interpolated line with its
#' whisker band, points, or density, each in the display's color. Map
#' `score`, `source` or any other column of the marks table to draw it your
#' own way. A line breaks into the `run`s the display breaks it into.
#' `geom_jbrowse_rules()` draws the reference lines `scales$y$rules`
#' declares, each labelled above its line as the display labels it.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_wiggle()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' bw <- system.file("extdata", "gapdh_phylop100way.bw",
#'   package = "ggjbrowse")
#' w <- jbrowse_wiggle(bw, "chr12:6,534,201-6,538,700")
#' library(ggplot2)
#' ggplot(w) + geom_jbrowse_wiggle(aes(fill = abs(score))) +
#'   scale_fill_viridis_c()
geom_jbrowse_wiggle <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "marks", "geom_jbrowse_wiggle")
  of <- function(rendering, pick) function(x) {
    d <- pick(x)
    if (x$rendering %in% rendering) d else d[0, , drop = FALSE]
  }
  marks <- function(x) in_rows(x, pixel_wide(x, !x$marks$band))
  band <- function(x) in_rows(x, pixel_wide(x, x$marks$band))
  steps <- function(x) {
    m <- marks(x)
    d <- m[rep(seq_len(nrow(m)), each = 2), , drop = FALSE]
    d$x <- as.vector(rbind(m$start, m$end))
    cut_at_origin(d, x$origin)
  }
  centers <- function(x) {
    d <- marks(x)
    d$x <- (d$start + d$end) / 2
    cut_at_origin(d, x$origin)
  }
  list(
    jbrowse_layer(ggplot2::GeomRect, data, band,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = low, ymax = score,
        fill = I(color), alpha = I(opacity)), list(...)),
    jbrowse_layer(ggplot2::GeomRect, data, of("xyplot", marks),
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = low, ymax = score,
        fill = I(color)), list(...)),
    jbrowse_layer(ggplot2::GeomRect, data, of("density", marks),
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = base, ymax = base + 1,
        fill = I(color)), list(...)),
    jbrowse_layer(ggplot2::GeomPoint, data, of("scatter", marks),
      for_geom(mapping, ggplot2::GeomPoint),
      ggplot2::aes(x = (start + end) / 2, y = score, colour = I(color)),
      list(size = 0.6, ...)),
    jbrowse_layer(ggplot2::GeomPath, data, of("line", steps),
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = score, colour = I(color), group = run),
      list(...)),
    jbrowse_layer(ggplot2::GeomPath, data, of("linecenter", centers),
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = score, colour = I(color), group = run),
      list(...)))
}

#' @rdname geom_jbrowse_wiggle
#' @export
geom_jbrowse_rules <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "rules", "geom_jbrowse_rules")
  rules <- function(x) {
    r <- if (x$rendering == "density" || isTRUE(x$row_layout)) x$rules[0, ]
      else x$rules
    r$x <- rep(xlim_of(x)[1], nrow(r))
    r
  }
  list(
    jbrowse_layer(ggplot2::GeomHline, data, rules,
      for_geom(mapping, ggplot2::GeomHline),
      ggplot2::aes(yintercept = value, colour = I(color)),
      list(linetype = "dashed", ...)),
    jbrowse_layer(ggplot2::GeomText, data,
      function(x) {
        r <- rules(x)
        r[r$label != "", , drop = FALSE]
      },
      for_geom(mapping, ggplot2::GeomText),
      ggplot2::aes(x = x, y = value, label = label, colour = I(color)),
      list(hjust = 0, vjust = -0.4, size = 2.6, ...)))
}

#' Draw a BigWig as JBrowse's quantitative track
#'
#' [geom_jbrowse_wiggle()] and [geom_jbrowse_rules()] on the score axis the
#' display autoscales to, a row per source where the display gives each one.
#'
#' @param object a [jbrowse_wiggle()] result
#' @param width the pixels the figure draws the regions across, which
#'   [jbrowse_width()] sizes a stack by; `NULL` takes the `width` the track
#'   was read at. Narrower than that warns, since the file's zoom level was
#'   picked for more pixels
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot, to extend or to stack over an alignments track with
#'   patchwork
#' @exportS3Method ggplot2::autoplot
#' @examples
#' bw <- system.file("extdata", "gapdh_phylop100way.bw",
#'   package = "ggjbrowse")
#' autoplot(jbrowse_wiggle(bw, "chr12:6,534,201-6,538,700"))
autoplot.jbrowse_wiggle <- function(object, width = NULL, ...) {
  rlang::check_dots_empty()
  object <- on_axis(drawn_at(object, width))
  density <- object$rendering == "density"
  rows <- isTRUE(object$row_layout)
  tree <- rows && !is.null(object$tree)
  n <- nrow(object$sources)
  p <- ggplot2::ggplot(object) +
    geom_jbrowse_wiggle() +
    geom_jbrowse_rules() +
    guides_jbrowse(object) +
    (if (!rows && !density) scale_y_jbrowse(object)) +
    ggplot2::coord_cartesian(
      xlim = xlim_of(object),
      ylim = if (rows) c(0, n) else if (density) c(0, 1)
        else object$y_scale$domain,
      expand = FALSE, clip = if (rows) "off" else "on") +
    genomic_x(object) +
    ggplot2::labs(x = axis_label(object), y = NULL) +
    ggplot2::theme_minimal() +
    (if (tree) geom_jbrowse_tree(gap = 3.2, span = TREE_PT)) +
    (if (rows) row_labels(object, if (tree) TREE_PT + 5.2 else 2.2))
  if (density && !rows) p + no_y else p
}

# A row layout's y axis: each source named at its row, `gap_pt` from the
# panel, and none where a row is under 3 of the display's pixels
row_labels <- function(x, gap_pt = 2.2) {
  n <- nrow(x$sources)
  y <- n - seq_len(n) + 0.5
  list(ggplot2::scale_y_continuous(breaks = y, labels = x$sources$label),
    ggplot2::theme(panel.grid = ggplot2::element_blank()),
    if (x$height / n >= 3) row_names(x$sources$label, y, gap_pt)
    else if (gap_pt > 2.2) list(no_y, axis_gap(gap_pt))
    else no_y)
}
