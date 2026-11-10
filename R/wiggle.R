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
#'   (`"bar"`, `"line"`, `"point"`, or `"rect"` for density, or
#'   `list(type =, interpolate =, extent =)` with `"step"` or `"linear"` for a
#'   line and `"min-max"` or `"none"` for the band around each bin), `origin`
#'   and `encoding`: its `color`
#'   (`list(value =)` for one CSS color, or the color object's slots) and its
#'   `y`, `list(field =, scale = list(type =, domainMin =, domainMax =,
#'   fenceQuantile =, zero =), rules =)`, whose `field` (`"score"`,
#'   `"minScore"` or `"maxScore"`) is which of a zoom bin's stored mean,
#'   minimum and maximum is plotted
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
#'     layers of their own, tinted as the display tints them; `clip_side` is
#'     1 on a bar whose score passes the axis's top, -1 its bottom, else 0,
#'     where the display paints `clip_strip` (its `color`, `px` tall) across
#'     the row's edge
#'   * `sources`: `name`, `label`, `group`, `color`; `bins` and `marks` name
#'     their source by `name`, and assigning `sources` drops the bins and
#'     marks of the sources it leaves out
#'   * `rules`: the reference lines `encoding$y$rules` declares, as the
#'     display draws them: none under density, which maps the score to color
#'
#'   `y_scale` is the score axis the display autoscales to, `ramps` the
#'   score's color bar where one ramp paints every row (density, or a
#'   gradient), `mark` the display's mark (`"bar"`, `"line"`, `"point"`,
#'   `"rect"`), `interpolate` how a line joins its scores (`"step"` or
#'   `"linear"`), and `summary_field` the bin column the display plots
#'   (`"score"`, `"minScore"` or `"maxScore"`).
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
  bp_per_px <- sum(regions$end - regions$start) / width
  t <- read_regions("wiggle", adapter, regions, assembly, display,
    bpPerPx = bp_per_px)
  marks <- table_of(t, "marks")
  names(marks)[names(marks) == "clipSide"] <- "clip_side"
  bins <- table_of(t, "bins")
  rules <- table_of(t, "rules")
  in_one_call(structure(
    list(bins = bins, marks = marks, sources = table_of(t, "sources"),
      rules = rules, display = t$config,
      keys = list(scores = engine_key(t$legend, t$legendTitle)),
      ramps = ramp_frames(t$ramps), y_scale = y_scale(t$scoreScale),
      fields = list(), mark = t$mark, interpolate = t$interpolate,
      summary_field = t$summaryField, origin = t$origin,
      row_layout = t$rowLayout, clip_strip = t$clipStrip,
      height = t$height, bp_per_px = bp_per_px),
    primary = "marks",
    class = c("jbrowse_wiggle", "jbrowse_result")),
    regions, t, width)
}

#' A quantitative track's sources as the display's clustering reads them
#'
#' The matrix the multi-row quantitative display's "Cluster rows" clusters: a
#' row per source of `x$sources`, in its order, named by its `name`, and a
#' column per pixel of the `width` the track was read at, holding the mean of
#' the bin summary the display's `y` names over the bins covering that pixel,
#' `NA` where none does. `hclust(dist(as.matrix(w)), "average")` is the R script the
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
#'   display = list(mark = "rect"))
#' dim(as.matrix(w))
#' w$tree <- hclust(dist(as.matrix(w)), "average")
as.matrix.jbrowse_wiggle <- function(x, ...) {
  r <- x$regions
  bp <- x$bp_per_px
  columns <- pmax(0, column_floor((r$end - r$start) / bp))
  offset <- cumsum(c(0, columns))[seq_along(columns)]
  b <- x$bins[x$bins$source %in% x$sources$name, , drop = FALSE]
  value <- as.numeric(switch(x$summary_field, minScore = b$min,
    maxScore = b$max, b$score))
  i <- as.integer(b$region)
  from <- pmax(0, column_floor((b$start - r$start[i]) / bp))
  to <- pmin(columns[i], pmax(column_floor((b$end - r$start[i]) / bp),
    from + 1))
  keep <- b$end > r$start[i] & to > from & !is.na(value)
  n <- (to - from)[keep]
  row <- rep(match(b$source[keep], x$sources$name), n)
  column <- rep(offset[i][keep] + from[keep], n) + sequence(n)
  cell <- (column - 1) * nrow(x$sources) + row
  m <- matrix(NA_real_, nrow(x$sources), sum(columns),
    dimnames = list(x$sources$name, unlist(lapply(seq_len(nrow(r)),
      function(k) sprintf("%s:%.0f", r$chrom[k],
        r$start[k] + floor((seq_len(columns[k]) - 1) * bp) + 1)))))
  if (length(cell)) {
    sums <- rowsum(rep(value[keep], n), cell)
    counts <- rowsum(rep(1, length(cell)), cell)
    m[as.numeric(rownames(sums))] <- sums / counts
  }
  m
}

# 1,000 bp over 228 columns is 227.99999999999997 columns
column_floor <- function(x) floor(x + 1e-9 * pmax(1, abs(x)))

#' @export
print.jbrowse_wiggle <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_wiggle> %s, %d bins from %d source(s), drawn as %s on %s axis %g to %g\n",
    loci_text(x$regions), nrow(x$bins),
    nrow(x$sources), x$mark, x$y_scale$type, x$y_scale$domain[1],
    x$y_scale$domain[2]))
  invisible(x)
}

# A display's score or depth axis as a transform: log, linear, or d3's symlog,
# sign(x) log(1 + |x| / c), which scales' asinh-based pseudo-log is not
score_transform <- function(scale) {
  c <- scale$symlog_constant
  switch(scale$type,
    log = scales::transform_log10(),
    symlog = scales::new_transform("symlog",
      function(x) sign(x) * log1p(abs(x) / c),
      function(y) sign(y) * c * expm1(abs(y))),
    scales::transform_identity())
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
#'   display = list(encoding = list(y = list(scale = list(type = "symlog")))))
#' library(ggplot2)
#' ggplot(w) + geom_jbrowse_wiggle() + scale_y_jbrowse(w)
scale_y_jbrowse <- function(object, ...) {
  scale_with(ggplot2::scale_y_continuous,
    list(transform = score_transform(object$y_scale)), ...)
}

# The display draws a mark narrower than a pixel a pixel wide, so a sparse
# track of 1 bp sites stays visible zoomed out.
pixel_wide <- function(x, keep) {
  m <- x$marks[keep, , drop = FALSE]
  m$end <- pmax(m$end, m$start + sum(x$regions$end - x$regions$start) / x$width)
  m
}

# Marks as the display stacks its rows: under a row layout each source's
# scores fill its own slice of the axis, first source on top, with no gap
# between rows, on the display's own scale; `base` is the slice's floor
in_rows <- function(x, d) {
  n <- nrow(x$sources)
  if (x$row_layout) {
    d$base <- n - 1 - d$row
    t <- score_transform(x$y_scale)
    lim <- t$transform(x$y_scale$domain)
    at <- function(v) {
      d$base + pmin(1, pmax(0, (t$transform(v) - lim[1]) / diff(lim)))
    }
    d$low <- at(d$low)
    d$score <- at(d$score)
  } else {
    d$base <- rep(0, nrow(d))
  }
  d
}

# Rules as a row layout draws them: each one inside the shared domain, in
# every row's slice of the axis
in_each_row <- function(x, rules) {
  lim <- x$y_scale$domain
  rules <- rules[rules$value >= lim[1] & rules$value <= lim[2], , drop = FALSE]
  n <- nrow(x$sources)
  d <- rules[rep(seq_len(nrow(rules)), each = n), , drop = FALSE]
  d$row <- rep(seq_len(n) - 1, nrow(rules))
  d$value <- in_rows(x, transform(d, low = value, score = value))$score
  rownames(d) <- NULL
  d[names(d) != "row"]
}

# The strip the display paints across a row's edge over each bar its axis
# clips, `clip_strip$px` of the row's height, on the axis's own scale
clip_strips <- function(x) {
  d <- pixel_wide(x, !x$marks$band & x$marks$clip_side != 0)
  t <- score_transform(x$y_scale)
  lim <- t$transform(x$y_scale$domain)
  rows <- if (x$row_layout) nrow(x$sources) else 1
  frac <- x$clip_strip$px * rows / x$height
  edge <- function(at) t$inverse(lim[1] + at * diff(lim))
  top <- d$clip_side > 0
  d$low <- ifelse(top, edge(1 - frac), x$y_scale$domain[1])
  d$score <- ifelse(top, x$y_scale$domain[2], edge(frac))
  d <- in_rows(x, d)
  d$color <- rep(x$clip_strip$color, nrow(d))
  d
}

# A line segment that changes color across the origin draws each side in its
# own color, as the display's shader does: a point at the crossing takes the
# next point's color, since a path segment takes its start point's.
cut_at_origin <- function(d, origin) {
  n <- nrow(d)
  if (n < 2) return(d)
  i <- which(d$run[-1] == d$run[-n] & d$color[-1] != d$color[-n] &
    (d$score[-n] - origin) * (d$score[-1] - origin) < 0)
  if (!length(i)) return(d)
  cut <- d[i + 1, , drop = FALSE]
  cut$x <- d$x[i] + (origin - d$score[i]) / (d$score[i + 1] - d$score[i]) *
    (d$x[i + 1] - d$x[i])
  cut$score <- origin
  rbind(d, cut)[order(c(seq_len(n), i + 0.5)), , drop = FALSE]
}

# Abutting bars of one style as one area, as overlapping bars paint it, so no
# antialiased seam shows between bars narrower than a device pixel
GeomJBrowseBars <- ggplot2::ggproto("GeomJBrowseBars", ggplot2::GeomRect,
  draw_panel = function(self, data, panel_params, coord, lineend = "butt",
                        linejoin = "mitre") {
    data <- vctrs::vec_slice(data,
      stats::complete.cases(data[c("xmin", "xmax", "ymin", "ymax")]))
    left <- pmin(data$xmin, data$xmax)
    data$xmax <- pmax(data$xmin, data$xmax)
    data$xmin <- left
    n <- nrow(data)
    style <- paste(data$group, data$fill, data$alpha, data$colour,
      data$linewidth, data$linetype)
    joined <- c(FALSE, style[-1] == style[-n] &
      data$xmin[-1] >= data$xmin[-n] & data$xmin[-1] <= data$xmax[-n] + 1e-9)
    data$run <- cumsum(!joined)
    edges <- vctrs::vec_unique(data.frame(run = rep(data$run, 2),
      x = c(data$xmin, data$xmax)))
    edges <- vctrs::vec_slice(edges, order(edges$run, edges$x))
    from <- vctrs::vec_match(data[c("run", "xmin")], setNames(edges,
      c("run", "xmin")))
    to <- vctrs::vec_match(data[c("run", "xmax")], setNames(edges,
      c("run", "xmax")))
    spans <- to - from
    bar <- rep(seq_len(n), spans)
    piece <- data.frame(at = rep(from, spans) + sequence(spans) - 1,
      hi = pmax(data$ymin, data$ymax)[bar],
      lo = pmin(data$ymin, data$ymax)[bar])
    hi <- vctrs::vec_slice(piece, order(piece$at, -piece$hi))
    hi <- vctrs::vec_slice(hi, !duplicated(hi$at))
    lo <- vctrs::vec_slice(piece, order(piece$at, piece$lo))
    lo <- vctrs::vec_slice(lo, !duplicated(lo$at))
    m <- nrow(hi)
    corner <- vctrs::vec_slice(data,
      rep(match(edges$run[hi$at], data$run), each = 4))
    corner$x <- vctrs::vec_interleave(edges$x[hi$at], edges$x[hi$at + 1],
      edges$x[hi$at], edges$x[hi$at + 1])
    corner$y <- vctrs::vec_interleave(hi$hi, hi$hi, lo$lo, lo$lo)
    corner$group <- corner$run
    bottom <- rep(c(FALSE, FALSE, TRUE, TRUE), m)
    along <- seq_len(4 * m)
    ggplot2::GeomPolygon$draw_panel(
      vctrs::vec_slice(corner,
        order(corner$run, bottom, ifelse(bottom, -along, along))),
      panel_params, coord, lineend = lineend, linejoin = linejoin)
  })

# Each density cell reaches under the next cell of its row, which paints over
# it, so no antialiased seam shows between cells narrower than a device pixel
under_next <- function(d) {
  n <- nrow(d)
  abuts <- c(d$row[-1] == d$row[-n] & d$start[-1] >= d$start[-n] &
    d$start[-1] <= d$end[-n], FALSE)
  d$end[abuts] <- d$end[which(abuts) + 1]
  d
}

#' A quantitative track's marks as ggplot layers
#'
#' `geom_jbrowse_wiggle()` draws the marks as the display's plot type draws
#' them: bars from the origin, a stepped or interpolated line with its
#' whisker band, points, or density, each in the display's color. Map
#' `score`, `source` or any other column of the marks table to draw it your
#' own way. A line breaks into the `run`s the display breaks it into.
#' `geom_jbrowse_rules()` draws the result's `rules`, each labelled above
#' its line as the display labels it, in every row of a row layout.
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
  of <- function(mark, pick, interpolate = NULL) function(x) {
    d <- pick(x)
    drawn <- x$mark == mark &&
      (is.null(interpolate) || identical(interpolate, x$interpolate))
    if (drawn) d else d[0, , drop = FALSE]
  }
  marks <- function(x) in_rows(x, pixel_wide(x, !x$marks$band))
  cells <- function(x) under_next(marks(x))
  band <- function(x) in_rows(x, pixel_wide(x, x$marks$band))
  steps <- function(x) {
    m <- pixel_wide(x, !x$marks$band)
    d <- m[rep(seq_len(nrow(m)), each = 2), , drop = FALSE]
    d$x <- as.vector(rbind(m$start, m$end))
    in_rows(x, cut_at_origin(d, x$origin))
  }
  centers <- function(x) {
    d <- pixel_wide(x, !x$marks$band)
    d$x <- (d$start + d$end) / 2
    in_rows(x, cut_at_origin(d, x$origin))
  }
  list(
    jbrowse_layer(GeomJBrowseBars, data, band,
      for_geom(mapping, GeomJBrowseBars),
      ggplot2::aes(xmin = start, xmax = end, ymin = low, ymax = score,
        fill = I(color), group = row, alpha = I(opacity)), list(...)),
    jbrowse_layer(GeomJBrowseBars, data, of("bar", marks),
      for_geom(mapping, GeomJBrowseBars),
      ggplot2::aes(xmin = start, xmax = end, ymin = low, ymax = score,
        fill = I(color), group = row), list(...)),
    jbrowse_layer(ggplot2::GeomRect, data, of("rect", cells),
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = base, ymax = base + 1,
        fill = I(color)), list(...)),
    jbrowse_layer(ggplot2::GeomPoint, data, of("point", marks),
      for_geom(mapping, ggplot2::GeomPoint),
      ggplot2::aes(x = (start + end) / 2, y = score, colour = I(color)),
      list(size = 0.6, ...)),
    jbrowse_layer(ggplot2::GeomPath, data, of("line", steps, "step"),
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = score, colour = I(color), group = run),
      list(...)),
    jbrowse_layer(ggplot2::GeomPath, data, of("line", centers, "linear"),
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = score, colour = I(color), group = run),
      list(...)),
    jbrowse_layer(ggplot2::GeomRect, data, of("bar", clip_strips),
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = low, ymax = score,
        fill = I(color)), list(...)))
}

#' @rdname geom_jbrowse_wiggle
#' @export
geom_jbrowse_rules <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "rules", "geom_jbrowse_rules")
  rules <- function(x) {
    r <- if (x$row_layout) in_each_row(x, x$rules) else x$rules
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
  density <- object$mark == "rect"
  rows <- object$row_layout
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
