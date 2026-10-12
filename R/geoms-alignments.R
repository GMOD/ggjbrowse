# The alignments track's marks, one geom per table. Each draws what the
# display draws by default and takes any column of its table in `mapping`.

# The pixels a base spans, which the display's sub-pixel fade reads: what a
# geom is told, else a figure 1,500 pixels wide
px_of <- function(x, px_per_bp) {
  px_per_bp %||% (10 * PX_PER_INCH / (x$region$end - x$region$start))
}

xlim_of <- function(x) c(x$region$start, x$region$end)

# The share of a pileup row left empty above and below a read: the display
# draws a read `featureHeight` pixels tall and, above 3, a 1 px gap under it
# (featureSpacingForHeight), so a compact read fills its row. A table drawn
# without its result takes the display's default of 7.
read_inset <- function(x) {
  h <- x$settings$featureHeight %||% 7
  (1 - h / (h + (h > 3))) / 2
}

with_inset <- function(d, model) {
  d$inset <- rep(read_inset(model), nrow(d))
  d
}

#' The alignments track's marks as ggplot layers
#'
#' One geom per table of a [jbrowse_alignments()] result. Each draws the
#' marks as the display does, in the display's colors (the table's `color`
#' column, an identity scale), and takes any column of its table in
#' `mapping`: `aes(fill = mapq)` fills the reads by mapping quality through
#' ggplot's own scales. With `data = NULL` a layer reads the result its plot
#' was made from, `ggplot(<jbrowse_alignments>)`.
#'
#' * `geom_jbrowse_reads()`: each read on its pileup row, split at its splice
#'   gaps
#' * `geom_jbrowse_qualities()`: under `baseColor` `"baseQuality"`, each
#'   aligned base in its quality's color; `aes(fill = score)` maps the Phred
#'   score through a scale of your own
#' * `geom_jbrowse_skips()`, `geom_jbrowse_deletions()`: the splice gaps'
#'   lines and the deletions' boxes
#' * `geom_jbrowse_mismatches()`, `geom_jbrowse_insertions()`,
#'   `geom_jbrowse_clips()`, `geom_jbrowse_modifications()`: the marks on the
#'   reads; `geom_jbrowse_letters()` writes the mismatched bases and the
#'   insertion lengths
#' * `geom_jbrowse_coverage()`, `geom_jbrowse_alleles()`,
#'   `geom_jbrowse_modification_coverage()`: the coverage band's depth and its
#'   stacks
#' * `geom_jbrowse_junctions()`, `geom_jbrowse_junction_counts()`: the
#'   sashimi arcs and their read counts
#' * `geom_jbrowse_arcs()`, `geom_jbrowse_arc_ticks()`: the read-connection
#'   band's arcs (or read-cloud bars, with a square at each foot) and its
#'   interchromosomal ticks
#'
#' The y aesthetics are the display's: pileup rows down from 0
#' (`scale_y_reverse()` draws them as the display stacks them), read depth,
#' and an arc's genomic radius or a read-cloud bar's insert size.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_alignments()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param fade fade a mark narrower than a pixel by its site's share of the
#'   depth, as the display does (`fadeLowFrequencyMismatches`)
#' @param px_per_bp how many pixels a base spans, which the fade and a
#'   mismatch's line width read; `NULL` takes a figure 1,500 pixels wide
#' @param min_allele_frequency leave out an allele segment whose share of the
#'   depth is below this, the display's `coverageSnpMinFrequency`
#' @param min_reads leave out a junction with fewer supporting reads, the
#'   display's `minSashimiScore`
#' @param ... passed to [ggplot2::layer()] params
#' @return a ggplot layer
#' @name geom_jbrowse_reads
#' @examples
#' bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "ctgA", 1000, 3000)
#' library(ggplot2)
#' ggplot(a) +
#'   geom_jbrowse_reads(aes(fill = mapq)) +
#'   geom_jbrowse_mismatches() +
#'   scale_fill_viridis_c() +
#'   scale_y_reverse()
NULL

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_reads <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "reads", "geom_jbrowse_reads")
  jbrowse_layer(ggplot2::GeomRect, data, read_blocks, mapping,
    ggplot2::aes(xmin = start, xmax = end, ymin = row + inset,
      ymax = row + 1 - inset,
      fill = color), list(...), rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_skips <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "gaps", "geom_jbrowse_skips")
  jbrowse_layer(ggplot2::GeomSegment, data,
    function(x) x$gaps[x$gaps$type == "skip", , drop = FALSE], mapping,
    ggplot2::aes(x = start, xend = end, y = row + 0.5, yend = row + 0.5,
      colour = color), list(linewidth = 0.3, ...), keyed = FALSE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_deletions <- function(mapping = NULL, data = NULL, fade = TRUE,
                                   px_per_bp = NULL, ...) {
  data <- replacing(data, "gaps", "geom_jbrowse_deletions")
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- x$gaps[x$gaps$type == "deletion", , drop = FALSE]
    d$alpha <- mark_alpha(((d$end - d$start) * px_of(x, px_per_bp))^2,
      d$fade, fade)
    d
  }, mapping, ggplot2::aes(xmin = start, xmax = end, ymin = row + 0.35,
    ymax = row + 0.65, fill = color, alpha = I(alpha)), list(...),
    keyed = FALSE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_modifications <- function(mapping = NULL, data = NULL,
                                       px_per_bp = NULL, ...) {
  data <- replacing(data, "modifications", "geom_jbrowse_modifications")
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- x$modifications
    d$xmax <- d$position + max(1, 1 / px_of(x, px_per_bp))
    d
  }, mapping,
    ggplot2::aes(xmin = position, xmax = xmax, ymin = row + inset,
      ymax = row + 1 - inset, fill = color, alpha = I(opacity)), list(...),
    rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_qualities <- function(mapping = NULL, data = NULL,
                                   px_per_bp = NULL, ...) {
  data <- replacing(data, "qualities", "geom_jbrowse_qualities")
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- x$qualities
    d$xmax <- d$position + max(1, 1 / px_of(x, px_per_bp))
    d
  }, mapping,
    ggplot2::aes(xmin = position, xmax = xmax, ymin = row + inset,
      ymax = row + 1 - inset, fill = color), list(...), keyed = FALSE,
    rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_mismatches <- function(mapping = NULL, data = NULL, fade = TRUE,
                                    px_per_bp = NULL, ...) {
  data <- replacing(data, "mismatches", "geom_jbrowse_mismatches")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    d <- x$mismatches
    d$alpha <- mark_alpha(px_of(x, px_per_bp), d$fade, fade)
    d$width <- rep(cell_width(xlim_of(x)), nrow(d))
    d
  }, mapping, ggplot2::aes(x = position + 0.5, xend = position + 0.5,
    y = row + inset, yend = row + 1 - inset, colour = color, alpha = I(alpha),
    linewidth = I(width)), list(lineend = "butt", ...), rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_insertions <- function(mapping = NULL, data = NULL, fade = TRUE,
                                    px_per_bp = NULL, ...) {
  data <- replacing(data, "insertions", "geom_jbrowse_insertions")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    d <- x$insertions
    d$alpha <- mark_alpha(px_of(x, px_per_bp)^2, d$fade, fade)
    d
  }, mapping, ggplot2::aes(x = position, xend = position,
    y = row + inset / 2, yend = row + 1 - inset / 2, colour = color,
    alpha = I(alpha)),
    list(linewidth = 0.4, ...), keyed = FALSE, rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_clips <- function(mapping = NULL, data = NULL, fade = TRUE,
                               px_per_bp = NULL, ...) {
  data <- replacing(data, "clips", "geom_jbrowse_clips")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    d <- x$clips
    d$alpha <- mark_alpha(px_of(x, px_per_bp), d$fade, fade)
    d
  }, mapping, ggplot2::aes(x = position, xend = position,
    y = row + inset / 2, yend = row + 1 - inset / 2, colour = color,
    alpha = I(alpha)),
    list(linewidth = 0.7, ...), keyed = FALSE, rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_letters <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, NULL, "geom_jbrowse_letters")
  list(
    jbrowse_layer(ggplot2::GeomText, data, function(x) x$mismatches,
      for_geom(mapping, ggplot2::GeomText),
      ggplot2::aes(x = position + 0.5, y = row + 0.5, label = base),
      list(size = 2.6, colour = "black", ...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomLabel, data, function(x) x$insertions,
      for_geom(mapping, ggplot2::GeomLabel),
      ggplot2::aes(x = position, y = row + 0.5, label = length, fill = color),
      list(size = 2.2, colour = "white",
        label.padding = grid::unit(0.1, "lines"), linewidth = 0, ...),
      keyed = FALSE))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_coverage <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "coverage", "geom_jbrowse_coverage")
  jbrowse_layer(ggplot2::GeomRect, data, function(x) x$coverage, mapping,
    ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = depth,
      fill = color), list(...), keyed = FALSE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_alleles <- function(mapping = NULL, data = NULL,
                                 min_allele_frequency = 0, ...) {
  data <- replacing(data, "allele_coverage", "geom_jbrowse_alleles")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    a <- x$allele_coverage
    a <- a[a$to - a$from >= min_allele_frequency, , drop = FALSE]
    a$width <- rep(cell_width(xlim_of(x)), nrow(a))
    a
  }, mapping, ggplot2::aes(x = position + 0.5, xend = position + 0.5,
    y = depth * from, yend = depth * to, colour = color,
    linewidth = I(width)), list(lineend = "butt", ...))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_modification_coverage <- function(mapping = NULL, data = NULL,
                                               ...) {
  data <- replacing(data, "modification_coverage", "geom_jbrowse_modification_coverage")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    m <- x$modification_coverage
    m$width <- rep(cell_width(xlim_of(x)), nrow(m))
    m
  }, mapping, ggplot2::aes(x = position + 0.5, xend = position + 0.5,
    y = depth * from, yend = depth * to, colour = color, alpha = I(opacity),
    linewidth = I(width)), list(lineend = "butt", ...), keyed = FALSE)
}

junction_paths <- function(x, min_reads) {
  j <- x$junctions[x$junctions$count >= min_reads, , drop = FALSE]
  arcs <- junction_arcs(j)
  arcs$color <- j$color[arcs$id]
  arcs$section <- j$section[arcs$id]
  arcs$count <- j$count[arcs$id]
  arcs$strand <- j$strand[arcs$id]
  arcs$width <- pmin(4, pmax(0.4, log(arcs$count + 1) / 2))
  arcs
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_junctions <- function(mapping = NULL, data = NULL,
                                   min_reads = 2, ...) {
  data <- replacing(data, "junctions", "geom_jbrowse_junctions")
  jbrowse_layer(ggplot2::GeomPath, data,
    function(x) junction_paths(x, min_reads), mapping,
    ggplot2::aes(x = x, y = y, group = id, colour = color,
      linewidth = I(width)), list(lineend = "round", ...))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_junction_counts <- function(mapping = NULL, data = NULL,
                                         min_reads = 2, ...) {
  data <- replacing(data, "junctions", "geom_jbrowse_junction_counts")
  jbrowse_layer(ggplot2::GeomLabel, data, function(x) {
    arcs <- junction_paths(x, min_reads)
    xlim <- xlim_of(x)
    shown <- arcs[arcs$x >= xlim[1] & arcs$x <= xlim[2], , drop = FALSE]
    top <- shown[order(shown$id, -shown$y), , drop = FALSE]
    top[!duplicated(top$id), , drop = FALSE]
  }, mapping, ggplot2::aes(x = x, y = y, label = count),
    list(size = 2.6, linewidth = 0,
      label.padding = grid::unit(0.15, "lines"), ...), keyed = FALSE)
}

# Each connection's arc as points: a dome from foot to foot, its height the
# genomic radius up to the band's reach `top`, as the display rises an arc as
# high as it is wide and flattens one wider than its band
connection_domes <- function(arcs, top, points = 48) {
  t <- rep(seq(0, pi, length.out = points), nrow(arcs))
  id <- rep(seq_len(nrow(arcs)), each = points)
  d <- arcs[id, , drop = FALSE]
  d$id <- id
  d$x <- d$start + (d$end - d$start) * (1 - cos(t)) / 2
  d$height <- pmin(d$y, top)
  d$y <- d$height * sin(t)
  d
}

is_cloud <- function(x) x$read_connections == "cloud"

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_arcs <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "arcs", "geom_jbrowse_arcs")
  domes <- function(x) connection_domes(x$arcs[x$arcs$shape == "arc", ,
    drop = FALSE], arc_top(x))
  bars <- function(x) {
    b <- x$arcs[x$arcs$shape != "arc", , drop = FALSE]
    b$y <- pmax(1, b$y)
    b
  }
  feet <- function(x) {
    b <- bars(x)
    rbind(transform(b, x = start), transform(b, x = end))
  }
  list(
    jbrowse_layer(ggplot2::GeomPath, data, domes,
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = y, group = id, colour = color,
        linewidth = I(width / ggplot2::.pt)), list(...)),
    jbrowse_layer(ggplot2::GeomSegment, data, bars,
      for_geom(mapping, ggplot2::GeomSegment),
      ggplot2::aes(x = start, xend = end, y = y, yend = y, colour = color,
        alpha = I(opacity), linewidth = I(width / ggplot2::.pt),
        linetype = I(ifelse(shape == "flatSplit", "22", "solid"))),
      list(...)),
    jbrowse_layer(ggplot2::GeomPoint, data, feet,
      for_geom(mapping, ggplot2::GeomPoint),
      ggplot2::aes(x = x, y = y, colour = color),
      list(shape = 15, size = 1.3, ...)))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_arc_ticks <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "arc_ticks", "geom_jbrowse_arc_ticks")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    t <- x$arc_ticks
    t$bottom <- rep(if (is_cloud(x)) 1 else 0, nrow(t))
    t$top <- rep(arc_top(x), nrow(t))
    t
  }, mapping, ggplot2::aes(x = position, xend = position, y = bottom,
    yend = top, colour = color, linewidth = I(width / ggplot2::.pt)),
    list(...))
}

# The arc band's reach: the read cloud's widest insert, or as many bp as the
# band's pixels span at the figure's scale, so an arc rises as high as it is
# wide and one wider than the band runs off it, as the display draws it
arc_top <- function(x) {
  if (is_cloud(x)) max(2, x$arc_span_max)
  else x$settings$readConnectionsHeight / px_of(x, NULL)
}
