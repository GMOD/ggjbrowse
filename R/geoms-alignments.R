# The alignments track's marks, one geom per table. Each draws what the
# display draws by default and takes any column of its table in `mapping`.

# The pixels a base spans, which the display's sub-pixel fade reads, in a
# figure `width` pixels across its region, else 1,500
px_of <- function(x, width) {
  (width %||% x$width %||% DEFAULT_WIDTH) / diff(xlim_of(x))
}

# The share of a pileup row left empty above and below a read, from the
# engine's read pitch: a read `featureHeight` pixels tall over a row
# `featureSpacing` pixels taller
read_inset <- function(x) {
  p <- x$read_pitch
  if (is.null(p)) return(0)
  (1 - p$featureHeight / (p$featureHeight + p$featureSpacing)) / 2
}

fades <- function(x) !isFALSE(x$display$pileup$fadeLowFrequency)

with_inset <- function(d, model) {
  d$inset <- rep(read_inset(model), nrow(d))
  d
}

#' The alignments track's marks as ggplot layers
#'
#' One geom per table of a [jbrowse_alignments()] result. Each draws the
#' marks as the display does, in the display's colors, the table's `color`
#' column painted as it is (`I(color)`), so it trains no scale. Map any
#' column of the table and ggplot's own scales take it: `aes(fill = mapq)`
#' fills the reads by mapping quality, beside the other geoms' display
#' colors. [guides_jbrowse()] adds the display's keys. With `data = NULL` a
#' layer reads the result its plot was made from,
#' `ggplot(<jbrowse_alignments>)`.
#'
#' Each geom draws what its table says: `geom_jbrowse_junctions()` the
#' junctions whose `shown` is true, `geom_jbrowse_alleles()` the allele
#' segments whose `shown` is true, and the mismatches, deletions,
#' insertions and clips each faded by its `fade` while the display's
#' `pileup$fadeLowFrequency` is on. An indel also fades as it shrinks below a
#' pixel, and an insertion of 10 bp or more draws as the display's wider bar,
#' never faded by its `fade`. Edit the table to draw otherwise:
#' `a$junctions$shown <- a$junctions$count >= 5`, or
#' `a$mismatches$fade <- 1` for no fade.
#'
#' * `geom_jbrowse_reads()`: each read on its pileup row, split at its splice
#'   gaps
#' * `geom_jbrowse_qualities()`: under `baseColor` `"baseQuality"`, each
#'   aligned base in its quality's color; `aes(fill = score)` maps the Phred
#'   score through a scale of your own
#' * `geom_jbrowse_skips()`, `geom_jbrowse_deletions()`: the splice gaps'
#'   lines and the deletions' boxes, and `geom_jbrowse_deletion_lengths()`
#'   writes a deletion's length in it where the box is wide enough
#' * `geom_jbrowse_overlaps()`: under `unit = "chain"`, the neutral fill over
#'   the span where two reads of a chain both align
#' * `geom_jbrowse_chains()`: under `unit = "chain"`, the hairline across a
#'   row from the first to the last read of a chain, a read's mate or split
#'   segments sharing it
#' * `geom_jbrowse_mismatches()`, `geom_jbrowse_insertions()`,
#'   `geom_jbrowse_clips()`, `geom_jbrowse_modifications()`: the marks on the
#'   reads; `geom_jbrowse_letters()` writes the mismatched bases and the
#'   insertion lengths
#' * `geom_jbrowse_coverage()`, `geom_jbrowse_alleles()`,
#'   `geom_jbrowse_modification_coverage()`: the coverage band's depth and its
#'   stacks
#' * `geom_jbrowse_interbase()`, `geom_jbrowse_indicators()`: the band's
#'   histogram of insertions and clips, hanging from its top as the display
#'   hangs it (a stack as deep as the depth axis's top reaches halfway down),
#'   and the triangles over the positions where they pile up
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
#'   result, such as a filtered copy of it. `reads` takes no data frame:
#'   assign the filtered reads into the result, which drops their marks too
#' @param width the pixels the figure spans across its region, which decide
#'   the fade, whether a one-base mark reaches a pixel and how high the
#'   read-connection band reaches; `NULL` takes the result's `width`
#' @param ... passed to [ggplot2::layer()] params
#' @return a ggplot layer
#' @name geom_jbrowse_reads
#' @examples
#' bam <- system.file("extdata", "gapdh_direct_rna.bam",
#'   package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "chr12:6,534,201-6,538,700")
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
      fill = I(color)), list(...), rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_deletion_lengths <- function(mapping = NULL, data = NULL,
                                          width = NULL, ...) {
  data <- replacing(data, "gaps", "geom_jbrowse_deletion_lengths")
  jbrowse_layer(ggplot2::GeomText, data, function(x) {
    rule <- x$glyph$label
    d <- x$gaps[x$gaps$type == "deletion", , drop = FALSE]
    d$length <- d$end - d$start
    tall <- (x$read_pitch$featureHeight %||% 0) >= (rule$minHeight %||% Inf)
    d$alpha <- label_fade((d$end - d$start) * px_of(x, width),
      nchar(d$length) * (rule$digitPx %||% Inf), rule)
    if (!tall) d$alpha[] <- 0
    d[d$alpha > 0, , drop = FALSE]
  }, mapping, ggplot2::aes(x = (start + end) / 2, y = row + 0.5,
    label = length, alpha = I(alpha)),
    list(size = 1.7, colour = "white", fontface = "bold", ...))
}

# A size label's opacity from the pixels its feature spans and the pixels its
# digits need: none below `fadeLo` times the need, then `fadeFloor` ramping to
# one at `fadeHi` times it, as the display fades a deletion's length
label_fade <- function(span_px, need_px, rule) {
  t <- pmin(1, pmax(0, (span_px - need_px) / (need_px * (rule$fadeHi %||% 2) -
    need_px)))
  ramp <- t * t * (3 - 2 * t)
  ifelse(span_px < need_px * (rule$fadeLo %||% Inf), 0,
    rule$fadeFloor + (1 - rule$fadeFloor) * ramp)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_overlaps <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "overlaps", "geom_jbrowse_overlaps")
  jbrowse_layer(ggplot2::GeomRect, data,
    function(x) summary_table(x, "overlaps"), mapping,
    ggplot2::aes(xmin = start, xmax = end, ymin = row + inset,
      ymax = row + 1 - inset, fill = I(color)), list(...),
    rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_chains <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "chains", "geom_jbrowse_chains")
  jbrowse_layer(ggplot2::GeomSegment, data,
    function(x) summary_table(x, "chains"), mapping,
    ggplot2::aes(x = start, xend = end, y = row + 0.5, yend = row + 0.5,
      colour = I(color), alpha = I(opacity)), list(linewidth = px_to_mm(1), ...))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_skips <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "gaps", "geom_jbrowse_skips")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    d <- x$gaps[x$gaps$type == "skip", , drop = FALSE]
    d$alpha <- rep(x$glyph$intronAlpha %||% 1, nrow(d))
    d
  }, mapping, ggplot2::aes(x = start, xend = end, y = row + 0.5,
    yend = row + 0.5, colour = I(color), alpha = I(alpha)),
    list(linewidth = px_to_mm(1), ...))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_deletions <- function(mapping = NULL, data = NULL,
                                   width = NULL, ...) {
  data <- replacing(data, "gaps", "geom_jbrowse_deletions")
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- x$gaps[x$gaps$type == "deletion", , drop = FALSE]
    px <- px_of(x, width)
    span <- (d$end - d$start) * px
    d$alpha <- mark_alpha(span^2, d$fade, fades(x)) * size_alpha(span)
    widen <- pmax(0, 1 / px - (d$end - d$start)) / 2
    d$xmin <- d$start - widen
    d$xmax <- d$end + widen
    d
  }, mapping, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = row + inset,
    ymax = row + 1 - inset, fill = I(color), alpha = I(alpha)), list(...), rows = TRUE)
}

# A per-position table drawn as rects one base wide, a pixel at least, over
# the read's row unless `defaults` says otherwise; `prep` adjusts the table
# first, and where `faded` the rects fade by the pixels per base
position_rects <- function(table, geom, data, mapping, width, params,
                           defaults, prep = NULL, faded = FALSE,
                           rows = TRUE) {
  data <- replacing(data, table, geom)
  defaults <- with_defaults(defaults, ggplot2::aes(xmin = position,
    xmax = xmax, ymin = row + inset, ymax = row + 1 - inset))
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- if (table %in% SUMMARIES) summary_table(x, table) else x[[table]]
    if (!is.null(prep)) d <- prep(d, x)
    d$xmax <- d$position + max(1, 1 / px_of(x, width))
    if (faded) d$alpha <- mark_alpha(px_of(x, width), d$fade, fades(x))
    d
  }, mapping, defaults, params, rows = rows)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_modifications <- function(mapping = NULL, data = NULL,
                                       width = NULL, ...) {
  position_rects("modifications", "geom_jbrowse_modifications", data, mapping,
    width, list(...), ggplot2::aes(fill = I(color), alpha = I(opacity)))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_qualities <- function(mapping = NULL, data = NULL,
                                   width = NULL, ...) {
  position_rects("qualities", "geom_jbrowse_qualities", data, mapping, width,
    list(...), ggplot2::aes(fill = I(color)))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_mismatches <- function(mapping = NULL, data = NULL,
                                    width = NULL, ...) {
  position_rects("mismatches", "geom_jbrowse_mismatches", data, mapping,
    width, list(...), ggplot2::aes(fill = I(color), alpha = I(alpha)),
    faded = TRUE)
}

# An indel's share of a pixel, which fades it on top of the frequency fade,
# none below one 8-bit alpha step: the shaders' sizeAlpha
size_alpha <- function(span_px) {
  a <- pmin(1, pmax(0, span_px))
  ifelse(a <= 1 / 255, 0, a)
}

LONG_INSERTION_BP <- 10
BOXED_INSERTION_PX <- 15

# An insertion marker's width in px, as insertionGlyph.slang sizes it: a 1 px
# tick, a long one's bar a third of its span (1 to 5 px), or a box for its
# length once that span reaches BOXED_INSERTION_PX in a row tall enough for
# text
insertion_px <- function(length, px, height, text_height) {
  long <- length >= LONG_INSERTION_BP
  boxed <- long & length * px >= BOXED_INSERTION_PX & height >= text_height
  digits <- pmin(5, floor(log10(pmax(1, length))) + 1)
  ifelse(boxed, 6 * digits + 10,
    ifelse(long, pmin(5, pmax(1, length * px / 3)), 1))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_insertions <- function(mapping = NULL, data = NULL,
                                    width = NULL, ...) {
  data <- replacing(data, "insertions", "geom_jbrowse_insertions")
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- x$insertions
    px <- px_of(x, width)
    long <- d$length >= LONG_INSERTION_BP
    d$alpha <- ifelse(long, 1, mark_alpha(px^2, d$fade, fades(x))) *
      size_alpha(d$length * px)
    half <- insertion_px(d$length, px, x$read_pitch$featureHeight %||% 0,
      x$glyph$label$minHeight %||% Inf) / px / 2
    d$xmin <- d$position - half
    d$xmax <- d$position + half
    d
  }, mapping, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = row + inset,
    ymax = row + 1 - inset, fill = I(color), alpha = I(alpha)), list(...),
    rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_clips <- function(mapping = NULL, data = NULL,
                               width = NULL, ...) {
  data <- replacing(data, "clips", "geom_jbrowse_clips")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    d <- x$clips
    d$alpha <- mark_alpha(px_of(x, width), d$fade, fades(x))
    d
  }, mapping, ggplot2::aes(x = position, xend = position, y = row + inset,
    yend = row + 1 - inset, colour = I(color), alpha = I(alpha)),
    list(linewidth = px_to_mm(1), ...), rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_letters <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, NULL, "geom_jbrowse_letters")
  list(
    jbrowse_layer(ggplot2::GeomText, data, function(x) x$mismatches,
      for_geom(mapping, ggplot2::GeomText),
      ggplot2::aes(x = position + 0.5, y = row + 0.5, label = base),
      list(size = 2.6, colour = "black", ...)),
    jbrowse_layer(ggplot2::GeomLabel, data, function(x) x$insertions,
      for_geom(mapping, ggplot2::GeomLabel),
      ggplot2::aes(x = position, y = row + 0.5, label = length, fill = I(color)),
      list(size = 2.2, colour = "white",
        label.padding = grid::unit(0.1, "lines"), linewidth = 0, ...)))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_coverage <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "coverage", "geom_jbrowse_coverage")
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- summary_table(x, "coverage")
    d$base <- rep(bar_base(x), nrow(d))
    d[d$depth > d$base, , drop = FALSE]
  }, mapping, ggplot2::aes(xmin = start, xmax = end, ymin = base, ymax = depth,
      fill = I(color)), list(...))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_alleles <- function(mapping = NULL, data = NULL, width = NULL,
                                 ...) {
  position_rects("allele_coverage", "geom_jbrowse_alleles", data, mapping,
    width, list(...), with_defaults(ggplot2::aes(fill = I(color)), stack_rect),
    rows = FALSE, prep = function(a, x) on_bar(x, a[shown(a), , drop = FALSE]))
}

stack_rect <- ggplot2::aes(ymin = low, ymax = high)

# The depth the band's bars grow from: the axis's bottom, which a log axis
# floors at one read
bar_base <- function(x) x$y_scale$domain[1] %||% 0

# A stack's `from` and `to` as the display places them, shares of its bar's
# drawn height from the bar's base, so a share means the same length on a
# log or symlog axis as on a linear one
on_bar <- function(x, d) {
  tr <- score_transform(x$y_scale)
  base <- tr$transform(bar_base(x))
  at <- function(share) tr$inverse(base + share * (tr$transform(d$depth) -
    base))
  d$low <- at(d$from)
  d$high <- at(d$to)
  d
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_interbase <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "interbase_coverage", "geom_jbrowse_interbase")
  jbrowse_layer(GeomInterbase, data, function(x) {
    d <- summary_table(x, "interbase_coverage")
    full <- 2 * x$y_scale$domain[2]
    n <- nrow(d)
    d$hang <- if (full > 0) d$from / full else rep(0, n)
    d$hang_end <- if (full > 0) d$to / full else rep(0, n)
    d$top <- rep(px_to_pt(x$glyph$interbaseTop), n)
    d$inset <- rep(px_to_pt(x$glyph$coverageInset), n)
    d$width <- rep(px_to_pt(1), n)
    d
  }, mapping, ggplot2::aes(x = position, hang = hang, hang_end = hang_end,
    top = I(top), inset = I(inset), width = I(width), fill = I(color)),
    list(...))
}

# Each position's stack as the display hangs it from the band's top, below
# the indicator triangles' tip row: `hang` and `hang_end` are fractions of
# the depth area, half of it for a stack as deep as the domain's top, so the
# hang is linear on any y scale and trains none
GeomInterbase <- ggplot2::ggproto("GeomInterbase", ggplot2::Geom,
  required_aes = c("x", "hang", "hang_end", "top", "inset", "width"),
  default_aes = ggplot2::aes(fill = "grey30", alpha = NA),
  draw_key = ggplot2::draw_key_rect,
  draw_panel = function(data, panel_params, coord) {
    d <- coord$transform(data, panel_params)
    if (!nrow(d)) return(grid::nullGrob())
    area <- grid::unit(1, "npc") - grid::unit(2 * d$inset, "pt")
    below <- function(f) {
      grid::unit(1, "npc") - grid::unit(d$top, "pt") - f * area
    }
    grid::rectGrob(
      x = grid::unit(d$x, "npc"), y = below(d$hang),
      width = grid::unit(d$width, "pt"),
      height = below(d$hang) - below(d$hang_end),
      just = c("left", "top"),
      gp = grid::gpar(fill = ggplot2::alpha(d$fill, d$alpha), col = NA))
  })

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_indicators <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "indicators", "geom_jbrowse_indicators")
  jbrowse_layer(GeomIndicator, data, function(x) {
    d <- summary_table(x, "indicators")
    d$width <- rep(px_to_pt(x$glyph$indicatorWidth), nrow(d))
    d$height <- rep(px_to_pt(x$glyph$indicatorHeight), nrow(d))
    d
  }, mapping, ggplot2::aes(x = position, fill = I(color), width = I(width),
    height = I(height)), list(...))
}

# A display's px in grid's points, which are 72.27 to the inch
px_to_pt <- function(px) px / PX_PER_INCH * 72.27

# Pixels at the figure's 150 an inch as ggplot's linewidth, in mm
px_to_mm <- function(px) px / PX_PER_INCH * 25.4

# A triangle `width` by `height` points at each `x`, its apex pointing down
# from the panel's top edge, whatever the y scale
GeomIndicator <- ggplot2::ggproto("GeomIndicator", ggplot2::Geom,
  required_aes = c("x", "width", "height"),
  default_aes = ggplot2::aes(fill = "grey30", alpha = NA),
  draw_key = ggplot2::draw_key_polygon,
  draw_panel = function(data, panel_params, coord) {
    d <- coord$transform(data, panel_params)
    n <- nrow(d)
    if (!n) return(grid::nullGrob())
    hw <- rep(d$width / 2, each = 3)
    grid::polygonGrob(
      x = grid::unit(rep(d$x, each = 3), "npc") +
        grid::unit(hw * rep(c(-1, 1, 0), n), "pt"),
      y = grid::unit(1, "npc") -
        grid::unit(rep(d$height, each = 3) * rep(c(0, 0, 1), n), "pt"),
      id = rep(seq_len(n), each = 3),
      gp = grid::gpar(fill = ggplot2::alpha(d$fill, d$alpha), col = NA))
  })

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_modification_coverage <- function(mapping = NULL, data = NULL,
                                               width = NULL, ...) {
  position_rects("modification_coverage", "geom_jbrowse_modification_coverage",
    data, mapping, width, list(...),
    with_defaults(ggplot2::aes(fill = I(color), alpha = I(opacity)), stack_rect),
    rows = FALSE, prep = on_bar_of)
}

on_bar_of <- function(d, x) on_bar(x, d)

# Each junction the band draws as its arc's points, its stroke the engine's
# px as a linewidth at the figure's PX_PER_INCH
junction_paths <- function(x) {
  j <- shown_junctions(x)
  arcs <- junction_arcs(j)
  arcs$color <- j$color[arcs$id]
  arcs$section <- j$section[arcs$id]
  arcs$count <- j$count[arcs$id]
  arcs$strand <- j$strand[arcs$id]
  arcs$width <- px_to_mm(j$strokeWidth[arcs$id])
  arcs
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_junctions <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "junctions", "geom_jbrowse_junctions")
  jbrowse_layer(ggplot2::GeomPath, data,
    junction_paths, mapping,
    ggplot2::aes(x = x, y = y, group = id, colour = I(color),
      linewidth = I(width)), list(lineend = "round", ...))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_junction_counts <- function(mapping = NULL, data = NULL,
                                         ...) {
  data <- replacing(data, "junctions", "geom_jbrowse_junction_counts")
  jbrowse_layer(ggplot2::GeomLabel, data, function(x) {
    arcs <- junction_paths(x)
    xlim <- xlim_of(x)
    shown <- arcs[arcs$x >= xlim[1] & arcs$x <= xlim[2], , drop = FALSE]
    top <- shown[order(shown$id, -shown$y), , drop = FALSE]
    top[!duplicated(top$id), , drop = FALSE]
  }, mapping, ggplot2::aes(x = x, y = y, label = count),
    list(size = 2.6, linewidth = 0,
      label.padding = grid::unit(0.15, "lines"), ...))
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

is_cloud <- function(x) {
  band <- x$display$readConnections
  isTRUE(band$show) && identical(band$mark, "cloud")
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_arcs <- function(mapping = NULL, data = NULL, width = NULL,
                              ...) {
  data <- replacing(data, "arcs", "geom_jbrowse_arcs")
  domes <- function(x) {
    a <- summary_table(x, "arcs")
    connection_domes(a[a$shape == "arc", , drop = FALSE], arc_top(x, width))
  }
  bars <- function(x) {
    b <- summary_table(x, "arcs")
    b <- b[b$shape != "arc", , drop = FALSE]
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
      ggplot2::aes(x = x, y = y, group = id, colour = I(color),
        linewidth = I(px_to_mm(width))), list(...)),
    jbrowse_layer(ggplot2::GeomSegment, data, bars,
      for_geom(mapping, ggplot2::GeomSegment),
      ggplot2::aes(x = start, xend = end, y = y, yend = y, colour = I(color),
        alpha = I(opacity), linewidth = I(px_to_mm(width)),
        linetype = I(ifelse(shape == "flatSplit", "22", "solid"))),
      list(...)),
    jbrowse_layer(ggplot2::GeomPoint, data, feet,
      for_geom(mapping, ggplot2::GeomPoint),
      ggplot2::aes(x = x, y = y, colour = I(color)),
      list(shape = 15, size = 1.3, ...)))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_arc_ticks <- function(mapping = NULL, data = NULL,
                                   width = NULL, ...) {
  data <- replacing(data, "arc_ticks", "geom_jbrowse_arc_ticks")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    t <- summary_table(x, "arc_ticks")
    t$bottom <- rep(if (is_cloud(x)) 1 else 0, nrow(t))
    t$top <- rep(arc_top(x, width), nrow(t))
    t
  }, mapping, ggplot2::aes(x = position, xend = position, y = bottom,
    yend = top, colour = I(color), linewidth = I(px_to_mm(width))),
    list(...))
}

# The arc band's reach: the read cloud's widest insert, or as many bp as the
# band's pixels span in a figure `width` pixels across, so an arc rises as
# high as it is wide and one wider than the band runs off it, as the display
# draws it
arc_top <- function(x, width = NULL) {
  if (is_cloud(x)) max(2, x$arc_span_max)
  else x$display$readConnections$height / px_of(x, width)
}
