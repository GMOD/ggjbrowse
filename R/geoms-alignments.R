# The alignments track's marks, one geom per table. Each draws what the
# display draws by default and takes any column of its table in `mapping`.

# The pixels a base spans, which the display's sub-pixel fade reads, in a
# figure `width` pixels across its region, else 1,500
px_of <- function(x, width) {
  (width %||% (10 * PX_PER_INCH)) / (x$region$end - x$region$start)
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
#'   result, such as a filtered copy of it
#' @param fade fade a mark narrower than a pixel by its site's share of the
#'   depth, as the display does (`fadeLowFrequencyMismatches`)
#' @param width the pixels the figure spans across its region, which decide
#'   the fade, whether a one-base mark reaches a pixel and how high the
#'   read-connection band reaches; `NULL` takes 1,500
#' @param min_allele_frequency leave out an allele segment whose share of the
#'   depth is below this, the display's `coverageSnpMinFrequency`
#' @param min_reads leave out a junction with fewer supporting reads, the
#'   display's `minSashimiScore`
#' @param ... passed to [ggplot2::layer()] params
#' @return a ggplot layer
#' @name geom_jbrowse_reads
#' @examples
#' bam <- system.file("extdata", "gapdh_direct_rna.bam",
#'   package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "chr12", 6534200, 6538700)
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
      fill = color), list(...), keyed = "complete", rows = TRUE)
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
    tall <- (x$settings$featureHeight %||% 7) >= (rule$minHeight %||% Inf)
    d$alpha <- label_fade((d$end - d$start) * px_of(x, width),
      nchar(d$length) * (rule$digitPx %||% Inf), rule)
    if (!tall) d$alpha[] <- 0
    d[d$alpha > 0, , drop = FALSE]
  }, mapping, ggplot2::aes(x = (start + end) / 2, y = row + 0.5,
    label = length, alpha = I(alpha)),
    list(size = 1.7, colour = "white", fontface = "bold", ...),
    keyed = FALSE)
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
  jbrowse_layer(ggplot2::GeomRect, data, function(x) x$overlaps, mapping,
    ggplot2::aes(xmin = start, xmax = end, ymin = row + inset,
      ymax = row + 1 - inset, fill = color), list(...), keyed = FALSE,
    rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_chains <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "chains", "geom_jbrowse_chains")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$chains, mapping,
    ggplot2::aes(x = start, xend = end, y = row + 0.5, yend = row + 0.5,
      colour = color, alpha = I(opacity)), list(linewidth = 0.3, ...),
    keyed = FALSE)
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
    yend = row + 0.5, colour = color, alpha = I(alpha)),
    list(linewidth = 0.3, ...), keyed = FALSE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_deletions <- function(mapping = NULL, data = NULL, fade = TRUE,
                                   width = NULL, ...) {
  data <- replacing(data, "gaps", "geom_jbrowse_deletions")
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- x$gaps[x$gaps$type == "deletion", , drop = FALSE]
    d$alpha <- mark_alpha(((d$end - d$start) * px_of(x, width))^2,
      d$fade, fade)
    d
  }, mapping, ggplot2::aes(xmin = start, xmax = end, ymin = row + inset,
    ymax = row + 1 - inset, fill = color, alpha = I(alpha)), list(...),
    keyed = FALSE, rows = TRUE)
}

# A per-position table drawn as rects one base wide, a pixel at least, over
# the read's row unless `defaults` says otherwise, painted by `rect`; `prep`
# adjusts the table first and `fade` fades the rects by the pixels per base
position_rects <- function(table, geom, data, mapping, width, params,
                           defaults, rect = ggplot2::GeomRect, prep = NULL,
                           fade = NULL, ..., rows = TRUE) {
  data <- replacing(data, table, geom)
  defaults <- with_defaults(defaults, ggplot2::aes(xmin = position,
    xmax = xmax, ymin = row + inset, ymax = row + 1 - inset))
  jbrowse_layer(rect, data, function(x) {
    d <- x[[table]]
    if (!is.null(prep)) d <- prep(d, x)
    d$xmax <- d$position + max(1, 1 / px_of(x, width))
    if (!is.null(fade)) d$alpha <- mark_alpha(px_of(x, width), d$fade, fade)
    d
  }, mapping, defaults, params, ..., rows = rows)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_modifications <- function(mapping = NULL, data = NULL,
                                       width = NULL, ...) {
  position_rects("modifications", "geom_jbrowse_modifications", data, mapping,
    width, list(...), ggplot2::aes(fill = color, alpha = I(opacity)))
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_qualities <- function(mapping = NULL, data = NULL,
                                   width = NULL, ...) {
  position_rects("qualities", "geom_jbrowse_qualities", data, mapping, width,
    list(...), ggplot2::aes(fill = color), keyed = FALSE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_mismatches <- function(mapping = NULL, data = NULL, fade = TRUE,
                                    width = NULL, ...) {
  position_rects("mismatches", "geom_jbrowse_mismatches", data, mapping,
    width, list(...), ggplot2::aes(colour = color, alpha = I(alpha)),
    rect = GeomSwatch, fade = fade)
}

# A per-position table drawn as vertical ticks, faded as the pixels per base
# raised to `power` shrink
position_ticks <- function(table, geom, data, mapping, fade, width, params,
                           power) {
  data <- replacing(data, table, geom)
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    d <- x[[table]]
    d$alpha <- mark_alpha(px_of(x, width)^power, d$fade, fade)
    d
  }, mapping, ggplot2::aes(x = position, xend = position,
    y = row + inset / 2, yend = row + 1 - inset / 2, colour = color,
    alpha = I(alpha)), params, keyed = FALSE, rows = TRUE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_insertions <- function(mapping = NULL, data = NULL, fade = TRUE,
                                    width = NULL, ...) {
  position_ticks("insertions", "geom_jbrowse_insertions", data, mapping, fade,
    width, list(linewidth = 0.4, ...), power = 2)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_clips <- function(mapping = NULL, data = NULL, fade = TRUE,
                               width = NULL, ...) {
  position_ticks("clips", "geom_jbrowse_clips", data, mapping, fade, width,
    list(linewidth = 0.7, ...), power = 1)
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
  jbrowse_layer(ggplot2::GeomRect, data, function(x) {
    d <- x$coverage
    d$base <- rep(bar_base(x), nrow(d))
    d[d$depth > d$base, , drop = FALSE]
  }, mapping, ggplot2::aes(xmin = start, xmax = end, ymin = base, ymax = depth,
      fill = color), list(...), keyed = FALSE)
}

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_alleles <- function(mapping = NULL, data = NULL,
                                 min_allele_frequency = 0, width = NULL, ...) {
  position_rects("allele_coverage", "geom_jbrowse_alleles", data, mapping,
    width, list(...), with_defaults(ggplot2::aes(colour = color), stack_rect),
    rect = GeomSwatch, rows = FALSE,
    prep = function(a, x) on_bar(x, a[a$to - a$from >= min_allele_frequency, ,
      drop = FALSE]))
}

stack_rect <- ggplot2::aes(ymin = low, ymax = high)

# The depth the band's bars grow from: the axis's bottom, which a log axis
# floors at one read
bar_base <- function(x) x$coverage_scale$domain[1] %||% 0

# A stack's `from` and `to` as the display places them, shares of its bar's
# drawn height from the bar's base, so a share means the same length on a
# log or symlog axis as on a linear one
on_bar <- function(x, d) {
  tr <- scales::as.transform(depth_transform(x$coverage_scale))
  base <- tr$transform(bar_base(x))
  at <- function(share) tr$inverse(base + share * (tr$transform(d$depth) -
    base))
  d$low <- at(d$from)
  d$high <- at(d$to)
  d
}

# Rects filled with their `colour`, as the segments these marks once were
# painted, so a pileup's fill scale keys its reads alone
GeomSwatch <- ggplot2::ggproto("GeomSwatch", ggplot2::GeomRect,
  default_aes = ggplot2::aes(colour = "grey30", alpha = NA, linewidth = 0,
    linetype = 1),
  draw_panel = function(self, data, panel_params, coord, ...) {
    data$fill <- data$colour
    data$colour <- NA
    ggplot2::ggproto_parent(ggplot2::GeomRect, self)$draw_panel(data,
      panel_params, coord, ...)
  },
  draw_key = function(data, params, size) {
    data$fill <- data$colour
    data$colour <- NA
    ggplot2::draw_key_polygon(data, params, size)
  })

#' @rdname geom_jbrowse_reads
#' @export
geom_jbrowse_interbase <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "interbase_coverage", "geom_jbrowse_interbase")
  jbrowse_layer(GeomInterbase, data, function(x) {
    d <- x$interbase_coverage
    full <- 2 * x$coverage_scale$domain[2]
    n <- nrow(d)
    d$hang <- if (full > 0) d$from / full else rep(0, n)
    d$hang_end <- if (full > 0) d$to / full else rep(0, n)
    d$top <- rep(px_to_pt(x$glyph$interbaseTop), n)
    d$inset <- rep(px_to_pt(x$glyph$coverageInset), n)
    d$width <- rep(px_to_pt(1), n)
    d
  }, mapping, ggplot2::aes(x = position, hang = hang, hang_end = hang_end,
    top = I(top), inset = I(inset), width = I(width), fill = color),
    list(...), keyed = "complete")
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
    d <- x$indicators
    d$width <- rep(px_to_pt(x$glyph$indicatorWidth), nrow(d))
    d$height <- rep(px_to_pt(x$glyph$indicatorHeight), nrow(d))
    d
  }, mapping, ggplot2::aes(x = position, fill = color, width = I(width),
    height = I(height)), list(...), keyed = FALSE)
}

# A display's px in grid's points, which are 72.27 to the inch
px_to_pt <- function(px) px / PX_PER_INCH * 72.27

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
    with_defaults(ggplot2::aes(fill = color, alpha = I(opacity)), stack_rect),
    rows = FALSE, keyed = FALSE, prep = on_bar_of)
}

on_bar_of <- function(d, x) on_bar(x, d)

junction_paths <- function(x, min_reads) {
  j <- shown_junctions(x, min_reads)
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
geom_jbrowse_arcs <- function(mapping = NULL, data = NULL, width = NULL,
                              ...) {
  data <- replacing(data, "arcs", "geom_jbrowse_arcs")
  domes <- function(x) connection_domes(x$arcs[x$arcs$shape == "arc", ,
    drop = FALSE], arc_top(x, width))
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
geom_jbrowse_arc_ticks <- function(mapping = NULL, data = NULL,
                                   width = NULL, ...) {
  data <- replacing(data, "arc_ticks", "geom_jbrowse_arc_ticks")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) {
    t <- x$arc_ticks
    t$bottom <- rep(if (is_cloud(x)) 1 else 0, nrow(t))
    t$top <- rep(arc_top(x, width), nrow(t))
    t
  }, mapping, ggplot2::aes(x = position, xend = position, y = bottom,
    yend = top, colour = color, linewidth = I(width / ggplot2::.pt)),
    list(...))
}

# The arc band's reach: the read cloud's widest insert, or as many bp as the
# band's pixels span in a figure `width` pixels across, so an arc rises as
# high as it is wide and one wider than the band runs off it, as the display
# draws it
arc_top <- function(x, width = NULL) {
  if (is_cloud(x)) max(2, x$arc_span_max)
  else x$settings$readConnectionsHeight / px_of(x, width)
}
