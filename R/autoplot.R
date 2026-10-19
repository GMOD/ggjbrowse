# The engine's palette's color for each base, which the sequence band tints
# the reference's bases by
base_colors <- local({
  colors <- NULL
  function() {
    if (is.null(colors)) {
      colors <<- unlist(jsonlite::fromJSON(in_engine(function(ctx)
        ctx$eval("jbrowse.alignmentColors()")))$bases)
    }
    colors
  }
})

# Each read as its aligned blocks: the read split at its splice skips
read_blocks <- function(x) {
  reads <- x$reads
  skips <- x$gaps[x$gaps$type == "skip", , drop = FALSE]
  on <- match(paste(skips$read, skips$region), paste(reads$read, reads$region))
  skips <- skips[!is.na(on), , drop = FALSE]
  on <- on[!is.na(on)]
  if (!length(on)) return(reads)
  skips <- skips[order(on, skips$start), , drop = FALSE]
  n <- tabulate(on, nrow(reads))
  spliced <- which(n > 0)
  rows <- rep(spliced, n[spliced] + 1)
  first <- cumsum(c(1, utils::head(n[spliced] + 1, -1)))
  last <- first + n[spliced]
  starts <- reads$start[rows]
  starts[-first] <- skips$end
  ends <- reads$end[rows]
  ends[-last] <- skips$start
  blocks <- reads[rows, , drop = FALSE]
  blocks$start <- starts
  blocks$end <- ends
  rbind(reads[n == 0, , drop = FALSE], blocks[ends > starts, , drop = FALSE])
}

# The points of each junction's arc: a half ellipse from its start to its end,
# as high as the engine's `arcHeight`, a share of the band
junction_arcs <- function(junctions, points = 48) {
  if (!nrow(junctions)) {
    return(data.frame(id = integer(0), x = numeric(0), y = numeric(0)))
  }
  span <- junctions$end - junctions$start
  t <- seq(0, pi, length.out = points)
  do.call(rbind, lapply(seq_len(nrow(junctions)), function(i) {
    data.frame(id = i,
      x = junctions$start[i] + span[i] * (1 - cos(t)) / 2,
      y = junctions$arcHeight[i] * sin(t))
  }))
}

# The pileup's y axis is packed rows, which mean nothing to a reader
no_y <- ggplot2::theme(axis.text.y = ggplot2::element_blank(),
  axis.ticks.y = ggplot2::element_blank(),
  panel.grid.major.y = ggplot2::element_blank(),
  panel.grid.minor.y = ggplot2::element_blank())

no_x <- ggplot2::theme(axis.text.x = ggplot2::element_blank(),
  axis.title.x = ggplot2::element_blank())

# A band's sections stacked as the display stacks them, one panel each; rows
# keep one height across sections where `rows` says the y axis counts rows
by_section <- function(x, rows = FALSE) {
  if (identical(x$sections$label, "")) return(NULL)
  ggplot2::facet_grid(rows = ggplot2::vars(section), scales = "free_y",
    space = if (rows) "free_y" else "fixed")
}

# Section labels read across, as the display writes them
across <- ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0,
  hjust = 0))

# Each row named left of the panel at its `y`, in text as tall as the row
# the panel gives it at draw time, up to ROW_LABEL_PT, and none under
# MIN_ROW_LABEL_PT. The axis writes the same labels at full size, unseen, so
# the margin fits them; the caller's y scale gives it those labels.
row_names <- function(labels, y, gap_pt = 2.2) {
  list(ggplot2::layer(geom = GeomRowLabels, stat = "identity",
      position = "identity", data = data.frame(y = y, label = labels),
      mapping = ggplot2::aes(x = -Inf, y = y, label = label),
      inherit.aes = FALSE, show.legend = FALSE,
      params = list(rows = length(y), gap_pt = gap_pt, colour = "grey30")),
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = ROW_LABEL_PT,
      colour = NA)),
    axis_gap(gap_pt))
}

# Room of `pt` between a panel and its y axis labels, held by unseen ticks:
# ggplot2 widens the axis by its ticks' length, not by its text's margin
axis_gap <- function(pt) {
  ggplot2::theme(axis.ticks.y = ggplot2::element_line(colour = NA),
    axis.ticks.length.y = grid::unit(pt, "pt"))
}

ROW_LABEL_PT <- 8
MIN_ROW_LABEL_PT <- 3

GeomRowLabels <- ggplot2::ggproto("GeomRowLabels", ggplot2::GeomText,
  extra_params = c("na.rm", "rows", "gap_pt"),
  draw_panel = function(self, data, panel_params, coord, rows, gap_pt) {
    text <- ggplot2::ggproto_parent(ggplot2::GeomText, self)$draw_panel(data,
      panel_params, coord)
    grid::gTree(children = grid::gList(text), rows = rows, gap_pt = gap_pt,
      cl = "ggjbrowse_row_labels")
  })

#' @exportS3Method grid::makeContent
makeContent.ggjbrowse_row_labels <- function(x) {
  row_pt <- grid::convertHeight(grid::unit(1, "npc"), "pt",
    valueOnly = TRUE) / x$rows
  size <- min(ROW_LABEL_PT, 0.85 * row_pt)
  if (size < MIN_ROW_LABEL_PT) return(grid::setChildren(x, grid::gList()))
  text <- x$children[[1]]
  text$gp$fontsize <- size
  text$x <- grid::unit(0, "npc") - grid::unit(x$gap_pt, "pt")
  text$hjust <- 1
  grid::setChildren(x, grid::gList(text))
}

# The display's fade of a sub-pixel mark: at `px` pixels wide it blends from
# its own width up to opaque by the share of the depth its site holds, which
# the engine zeroes under the display's threshold (the shaders'
# frequencyAlpha); a mark a pixel or wider never fades, nor does anything
# while the display's pileup$fadeLowFrequency is off.
mark_alpha <- function(px, fade, on) {
  if (!on) return(rep(1, length(fade)))
  px <- rep_len(px, length(fade))
  ifelse(px < 1, px + fade * (1 - px), 1)
}

# Each band is the geoms of its tables, the display's keys as scales, its y
# axis, and its sections as facets
band <- function(p, x, rows = FALSE, ylim = NULL, expand = FALSE) {
  p +
    ggplot2::coord_cartesian(xlim = xlim_of(x), ylim = ylim, expand = expand) +
    by_section(x, rows) +
    genomic_x(x) +
    ggplot2::labs(x = axis_label(x), y = NULL) +
    ggplot2::theme_minimal() +
    across
}

# The depth axis the display draws: its autoscaled domain and ticks, shared
# by every section, inset at both ends by the room the display keeps for its
# labels, the triangles sitting in the top one
coverage_plot <- function(x, o) {
  scale <- x$y_scale
  interbase <- !isFALSE(x$display$coverage$showInterbaseIndicators)
  p <- ggplot2::ggplot(x) +
    geom_jbrowse_coverage() +
    geom_jbrowse_alleles(width = o$width) +
    geom_jbrowse_modification_coverage(width = o$width) +
    guides_jbrowse(x, c("bases", if (interbase) "interbase"), ramps = FALSE)
  if (interbase) p <- p + geom_jbrowse_interbase() + geom_jbrowse_indicators()
  inset <- x$glyph$coverageInset
  band(p, x, ylim = scale$domain,
    expand = c(top = TRUE, right = FALSE, bottom = TRUE, left = FALSE)) +
    ggplot2::scale_y_continuous(breaks = scale$ticks,
      transform = score_transform(scale),
      expand = ggplot2::expansion(mult = inset /
        (x$display$coverage$height - 2 * inset))) +
    ggplot2::labs(y = "depth")
}

# The reference under the pileup, each base a lettered cell in its color
sequence_plot <- function(x, o) {
  if (is.null(x$sequence)) {
    stop("the sequence band draws the reference's bases: read the alignments ",
      "with `assembly =`", call. = FALSE)
  }
  xlim <- xlim_of(x)
  r <- x$regions
  bases <- do.call(rbind, lapply(seq_len(nrow(r)), function(i) {
    base <- toupper(strsplit(x$sequence[i] %|NA|% "", "")[[1]])
    data.frame(position = r$offset[i] + seq_along(base) - 1, base = base)
  }))
  bases$cell <- max(1, 1 / px_of(x, o$width))
  bases$tint <- ggplot2::alpha(base_colors()[bases$base], 0.35)
  ggplot2::ggplot(bases) +
    ggplot2::geom_rect(ggplot2::aes(xmin = position, xmax = position + cell,
      ymin = 0, ymax = 1, fill = I(tint))) +
    ggplot2::geom_text(ggplot2::aes(x = position + 0.5, y = 0.5, label = base),
      size = 2.6) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(x = axis_label(x), y = NULL) +
    ggplot2::theme_minimal() +
    no_y
}

pileup_plot <- function(x, o) {
  px <- o$width
  p <- ggplot2::ggplot(x) +
    geom_jbrowse_reads() +
    geom_jbrowse_overlaps() +
    guides_jbrowse(x, c("reads", "bases")) +
    geom_jbrowse_qualities(width = px) +
    geom_jbrowse_skips() +
    geom_jbrowse_deletions(width = px) +
    geom_jbrowse_modifications(width = px) +
    geom_jbrowse_mismatches(width = px) +
    geom_jbrowse_insertions(width = px) +
    geom_jbrowse_clips(width = px) +
    geom_jbrowse_chains() +
    geom_jbrowse_deletion_lengths(width = px)
  if (o$letters) p <- p + geom_jbrowse_letters()
  band(p + ggplot2::scale_y_reverse(), x, rows = TRUE) + no_y
}

shown <- function(d) d$shown %||% rep(TRUE, nrow(d))

shown_junctions <- function(x) {
  j <- summary_table(x, "junctions")
  j[shown(j), , drop = FALSE]
}

junctions_plot <- function(x, o) {
  p <- ggplot2::ggplot(x) +
    geom_jbrowse_junctions() +
    guides_jbrowse(x, "junctions", ramps = FALSE)
  if (isTRUE(x$display$sashimi$showCounts)) {
    p <- p + geom_jbrowse_junction_counts()
  }
  band(p, x, ylim = c(0, 1.08)) + no_y
}

# The read-connection band as the display draws it under the coverage: arcs
# hanging from the top, or the read cloud on a log insert-size axis
arcs_plot <- function(x, o) {
  cloud <- is_cloud(x)
  p <- band(ggplot2::ggplot(x) +
    geom_jbrowse_arc_ticks(width = o$width) +
    geom_jbrowse_arcs(width = o$width) +
    guides_jbrowse(x, "arcs", ramps = FALSE) +
    if (cloud) ggplot2::scale_y_log10() else ggplot2::scale_y_reverse(),
    x, ylim = c(if (cloud) 1 else 0, arc_top(x, o$width)))
  if (cloud) p + ggplot2::labs(y = "insert size (bp)") else p + no_y
}

BAND_PLOTS <- list(junctions = junctions_plot, coverage = coverage_plot,
  arcs = arcs_plot, sequence = sequence_plot, pileup = pileup_plot)

#' Draw alignments as JBrowse's alignments track, or one band of it
#'
#' Each band draws what the result says: the junctions whose `shown` is
#' true, the allele segments whose `shown` is true, and each mark faded by
#' its `fade` under the display's `pileup$fadeLowFrequency`. Edit those
#' columns to draw otherwise, as in
#' `a$junctions$shown <- a$junctions$count >= 5`. [jbrowse_stack()] stacks
#' other tracks over the track.
#'
#' @param object a [jbrowse_alignments()] result
#' @param what `"track"`: the display's stack, splice-junction arcs where
#'   the reads splice, the coverage band, the read-connection band under
#'   `display$readConnections` and the pileup, one above the other on a
#'   shared axis, each band split into the `facet`'s sections, and each
#'   drawn where the display's `coverage`, `pileup` and `sashimi` bands
#'   `show`; `"pileup"`: each read on its row, split at
#'   its splice gaps, filled as the display's `color` object paints it and
#'   keyed as the display keys it, with its mismatches, deletions,
#'   insertions, clips and base modifications; `"coverage"`: the depth, each
#'   mismatch position's alleles stacked in its bases' colors;
#'   `"junctions"`: the arcs alone; `"arcs"`: the read-connection band
#'   alone, mate-pair and split-read arcs or the read cloud on a log
#'   insert-size axis; `"sequence"`: the reference's bases alone
#' @param letters write each mismatch's base and each insertion's length on
#'   the pileup; `NULL` does so where the display does, at 6.5 px a base or
#'   more on reads tall enough for text
#' @param width the pixels the figure spans across its region, which the
#'   sub-pixel fade and the read-connection band's reach read; `NULL` takes
#'   the `width` the track was read at
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot, to extend with layers, scales and themes or to stack with
#'   patchwork; a `"track"` is a patchwork of its bands
#' @exportS3Method ggplot2::autoplot
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' a <- jbrowse_alignments(e("gapdh_direct_rna.bam"),
#'   "chr12:6,534,201-6,538,700",
#'   display = list(color = list(field = "strand")))
#' autoplot(a)
#' autoplot(a, "coverage") + ggplot2::theme_classic()
autoplot.jbrowse_alignments <- function(object,
                                        what = c("track", "pileup", "coverage",
                                          "junctions", "arcs", "sequence"),
                                        letters = NULL, width = NULL, ...) {
  rlang::check_dots_empty()
  what <- match.arg(what)
  object <- on_axis(drawn_at(object, width))
  o <- alignment_options(object, letters)
  if (what != "track") return(BAND_PLOTS[[what]](object, o))
  stack_bands(list(alignment_bands(object, o)))
}

MIN_PX_PER_BP_FOR_TEXT <- 6.5

alignment_options <- function(object, letters = NULL) {
  s <- object$display
  tall <- (object$read_pitch$featureHeight %||% 0) >=
    (object$glyph$label$minHeight %||% Inf)
  list(
    letters = letters %||% (tall && px_of(object, NULL) >= MIN_PX_PER_BP_FOR_TEXT),
    coverage = !isFALSE(s$coverage$show),
    pileup = !isFALSE(s$pileup$show),
    junctions = !isFALSE(s$sashimi$show))
}

# Each band an alignments track draws, top to bottom, by its height in the
# display's pixels: the coverage, sashimi and read-connection bands at their
# configured heights, the pileup its rows at the read pitch
alignment_px <- function(object, o) {
  s <- object$display
  pitch <- object$read_pitch$featureHeight + object$read_pitch$featureSpacing
  rows <- sum(vapply(split(object$reads$row, object$reads$section),
    function(r) if (length(r)) max(r) + 1 else 0, 0))
  shown <- c(
    junctions = o$junctions && nrow(shown_junctions(object)) > 0,
    coverage = o$coverage,
    arcs = nrow(object$arcs) + nrow(object$arc_ticks) > 0,
    pileup = o$pileup)
  c(junctions = s$sashimi$height,
    coverage = s$coverage$height,
    arcs = s$readConnections$height,
    pileup = max(30, rows * pitch))[shown]
}

# A section's bands one after another, as the display stacks each section's
# coverage over its own reads, the section named once, on its first band
alignment_bands <- function(object, o = alignment_options(object)) {
  object <- on_axis(object)
  labels <- object$sections$label
  sections <- if (length(labels) > 1) {
    lapply(labels, function(l) section_of(object, l))
  } else {
    list(object)
  }
  bands <- lapply(sections, function(x) {
    px <- alignment_px(x, o)
    plots <- lapply(names(px), function(b) BAND_PLOTS[[b]](x, o))
    plots[-1] <- lapply(plots[-1], function(p)
      p + ggplot2::theme(strip.text.y = ggplot2::element_blank()))
    list(plots = plots, heights = unname(px))
  })
  list(plots = do.call(c, lapply(bands, `[[`, "plots")),
    heights = unlist(lapply(bands, `[[`, "heights")))
}

# The result with one section's rows: its reads, their marks and its coverage
section_of <- function(object, label) {
  cls <- class(object)
  object <- unclass(object)
  for (name in setdiff(names(object), "sections")) {
    d <- object[[name]]
    if (!is.data.frame(d) || is.null(d$section)) next
    object[[name]] <- d[d$section == label, , drop = FALSE]
  }
  object$sections <- object$sections[object$sections$label == label, ,
    drop = FALSE]
  structure(object, class = cls)
}
