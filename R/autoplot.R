# The colors JBrowse's alignments display draws in, from the engine's own
# palette: `read`, `forward`, `reverse`, `skip`, `deletion`, `insertion`,
# `softclip`, `hardclip`, `coverage` and `bases`, one per mismatch base
jbrowse_colors <- local({
  colors <- NULL
  function() {
    if (is.null(colors)) {
      read <- jsonlite::fromJSON(in_engine(function(ctx)
        ctx$eval("jbrowse.alignmentColors()")))
      read$bases <- unlist(read$bases)
      colors <<- read
    }
    colors
  }
})

# Each read as its aligned blocks: the read split at its splice skips
read_blocks <- function(x) {
  skips <- x$gaps[x$gaps$type == "skip", , drop = FALSE]
  reads <- x$reads
  if (!nrow(skips)) return(reads)
  by_read <- split(skips, skips$read)
  spliced <- reads$read %in% names(by_read)
  blocks <- lapply(which(spliced), function(i) {
    s <- by_read[[as.character(reads$read[i])]]
    s <- s[order(s$start), ]
    starts <- c(reads$start[i], s$end)
    ends <- c(s$start, reads$end[i])
    keep <- ends > starts
    r <- reads[rep(i, sum(keep)), , drop = FALSE]
    r$start <- starts[keep]
    r$end <- ends[keep]
    r
  })
  do.call(rbind, c(list(reads[!spliced, , drop = FALSE]), blocks))
}

# The points of each junction's arc: a half ellipse from its start to its end,
# its height following the junction's genomic span on the display's fixed log
# scale, so a long intron arcs higher than a short one whatever else is drawn.
junction_arcs <- function(junctions, points = 48) {
  if (!nrow(junctions)) {
    return(data.frame(id = integer(0), x = numeric(0), y = numeric(0)))
  }
  span <- junctions$end - junctions$start
  norm <- pmin(1, pmax(0, (log(pmax(1, span)) - log(50)) / (log(1e5) - log(50))))
  height <- 0.3 + 0.65 * norm
  t <- seq(0, pi, length.out = points)
  do.call(rbind, lapply(seq_len(nrow(junctions)), function(i) {
    data.frame(id = i,
      x = junctions$start[i] + span[i] * (1 - cos(t)) / 2,
      y = height[i] * sin(t))
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
# while the display's fadeLowFrequencyMismatches is off.
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
  scale <- x$coverage_scale
  p <- ggplot2::ggplot(x) +
    geom_jbrowse_coverage() +
    geom_jbrowse_alleles(min_allele_frequency = o$min_allele_frequency,
      width = o$width) +
    geom_jbrowse_modification_coverage(width = o$width)
  if (!isFALSE(x$settings$showInterbaseIndicators)) {
    p <- p + geom_jbrowse_interbase() + geom_jbrowse_indicators() +
      scale_fill_jbrowse(x, "interbase")
  }
  inset <- x$glyph$coverageInset
  band(p + scale_colour_jbrowse(x, "bases"), x, ylim = scale$domain,
    expand = c(top = TRUE, right = FALSE, bottom = TRUE, left = FALSE)) +
    ggplot2::scale_y_continuous(breaks = scale$ticks,
      transform = depth_transform(scale),
      expand = ggplot2::expansion(mult = inset /
        (x$settings$coverageHeight - 2 * inset))) +
    ggplot2::labs(y = "depth")
}

# The display's depth scale: log, linear, or its symlog,
# sign(x) log(1 + |x| / C), which scales' asinh-based pseudo-log is not
depth_transform <- function(scale) {
  c <- scale$symlogConstant
  switch(scale$type,
    log = "log10",
    symlog = scales::new_transform("symlog",
      function(x) sign(x) * log1p(abs(x) / c),
      function(y) sign(y) * c * expm1(abs(y))),
    "identity")
}

# The reference under the pileup, each base a lettered cell in its color
sequence_plot <- function(x, o) {
  xlim <- xlim_of(x)
  bases <- data.frame(position = seq(xlim[1], xlim[2] - 1))
  bases$base <- toupper(strsplit(jbrowse_sequence(x$reference$file,
    x$reference$chrom, xlim[1], xlim[2]), "")[[1]])
  bases$cell <- max(1, 1 / px_of(x, o$width))
  bases$tint <- ggplot2::alpha(jbrowse_colors()$bases[bases$base], 0.35)
  ggplot2::ggplot(bases) +
    ggplot2::geom_rect(ggplot2::aes(xmin = position, xmax = position + cell,
      ymin = 0, ymax = 1, fill = tint)) +
    ggplot2::scale_fill_identity() +
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
    scale_fill_jbrowse(x, "reads") +
    guides_jbrowse(x) +
    geom_jbrowse_qualities(width = px) +
    geom_jbrowse_skips() +
    geom_jbrowse_deletions(fade = o$fade, width = px) +
    geom_jbrowse_modifications(width = px) +
    geom_jbrowse_mismatches(fade = o$fade, width = px) +
    scale_colour_jbrowse(x, "bases") +
    geom_jbrowse_insertions(fade = o$fade, width = px) +
    geom_jbrowse_clips(fade = o$fade, width = px)
  if (o$letters) p <- p + geom_jbrowse_letters()
  band(p + ggplot2::scale_y_reverse(), x, rows = TRUE) + no_y
}

junctions_plot <- function(x, o) {
  shown <- x$junctions$color[x$junctions$count >= o$min_reads]
  key <- x$keys$junctions
  band(ggplot2::ggplot(x) +
    geom_jbrowse_junctions(min_reads = o$min_reads) +
    key_scale(ggplot2::scale_colour_identity,
      key[key$color %in% shown, , drop = FALSE], attr(key, "title")) +
    geom_jbrowse_junction_counts(min_reads = o$min_reads), x,
    ylim = c(0, 1.08)) +
    no_y
}

# The read-connection band as the display draws it under the coverage: arcs
# hanging from the top, or the read cloud on a log insert-size axis
arcs_plot <- function(x, o) {
  cloud <- is_cloud(x)
  p <- band(ggplot2::ggplot(x) +
    geom_jbrowse_arc_ticks(width = o$width) +
    geom_jbrowse_arcs(width = o$width) +
    scale_colour_jbrowse(x, "arcs") +
    if (cloud) ggplot2::scale_y_log10() else ggplot2::scale_y_reverse(),
    x, ylim = c(if (cloud) 1 else 0, arc_top(x, o$width)))
  if (cloud) p + ggplot2::labs(y = "insert size (bp)") else p + no_y
}

BAND_PLOTS <- list(junctions = junctions_plot, coverage = coverage_plot,
  arcs = arcs_plot, sequence = sequence_plot, pileup = pileup_plot)

#' Draw alignments as JBrowse's alignments track, or one band of it
#'
#' @param object a [jbrowse_alignments()] result
#' @param what `"track"`: the display's stack, any `above` on top,
#'   then splice-junction arcs where the reads splice, the coverage band,
#'   the read-connection band under `display$readConnections` and the
#'   pileup, one above the other, then any `below`, all on a shared axis, each band split into the `facet`'s sections, and each
#'   drawn where the display's `showCoverage`, `showPileup` and
#'   `showSashimiArcs` draw it; `"pileup"`: each read on its row, split at
#'   its splice gaps, filled as the display's `color` object paints it and
#'   keyed as the display keys it, with its mismatches, deletions,
#'   insertions, clips and base modifications; `"coverage"`: the depth, each
#'   mismatch position's alleles stacked in its bases' colors;
#'   `"junctions"`: the arcs alone; `"arcs"`: the read-connection band
#'   alone, mate-pair and split-read arcs or the read cloud on a log
#'   insert-size axis; `"sequence"`: the reference's bases alone
#' @param min_allele_frequency hide a coverage-band allele segment whose share
#'   of the depth there is below this; `NULL` takes the display's
#'   `coverageSnpMinFrequency`
#' @param min_reads hide a junction arc with fewer supporting reads;
#'   `NULL` takes the display's `minSashimiScore`
#' @param letters write each mismatch's base and each insertion's length on
#'   the pileup; `NULL` does so where the region is at most 300 bp, where a
#'   base is wide enough to read
#' @param fade_low_frequency fade a mismatch, indel or clip narrower than a
#'   pixel by its site's share of the depth, so sequencing errors recede on a
#'   deep pileup; `NULL` takes the display's `fadeLowFrequencyMismatches`
#' @param width the pixels the figure spans across its region, which that
#'   fade and the read-connection band's reach read; `NULL` takes 1,500
#' @param above,below tracks to draw above and below a `"track"`, each a
#'   ggjbrowse result over the same regions or a ggplot, or a list of them, as
#'   [jbrowse_stack()] stacks them
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot, to extend with layers, scales and themes or to stack with
#'   patchwork; a `"track"` is a patchwork of its bands
#' @exportS3Method ggplot2::autoplot
#' @examples
#' bam <- system.file("extdata", "spliced.bam", package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "ctgA:401-1,300",
#'   display = list(color = list(field = "strand")))
#' autoplot(a)
#' autoplot(a, "pileup") + ggplot2::theme_classic()
#' gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
#' autoplot(a, above = jbrowse_features(gff, "ctgA:401-1,300"))
autoplot.jbrowse_alignments <- function(object,
                                        what = c("track", "pileup", "coverage",
                                          "junctions", "arcs", "sequence"),
                                        min_allele_frequency = NULL,
                                        min_reads = NULL,
                                        letters = NULL,
                                        fade_low_frequency = NULL,
                                        width = NULL, above = NULL,
                                        below = NULL, ...) {
  rlang::check_dots_empty()
  what <- match.arg(what)
  object <- on_axis(object)
  o <- alignment_options(object, min_allele_frequency, min_reads,
    letters, fade_low_frequency, width)
  if (what != "track") return(BAND_PLOTS[[what]](object, o))
  above <- listed_tracks(above)
  below <- listed_tracks(below)
  same_regions(Filter(function(t) inherits(t, "jbrowse_result"),
    c(above, list(object), below)))
  stack_bands(c(lapply(above, band_of, object), list(alignment_bands(object,
    o)), lapply(below, band_of, object)))
}

# What the alignments bands draw: the caller's settings, else the display's
alignment_options <- function(object, min_allele_frequency = NULL,
                              min_reads = NULL, letters = NULL,
                              fade_low_frequency = NULL, width = NULL) {
  s <- object$settings %||% list()
  xlim <- xlim_of(object)
  list(
    min_allele_frequency = min_allele_frequency %||%
      s$coverageSnpMinFrequency %||% 0,
    min_reads = min_reads %||% s$minSashimiScore %||% 2,
    letters = letters %||% (xlim[2] - xlim[1] <= 300),
    fade = fade_low_frequency %||% s$fadeLowFrequencyMismatches %||% TRUE,
    width = width,
    coverage = !isFALSE(s$showCoverage),
    pileup = !isFALSE(s$showPileup),
    junctions = !isFALSE(s$showSashimiArcs))
}

# Each band an alignments track draws, top to bottom, by its height in the
# display's pixels: the coverage, sashimi and read-connection bands at their
# configured heights per section, the pileup its rows at the read pitch
alignment_px <- function(object, o) {
  s <- object$settings
  sections <- nrow(object$sections)
  h <- s$featureHeight
  rows <- sum(vapply(split(object$reads$row, object$reads$section),
    function(r) if (length(r)) max(r) + 1 else 0, 0))
  shown <- c(
    junctions = o$coverage && o$junctions &&
      any(object$junctions$count >= o$min_reads),
    coverage = o$coverage,
    arcs = nrow(object$arcs) + nrow(object$arc_ticks) > 0,
    pileup = o$pileup)
  c(junctions = s$sashimiArcsHeight * sections,
    coverage = s$coverageHeight * sections,
    arcs = s$readConnectionsHeight * sections,
    pileup = max(30, rows * (h + (h > 3))))[shown]
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
  for (name in setdiff(names(object), "sections")) {
    d <- object[[name]]
    if (!is.data.frame(d) || is.null(d$section)) next
    object[[name]] <- d[d$section == label, , drop = FALSE]
  }
  object$sections <- object$sections[object$sections$label == label, ,
    drop = FALSE]
  object
}
