#' The colors JBrowse's alignments display draws in
#'
#' `read` (the plain fill), `forward` and `reverse` (each strand, a splice
#' junction's too), `skip`, `deletion`, `insertion`, `softclip`, `hardclip`,
#' `coverage` and `bases` (a named vector, one per mismatch base), read from
#' the engine's own palette.
#' @export
jbrowse_colors <- local({
  colors <- NULL
  function() {
    if (is.null(colors)) {
      read <- jsonlite::fromJSON(engine()$eval("jbrowse.alignmentColors()"))
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
  reads$read <- seq_len(nrow(reads))
  if (!nrow(skips)) return(reads)
  by_read <- split(skips, skips$read)
  spliced <- reads$read %in% names(by_read)
  blocks <- lapply(which(spliced), function(i) {
    s <- by_read[[as.character(i)]]
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

# A base's width on a figure about 170 mm wide, never thinner than a tick
cell_width <- function(xlim) max(0.3, 170 / (xlim[2] - xlim[1]))

# One scale for every base drawn, the same in every band so they share one
# key: the modification coverage's own colors ride along off the key
base_scale <- function(colors, x) {
  own <- unique(x$modification_coverage$color)
  ggplot2::scale_colour_manual(values = c(colors$bases, stats::setNames(own, own)),
    name = "Base", breaks = names(colors$bases))
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

# The coverage band: the depth, with the display's own stacks over it, each
# mismatch position's alleles in their base colors (a segment under
# `min_allele_frequency` left out, as the display's coverageSnpMinFrequency
# hides it) and, under the modifications layer, each position's modification
# states in the colors and opacities the display paints them.
coverage_plot <- function(x, colors, min_allele_frequency, xlim, xlab) {
  a <- x$allele_coverage
  alleles <- a[a$to - a$from >= min_allele_frequency, , drop = FALSE]
  mods <- x$modification_coverage
  ggplot2::ggplot(x$coverage) +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = 0,
      ymax = depth), fill = colors$coverage) +
    ggplot2::geom_segment(ggplot2::aes(x = position + 0.5, xend = position + 0.5,
      y = depth * from, yend = depth * to, colour = base), alleles,
      linewidth = cell_width(xlim), lineend = "butt") +
    base_scale(colors, x) +
    ggplot2::geom_segment(ggplot2::aes(x = position + 0.5, xend = position + 0.5,
      y = depth * from, yend = depth * to, colour = color, alpha = opacity),
      mods, linewidth = cell_width(xlim), lineend = "butt",
      show.legend = FALSE) +
    ggplot2::scale_alpha_identity() +
    ggplot2::coord_cartesian(xlim = xlim, expand = FALSE) +
    by_section(x) +
    ggplot2::labs(x = xlab, y = "depth") +
    ggplot2::theme_minimal()
}

# The reference under the pileup, each base a lettered cell in its color
sequence_plot <- function(x, colors, xlim, xlab) {
  bases <- data.frame(position = seq(xlim[1], xlim[2] - 1))
  bases$base <- toupper(strsplit(jbrowse_sequence(x$reference$file,
    x$reference$chrom, xlim[1], xlim[2]), "")[[1]])
  bases$tint <- ggplot2::alpha(colors$bases[bases$base], 0.35)
  ggplot2::ggplot(bases) +
    ggplot2::geom_segment(ggplot2::aes(x = position + 0.5, xend = position + 0.5,
      y = 0, yend = 1, colour = tint), linewidth = cell_width(xlim),
      lineend = "butt") +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(ggplot2::aes(x = position + 0.5, y = 0.5, label = base),
      size = 2.6) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(x = xlab, y = NULL) +
    ggplot2::theme_minimal() +
    no_y
}

pileup_plot <- function(x, colors, letters, xlim, xlab, px, fade) {
  blocks <- read_blocks(x)
  legend <- x$legend
  keyed <- legend[!duplicated(legend$color), , drop = FALSE]
  mods <- x$modifications
  gaps <- x$gaps
  cell <- cell_width(xlim)
  deletions <- gaps[gaps$type == "deletion", , drop = FALSE]
  deletions$alpha <- mark_alpha(((deletions$end - deletions$start) * px)^2,
    deletions$fade, fade)
  mismatches <- x$mismatches
  mismatches$alpha <- mark_alpha(px, mismatches$fade, fade)
  ins <- x$insertions
  ins$alpha <- mark_alpha(px^2, ins$fade, fade)
  clips <- x$clips
  clips$alpha <- mark_alpha(px, clips$fade, fade)
  soft <- clips[clips$type == "soft", , drop = FALSE]
  hard <- clips[clips$type == "hard", , drop = FALSE]
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = row + 0.1,
      ymax = row + 0.9, fill = color), blocks) +
    ggplot2::scale_fill_identity(name = attr(legend, "title"),
      breaks = keyed$color, labels = keyed$label, guide = "legend") +
    ggplot2::geom_segment(ggplot2::aes(x = start, xend = end, y = row + 0.5,
      yend = row + 0.5), gaps[gaps$type == "skip", , drop = FALSE],
      colour = colors$skip, linewidth = 0.3) +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = row + 0.35,
      ymax = row + 0.65, alpha = alpha), deletions, fill = colors$deletion) +
    ggplot2::geom_rect(ggplot2::aes(xmin = position, xmax = position + 1,
      ymin = row + 0.1, ymax = row + 0.9, fill = color, alpha = opacity), mods) +
    ggplot2::geom_segment(ggplot2::aes(x = position + 0.5, xend = position + 0.5,
      y = row + 0.1, yend = row + 0.9, colour = base, alpha = alpha),
      mismatches, linewidth = cell, lineend = "butt") +
    base_scale(colors, x) +
    ggplot2::geom_segment(ggplot2::aes(x = position, xend = position,
      y = row + 0.05, yend = row + 0.95, alpha = alpha), ins,
      colour = colors$insertion, linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = position, xend = position,
      y = row + 0.05, yend = row + 0.95, alpha = alpha), soft,
      colour = colors$softclip, linewidth = 0.7) +
    ggplot2::geom_segment(ggplot2::aes(x = position, xend = position,
      y = row + 0.05, yend = row + 0.95, alpha = alpha), hard,
      colour = colors$hardclip, linewidth = 0.7) +
    ggplot2::scale_alpha_identity()
  if (letters) {
    p <- p +
      ggplot2::geom_text(ggplot2::aes(x = position + 0.5, y = row + 0.5,
        label = base), mismatches, size = 2.6, colour = "black") +
      ggplot2::geom_label(ggplot2::aes(x = position, y = row + 0.5,
        label = length), ins, size = 2.2, fill = colors$insertion,
        colour = "white", label.padding = grid::unit(0.1, "lines"),
        linewidth = 0)
  }
  p +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(xlim = xlim, expand = FALSE) +
    by_section(x, rows = TRUE) +
    ggplot2::labs(x = xlab, y = NULL) +
    ggplot2::theme_minimal() +
    no_y
}

junctions_plot <- function(x, colors, min_junction_reads, xlim, xlab) {
  j <- x$junctions[x$junctions$count >= min_junction_reads, , drop = FALSE]
  arcs <- junction_arcs(j)
  strand <- factor(j$strand, levels = c(1, -1, 0))
  arcs$strand <- strand[arcs$id]
  arcs$section <- j$section[arcs$id]
  arcs$width <- pmin(4, pmax(0.4, log(j$count[arcs$id] + 1) / 2))
  # the count sits at the highest point of the arc that is in view
  shown <- arcs[arcs$x >= xlim[1] & arcs$x <= xlim[2], , drop = FALSE]
  top <- shown[order(shown$id, -shown$y), ]
  top <- top[!duplicated(top$id), , drop = FALSE]
  apex <- data.frame(x = top$x, y = top$y, label = j$count[top$id],
    section = top$section)
  ggplot2::ggplot() +
    ggplot2::geom_path(ggplot2::aes(x = x, y = y, group = id, colour = strand,
      linewidth = width), arcs, lineend = "round") +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_colour_manual(values = c("1" = colors$forward,
      "-1" = colors$reverse, "0" = colors$read),
      labels = c("1" = "Forward", "-1" = "Reverse", "0" = "Unknown"),
      name = "Junction strand", drop = TRUE) +
    ggplot2::geom_label(ggplot2::aes(x = x, y = y, label = label), apex,
      size = 2.6, linewidth = 0, label.padding = grid::unit(0.15, "lines")) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1.08), expand = FALSE) +
    by_section(x) +
    ggplot2::labs(x = xlab, y = NULL) +
    ggplot2::theme_minimal() +
    no_y
}

# Each connection's arc as points: a half circle from foot to foot, its
# height the genomic radius, as the display's band rises an arc as high as it
# is wide
connection_domes <- function(arcs, points = 48) {
  t <- seq(0, pi, length.out = points)
  do.call(rbind, lapply(seq_len(nrow(arcs)), function(i) {
    data.frame(id = i, x = arcs$start[i] + (arcs$end[i] - arcs$start[i]) *
      (1 - cos(t)) / 2, y = arcs$y[i] * sin(t), color = arcs$color[i],
      width = arcs$width[i], section = arcs$section[i])
  })) %||% data.frame(id = integer(0), x = numeric(0), y = numeric(0),
    color = character(0), width = numeric(0), section = arcs$section)
}

# The read-connection band as the display draws it under the coverage: each
# mate pair's or split read's arc hanging from the top, or, as the read
# cloud, each pair a bar at its insert size on a log axis with a square at
# each foot, split reads dashed; interchromosomal ticks span the band
arcs_plot <- function(x, xlim, xlab) {
  arcs <- x$arcs
  ticks <- x$arc_ticks
  legend <- x$arc_legend
  keyed <- legend[!duplicated(legend$color), , drop = FALSE]
  cloud <- x$read_connections == "cloud"
  top <- if (cloud) max(2, x$arc_span_max) else max(1, arcs$y)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = position, xend = position,
      y = if (cloud) 1 else 0, yend = top, colour = color,
      linewidth = width / ggplot2::.pt), ticks)
  if (cloud) {
    arcs$y <- pmax(1, arcs$y)
    feet <- rbind(transform(arcs, x = start), transform(arcs, x = end))
    p <- p +
      ggplot2::geom_segment(ggplot2::aes(x = start, xend = end, y = y,
        yend = y, colour = color, alpha = opacity,
        linewidth = width / ggplot2::.pt,
        linetype = shape == "flatSplit"), arcs) +
      ggplot2::scale_linetype_manual(values = c("FALSE" = "solid",
        "TRUE" = "22"), guide = "none") +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = color), feet,
        shape = 15, size = 1.3) +
      ggplot2::scale_y_log10() +
      ggplot2::scale_alpha_identity()
  } else {
    p <- p +
      ggplot2::geom_path(ggplot2::aes(x = x, y = y, group = id,
        colour = color, linewidth = width / ggplot2::.pt),
        connection_domes(arcs)) +
      ggplot2::scale_y_reverse()
  }
  p +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_colour_identity(name = attr(legend, "title"),
      breaks = keyed$color, labels = keyed$label, guide = "legend") +
    ggplot2::coord_cartesian(xlim = xlim,
      ylim = if (cloud) c(1, top) else c(0, top), expand = FALSE) +
    by_section(x) +
    ggplot2::labs(x = xlab, y = if (cloud) "insert size (bp)" else NULL) +
    ggplot2::theme_minimal() +
    if (cloud) NULL else no_y
}

#' Draw alignments as JBrowse's alignments track, or one band of it
#'
#' @param object a [jbrowse_alignments()] result
#' @param what `"track"`: the display's stack, any `features` tracks on top,
#'   then splice-junction arcs where the reads splice, the coverage band,
#'   the read-connection band under `read_connections`, the reference's bases
#'   where they are wide enough to read and the pileup, one above the other
#'   on a shared axis, each band split into the `facet`'s sections; `"pileup"`: each read on its row, split at its splice gaps, filled
#'   as the display's `color` object paints it and keyed as the display keys
#'   it, with its mismatches, deletions, insertions, clips and base
#'   modifications; `"coverage"`: the depth, each mismatch position's alleles
#'   stacked in its bases' colors; `"junctions"`: the arcs alone;
#'   `"arcs"`: the read-connection band alone, mate-pair and split-read arcs
#'   or the read cloud on a log insert-size axis; `"sequence"`: the
#'   reference's bases alone
#' @param min_allele_frequency hide a coverage-band allele segment whose share
#'   of the depth there is below this, the display's `coverageSnpMinFrequency`
#' @param min_junction_reads hide a junction arc with fewer supporting reads,
#'   the display's `minSashimiScore`
#' @param letters write each mismatch's base and each insertion's length on
#'   the pileup; `NULL` does so where the region is at most 300 bp, where a
#'   base is wide enough to read
#' @param fade_low_frequency fade a mismatch, indel or clip narrower than a
#'   pixel by its site's share of the depth, as the display does
#'   (`fadeLowFrequencyMismatches`), so sequencing errors recede on a deep
#'   pileup; `FALSE` draws every mark opaque
#' @param px_per_bp how many pixels a base spans, which that fade reads;
#'   `NULL` takes a figure about 1200 pixels wide
#' @param features feature tracks to draw above a `"track"`, each a
#'   [jbrowse_features()] result over the same region, or a list of them
#' @param ... unused
#' @return a ggplot, to extend with layers, scales and themes or to stack with
#'   patchwork; a `"track"` is a patchwork of its bands
#' @exportS3Method ggplot2::autoplot
#' @examples
#' bam <- system.file("extdata", "spliced.bam", package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "ctgA", 400, 1300, color = "strand")
#' autoplot(a)
#' autoplot(a, "pileup") + ggplot2::theme_classic()
#' gff <- system.file("extdata", "volvox.sort.gff3.gz", package = "ggjbrowse")
#' autoplot(a, features = jbrowse_features(gff, "ctgA", 400, 1300))
autoplot.jbrowse_alignments <- function(object,
                                        what = c("track", "pileup", "coverage",
                                          "junctions", "arcs", "sequence"),
                                        min_allele_frequency = 0,
                                        min_junction_reads = 2,
                                        letters = NULL,
                                        fade_low_frequency = TRUE,
                                        px_per_bp = NULL, features = NULL,
                                        ...) {
  what <- match.arg(what)
  colors <- jbrowse_colors()
  xlab <- sprintf("%s (bp)", object$region$chrom)
  xlim <- c(object$region$start, object$region$end)
  if (is.null(letters)) letters <- xlim[2] - xlim[1] <= 300
  if (is.null(px_per_bp)) px_per_bp <- 1200 / (xlim[2] - xlim[1])
  pileup <- function() pileup_plot(object, colors, letters, xlim, xlab,
    px_per_bp, fade_low_frequency)
  if (what == "coverage") {
    return(coverage_plot(object, colors, min_allele_frequency, xlim, xlab))
  }
  if (what == "pileup") {
    return(pileup())
  }
  if (what == "junctions") {
    return(junctions_plot(object, colors, min_junction_reads, xlim, xlab))
  }
  if (what == "arcs") {
    return(arcs_plot(object, xlim, xlab))
  }
  if (what == "sequence") {
    return(sequence_plot(object, colors, xlim, xlab))
  }
  connected <- nrow(object$arcs) + nrow(object$arc_ticks) > 0
  sections <- nrow(object$sections)
  spliced <- any(object$junctions$count >= min_junction_reads)
  sequenced <- letters && !is.null(object$reference)
  if (inherits(features, "jbrowse_features")) features <- list(features)
  bands <- c(
    lapply(features, function(f) autoplot(f) + no_x),
    if (spliced) list(junctions_plot(object, colors, min_junction_reads, xlim,
      xlab) + no_x),
    list(coverage_plot(object, colors, min_allele_frequency, xlim, xlab) + no_x),
    if (connected) list(arcs_plot(object, xlim, xlab) + no_x),
    if (sequenced) list(sequence_plot(object, colors, xlim, xlab) + no_x),
    list(pileup()))
  # a feature track's height in pixels, against a pileup band of about 300
  heights <- c(vapply(features, function(f) max(0.3, f$height / 100), 1),
    if (spliced) 0.8 * sections, sections, if (connected) 0.8 * sections,
    if (sequenced) 0.25, 3)
  patchwork::wrap_plots(bands, ncol = 1, heights = heights, guides = "collect")
}
