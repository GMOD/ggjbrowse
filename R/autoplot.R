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
band <- function(p, x, xlim, xlab, rows = FALSE, ylim = NULL) {
  p +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    by_section(x, rows) +
    ggplot2::labs(x = xlab, y = NULL) +
    ggplot2::theme_minimal()
}

coverage_plot <- function(x, min_allele_frequency, xlim, xlab) {
  band(ggplot2::ggplot(x) +
    geom_jbrowse_coverage() +
    geom_jbrowse_alleles(min_allele_frequency = min_allele_frequency) +
    geom_jbrowse_modification_coverage() +
    scale_colour_jbrowse(x, "bases"), x, xlim, xlab) +
    ggplot2::labs(y = "depth")
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

pileup_plot <- function(x, letters, xlim, xlab, px, fade) {
  p <- ggplot2::ggplot(x) +
    geom_jbrowse_reads() +
    scale_fill_jbrowse(x, "reads") +
    geom_jbrowse_skips() +
    geom_jbrowse_deletions(fade = fade, px_per_bp = px) +
    geom_jbrowse_modifications() +
    geom_jbrowse_mismatches(fade = fade, px_per_bp = px) +
    scale_colour_jbrowse(x, "bases") +
    geom_jbrowse_insertions(fade = fade, px_per_bp = px) +
    geom_jbrowse_clips(fade = fade, px_per_bp = px)
  if (letters) p <- p + geom_jbrowse_letters()
  band(p + ggplot2::scale_y_reverse(), x, xlim, xlab, rows = TRUE) + no_y
}

junctions_plot <- function(x, min_junction_reads, xlim, xlab) {
  shown <- x$junctions$color[x$junctions$count >= min_junction_reads]
  key <- x$keys$junctions
  band(ggplot2::ggplot(x) +
    geom_jbrowse_junctions(min_reads = min_junction_reads) +
    key_scale(ggplot2::scale_colour_identity,
      key[key$color %in% shown, , drop = FALSE], attr(key, "title")) +
    geom_jbrowse_junction_counts(min_reads = min_junction_reads), x, xlim,
    xlab, ylim = c(0, 1.08)) +
    no_y
}

# The read-connection band as the display draws it under the coverage: arcs
# hanging from the top, or the read cloud on a log insert-size axis
arcs_plot <- function(x, xlim, xlab) {
  cloud <- is_cloud(x)
  p <- band(ggplot2::ggplot(x) +
    geom_jbrowse_arc_ticks() +
    geom_jbrowse_arcs() +
    scale_colour_jbrowse(x, "arcs") +
    if (cloud) ggplot2::scale_y_log10() else ggplot2::scale_y_reverse(),
    x, xlim, xlab, ylim = c(if (cloud) 1 else 0, arc_top(x)))
  if (cloud) p + ggplot2::labs(y = "insert size (bp)") else p + no_y
}

#' Draw alignments as JBrowse's alignments track, or one band of it
#'
#' @param object a [jbrowse_alignments()] result
#' @param what `"track"`: the display's stack, any `tracks` on top,
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
#' @param tracks tracks to draw above a `"track"`, each a
#'   [jbrowse_features()] or [jbrowse_wiggle()] result over the same region,
#'   or a list of them
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
#' autoplot(a, tracks = jbrowse_features(gff, "ctgA", 400, 1300))
autoplot.jbrowse_alignments <- function(object,
                                        what = c("track", "pileup", "coverage",
                                          "junctions", "arcs", "sequence"),
                                        min_allele_frequency = 0,
                                        min_junction_reads = 2,
                                        letters = NULL,
                                        fade_low_frequency = TRUE,
                                        px_per_bp = NULL, tracks = NULL,
                                        ...) {
  what <- match.arg(what)
  xlab <- sprintf("%s (bp)", object$region$chrom)
  xlim <- c(object$region$start, object$region$end)
  if (is.null(letters)) letters <- xlim[2] - xlim[1] <= 300
  pileup <- function() pileup_plot(object, letters, xlim, xlab,
    px_per_bp, fade_low_frequency)
  if (what == "coverage") {
    return(coverage_plot(object, min_allele_frequency, xlim, xlab))
  }
  if (what == "pileup") {
    return(pileup())
  }
  if (what == "junctions") {
    return(junctions_plot(object, min_junction_reads, xlim, xlab))
  }
  if (what == "arcs") {
    return(arcs_plot(object, xlim, xlab))
  }
  if (what == "sequence") {
    return(sequence_plot(object, jbrowse_colors(), xlim, xlab))
  }
  connected <- nrow(object$arcs) + nrow(object$arc_ticks) > 0
  sections <- nrow(object$sections)
  spliced <- any(object$junctions$count >= min_junction_reads)
  sequenced <- letters && !is.null(object$reference)
  if (!is.null(tracks$region)) tracks <- list(tracks)
  bands <- c(
    lapply(tracks, function(t) autoplot(t) + no_x),
    if (spliced) list(junctions_plot(object, min_junction_reads, xlim,
      xlab) + no_x),
    list(coverage_plot(object, min_allele_frequency, xlim, xlab) + no_x),
    if (connected) list(arcs_plot(object, xlim, xlab) + no_x),
    if (sequenced) list(sequence_plot(object, jbrowse_colors(), xlim, xlab) +
      no_x),
    list(pileup()))
  # a feature track's height in pixels, against a pileup band of about 300
  heights <- c(vapply(tracks, function(t)
      if (is.null(t$height)) 1 else max(0.3, t$height / 100), 1),
    if (spliced) 0.8 * sections, sections, if (connected) 0.8 * sections,
    if (sequenced) 0.25, 3)
  patchwork::wrap_plots(bands, ncol = 1, heights = heights, guides = "collect")
}
