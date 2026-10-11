#' Stack tracks as a JBrowse view stacks them
#'
#' Each track's plot, or an alignments track's bands, one above the other in
#' the order given, on one x axis written under the last, with the keys
#' collected beside them: the figure a JBrowse view of these tracks shows.
#' Each piece is an ordinary ggplot, so `autoplot(reads, "pileup")` or a plot
#' of your own built from the same tables restyles one band.
#'
#' @param ... [jbrowse_features()], [jbrowse_wiggle()], [jbrowse_variants()]
#'   and [jbrowse_alignments()] results over the same regions, top to bottom,
#'   or lists of them
#' @return a patchwork
#' @export
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' locus <- "ctgA:1,001-9,000"
#' genes <- jbrowse_features(e("volvox.sort.gff3.gz"), locus)
#' signal <- jbrowse_wiggle(e("volvox_microarray_posneg.bw"), locus)
#' reads <- jbrowse_alignments(e("volvox-sorted.bam"), locus, color = "strand")
#' jbrowse_stack(genes, signal, reads)
jbrowse_stack <- function(...) {
  tracks <- list(...)
  tracks <- do.call(c, lapply(tracks, function(t)
    if (inherits(t, "jbrowse_result")) list(t) else t))
  stack_bands(lapply(tracks, track_bands))
}

# A track's bands and their relative heights: a feature track as tall as its
# laid-out pixels against a pileup band of about 300
track_bands <- function(x) {
  if (inherits(x, "jbrowse_alignments")) return(alignment_bands(x))
  height <- if (inherits(x, "jbrowse_features")) max(0.3, x$height / 100)
    else if (inherits(x, "jbrowse_variants")) max(1, nrow(x$rows) / 15)
    else if (isTRUE(x$row_layout)) 0.6 * nrow(x$sources)
    else 1
  list(plots = list(ggplot2::autoplot(x)), heights = height)
}

stack_bands <- function(bands) {
  plots <- do.call(c, lapply(bands, `[[`, "plots"))
  heights <- unlist(lapply(bands, `[[`, "heights"))
  last <- length(plots)
  plots[-last] <- lapply(plots[-last], function(p) p + no_x)
  patchwork::wrap_plots(plots, ncol = 1, heights = heights,
    guides = "collect")
}

# The x axis's title: the contig, or every region of an axis that lays
# several end to end
axis_label <- function(x) {
  r <- x$regions
  if (is.null(r)) return(sprintf("%s (bp)", x$region$chrom))
  sprintf("%s (bp, end to end)", paste(sprintf("%s:%s-%s", r$chrom,
    format(r$start + 1, big.mark = ",", trim = TRUE),
    format(r$end, big.mark = ",", trim = TRUE)), collapse = ", "))
}

# The x scale of an axis laying several regions end to end: round genomic
# positions in each region, about as many as its share of the axis earns,
# each labelled with its position in its own region, and none on an edge two
# regions face across their gap, where two labels would collide
genomic_x <- function(x) {
  r <- x$regions
  if (is.null(r)) return(NULL)
  w <- r$end - r$start
  at <- lapply(seq_len(nrow(r)), function(i) {
    b <- pretty(c(r$start[i], r$end[i]), n = max(2, round(6 * w[i] / sum(w))))
    lo <- if (i == 1) r$start[i] else r$start[i] + 1
    hi <- if (i == nrow(r)) r$end[i] else r$end[i] - 1
    b[b >= lo & b <= hi]
  })
  ggplot2::scale_x_continuous(
    breaks = unlist(lapply(seq_len(nrow(r)), function(i)
      at[[i]] - r$start[i] + r$offset[i])),
    labels = format(unlist(at), big.mark = ",", trim = TRUE,
      scientific = FALSE))
}
