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
#'   or lists of them; or a ggplot of your own, drawn on the stack's x axis:
#'   genomic positions over one region, and over several the axis that
#'   `ggplot(<result>)` draws on. A track's name, as an argument's or a list
#'   element's, titles its first band, as the view labels a track:
#'   `jbrowse_stack(genes, CHH = reads)`
#' @return a patchwork
#' @export
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' locus <- "chr12:6,534,201-6,538,700"
#' genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"), locus)
#' signal <- jbrowse_wiggle(e("gapdh_phylop100way.bw"), locus)
#' reads <- jbrowse_alignments(e("gapdh_direct_rna.bam"), locus)
#' jbrowse_stack(genes, signal, reads)
#' jbrowse_stack(GENCODE = genes, phyloP = signal, reads)
jbrowse_stack <- function(...) {
  args <- list(...)
  tracks <- do.call(c, unname(Map(named_tracks, args,
    names(args) %||% rep("", length(args)))))
  results <- Filter(function(t) inherits(t, "jbrowse_result"), tracks)
  if (!length(results)) {
    stop("jbrowse_stack needs at least one ggjbrowse reader's result",
      call. = FALSE)
  }
  same_regions(results)
  px <- stack_px(results)
  stack_bands(unname(Map(function(track, name)
    titled(band_of(track, results[[1]], px), name), tracks,
    names(tracks) %||% rep("", length(tracks)))))
}

# An argument's tracks, a lone track named by the argument and a list's by
# its own names
named_tracks <- function(tracks, name) {
  listed <- listed_tracks(tracks)
  if (length(listed) == 1 && is.null(names(listed))) names(listed) <- name
  if (is.null(names(listed))) names(listed) <- rep("", length(listed))
  listed
}

# A band titled on its first plot, as the view labels a track above it
titled <- function(band, name) {
  if (!nzchar(name)) return(band)
  band$plots[[1]] <- band$plots[[1]] + ggplot2::labs(title = name) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, hjust = 0,
      margin = ggplot2::margin(0, 0, 2, 0)))
  band
}

# Tracks as a list: one result or ggplot, or a list of them
listed_tracks <- function(tracks) {
  if (inherits(tracks, "jbrowse_result") || ggplot2::is_ggplot(tracks)) {
    list(tracks)
  } else as.list(tracks)
}

# A track's bands on the axis of `axis`, a result: a ggplot of one's own as
# one band, a result as its display's at the stack's `px`
band_of <- function(track, axis, px = NULL) {
  if (ggplot2::is_ggplot(track)) own_band(track, axis)
  else track_bands(track, px)
}

#' The figure height that draws a stack at its display's sizes
#'
#' The inches a [jbrowse_stack()] figure needs to give each band the pixels
#' the JBrowse display gives it, at 150 pixels an inch, the scale
#' [jbrowse_width()] gives the panels across, with room for the x axis and for the keys beside
#' the stack, measured as ggplot2 lays them out: a compact pileup makes a
#' shorter figure than a normal one. Pass it as `ggsave()`'s `height`.
#' @param fig a [jbrowse_stack()] figure, restyled or not
#' @return inches, from 2 to 49
#' @export
#' @examples
#' bam <- system.file("extdata", "gapdh_direct_rna.bam",
#'   package = "ggjbrowse")
#' reads <- jbrowse_alignments(bam, "chr12:6,534,201-6,538,700",
#'   display = list(featureHeight = 3))
#' jbrowse_height(jbrowse_stack(reads))
jbrowse_height <- function(fig) {
  heights <- fig$patches$layout$heights
  if (is.null(heights)) {
    stop("jbrowse_height() measures a jbrowse_stack() figure", call. = FALSE)
  }
  bands <- sum(heights) / PX_PER_INCH + 0.08 * length(heights) + 0.7
  grob <- patchwork::patchworkGrob(fig)
  keys <- vapply(grob$grobs[startsWith(grob$layout$name, "guide-box")],
    function(g) grid::convertHeight(grid::grobHeight(g), "in",
      valueOnly = TRUE), 0)
  min(49, max(2, bands, max(0, keys) + 0.3))
}

#' The figure width that draws a stack at its display's sizes
#'
#' The inches a [jbrowse_stack()] figure needs to give its panels the pixels
#' the tracks were read at, at the 150 pixels an inch [jbrowse_height()]
#' gives each band, with room for the y axis and the keys beside the stack.
#' A track read at a `width` of 2,500 pixels then draws its labels as large
#' as one read at 800, where a fixed figure width would shrink them to fit.
#' Pass it as `ggsave()`'s `width`.
#' @inheritParams jbrowse_height
#' @return inches, at most 49; a stack of tracks read with no `width` of
#'   their own spans 1,500 pixels
#' @export
#' @examples
#' gff <- system.file("extdata", "gapdh_gencode.gff3.gz", package = "ggjbrowse")
#' genes <- jbrowse_features(gff, "chr12:6,534,201-6,538,700", width = 2500)
#' jbrowse_width(jbrowse_stack(genes))
jbrowse_width <- function(fig) {
  if (is.null(fig$patches$layout$heights)) {
    stop("jbrowse_width() measures a jbrowse_stack() figure", call. = FALSE)
  }
  px <- max(vapply(c(fig$patches$plots, list(fig)), function(p)
    result_px(attr(p$data, "jbrowse")), 0))
  if (px == 0) px <- 10 * PX_PER_INCH
  grob <- patchwork::patchworkGrob(fig)
  margins <- grid::convertWidth(sum(grob$widths), "in", valueOnly = TRUE)
  min(49, px / PX_PER_INCH + margins)
}

# The display pixels a result was read across, 0 for none or a track read
# with no width
result_px <- function(x) {
  if (is.null(x)) return(0)
  if (!is.null(x$width)) return(x$width)
  if (is.null(x$bp_per_px)) return(0)
  r <- regions_read(x)
  sum(r$end - r$start) / x$bp_per_px
}

# The widest pixels any of `results` was read across, NULL if none was read
# at a width, so an alignments track stacked with them draws at that width
stack_px <- function(results) {
  px <- max(0, vapply(results, result_px, 0))
  if (px > 0) px
}

# The display pixels a figure gives an inch: a 1,500 pixel view across 10
# inches
PX_PER_INCH <- 150

regions_read <- function(x) x$regions %||% as.data.frame(x$region)

# Regions as JBrowse's location box writes them
loci_of_layout <- function(layout) {
  sprintf("%s:%s-%s", layout$chrom, format(layout$start + 1, big.mark = ",",
    trim = TRUE), format(layout$end, big.mark = ",", trim = TRUE))
}

# The loci a result reads
loci_of <- function(x) loci_of_layout(regions_read(x))

# The loci a result reads, for its print header
print_loci <- function(x) paste(loci_of(x), collapse = ", ")

# Each region's round breaks in its own coordinates, about as many as its
# share of the axis earns, none on an edge facing another region
region_breaks <- function(r) {
  w <- r$end - r$start
  lapply(seq_len(nrow(r)), function(i) {
    b <- pretty(c(r$start[i], r$end[i]), n = max(2, round(6 * w[i] / sum(w))))
    lo <- if (i == 1) r$start[i] else r$start[i] + 1
    hi <- if (i == nrow(r)) r$end[i] else r$end[i] - 1
    b[b >= lo & b <= hi]
  })
}

# Stops where the results read different spans: tracks on one axis must read
# the same regions. A contig may be spelt differently by each file, as an
# exported view spells each track's (`1` beside `chr1`), so a `chr` prefix
# passes, another spelling warns, and contigs that do not pair one to one
# across the tracks stop.
same_regions <- function(results) {
  span <- function(x) {
    r <- regions_read(x)
    paste(r$start, r$end, collapse = ";")
  }
  spans <- vapply(results, span, "")
  differ <- which(spans != spans[1])
  if (length(differ)) {
    stop("jbrowse_stack draws its tracks over the same regions, but track ",
      differ[1], " reads ", paste(loci_of(results[[differ[1]]]),
        collapse = ", "), " and track 1 ", paste(loci_of(results[[1]]),
        collapse = ", "), call. = FALSE)
  }
  bare <- function(chrom) sub("^chr", "", chrom, ignore.case = TRUE)
  first <- regions_read(results[[1]])$chrom
  for (i in seq_along(results)[-1]) {
    chrom <- regions_read(results[[i]])$chrom
    pairs <- unique(data.frame(a = first, b = chrom))
    if (anyDuplicated(pairs$a) || anyDuplicated(pairs$b)) {
      stop("jbrowse_stack draws its tracks over the same contigs, but track ",
        i, " reads ", paste(loci_of(results[[i]]), collapse = ", "),
        " and track 1 ", paste(loci_of(results[[1]]), collapse = ", "),
        call. = FALSE)
    }
    renamed <- pairs[bare(pairs$a) != bare(pairs$b), , drop = FALSE]
    if (nrow(renamed)) {
      warning("track ", i, " reads ", renamed$b[1], " where track 1 reads ",
        renamed$a[1], "; jbrowse_stack draws them as one contig, as it ",
        "would two spellings of one", call. = FALSE)
    }
  }
}

# A ggplot of the caller's own as a band on the stack's x range and region
# labels, a result's table in it moved onto the axis, and its own y limits and
# x scale kept
own_band <- function(p, axis) {
  ylim <- if (inherits(p$coordinates, "CoordCartesian")) {
    p$coordinates$limits$y
  }
  p <- suppressMessages(p + ggplot2::coord_cartesian(xlim = xlim_of(axis),
    ylim = ylim, expand = FALSE))
  list(plots = list(onto_axis(p, axis)), heights = 100)
}

# A plot over a result's regions with every table in it on the axis that lays
# them end to end, the plot's own and each layer's, and that axis's scale
# where the plot sets none
onto_axis <- function(p, axis) {
  if (is.null(axis$regions)) return(p)
  if (is.data.frame(p$data)) p$data <- frame_on_axis(p$data, axis$regions)
  p$layers <- lapply(p$layers, function(layer) {
    if (!is.data.frame(layer$data)) return(layer)
    moved <- ggplot2::ggproto(NULL, layer)
    moved$data <- frame_on_axis(layer$data, axis$regions)
    moved
  })
  if (is.null(p$scales$get_scales("x"))) p <- p + genomic_x(axis)
  p
}

#' @exportS3Method ggplot2::ggplot
ggplot.jbrowse_result <- function(data = NULL, mapping = ggplot2::aes(), ...,
                                  environment = parent.frame()) {
  p <- NextMethod()
  class(p) <- c("ggjbrowse_plot", class(p))
  p
}

#' @exportS3Method ggplot2::ggplot_build
ggplot_build.ggjbrowse_plot <- function(plot, ...) {
  plot <- onto_axis(plot, attr(plot$data, "jbrowse"))
  NextMethod()
}

#' The x axis of a result's regions as a ggplot scale
#'
#' The scale a result's regions draw on: genomic positions written with
#' commas, or, for several regions laid end to end, each tick labelled with
#' its position in its own region and a bar at each seam. `ggplot(<result>)`
#' adds it where the plot sets no x scale of its own, and moves the result's
#' tables in other layers onto it.
#' @param object a ggjbrowse reader's result
#' @param ... passed to [ggplot2::scale_x_continuous()], such as `name`, in
#'   place of the regions' own `breaks` and `labels`
#' @return a list of a ggplot scale and the seams' layer
#' @export
#' @examples
#' gff <- system.file("extdata", "gapdh_gencode.gff3.gz", package = "ggjbrowse")
#' g <- jbrowse_features(gff, c("chr12:6,534,201-6,535,700",
#'   "chr12:6,537,201-6,538,700"))
#' library(ggplot2)
#' ggplot(g) + geom_jbrowse_glyphs() + scale_x_jbrowse(g)
scale_x_jbrowse <- function(object, ...) genomic_x(object, ...)

# A track's bands and their heights in display pixels: a feature track as
# tall as its laid-out pixels, an alignments track's bands as its display
# sizes them at `px`
track_bands <- function(x, px = NULL) {
  if (inherits(x, "jbrowse_alignments")) {
    return(alignment_bands(x, alignment_options(x, width = px)))
  }
  height <- if (inherits(x, "jbrowse_features")) max(30, x$height)
    else if (inherits(x, "jbrowse_variants")) row_px(x) * nrow(x$rows)
    else x$height
  list(plots = list(ggplot2::autoplot(x)), heights = height)
}

# Bands sit close, as the display stacks its tracks and a track its bands.
# The plot at index `axis_at` keeps its x axis: the last by default, none for 0.
stack_bands <- function(bands, axis_at = NULL) {
  plots <- do.call(c, lapply(bands, `[[`, "plots"))
  heights <- unlist(lapply(bands, `[[`, "heights"))
  bare <- seq_along(plots) != (axis_at %||% length(plots))
  plots <- lapply(plots, function(p)
    p + ggplot2::theme(plot.margin = ggplot2::margin(1, 26, 1, 5.5)))
  plots[bare] <- lapply(plots[bare], function(p) p + no_x)
  patchwork::wrap_plots(plots, ncol = 1, heights = heights,
    guides = "collect")
}

# The x axis's title: the contig, the two regions of an axis that lays two end
# to end, or the contigs and the count of more
axis_label <- function(x) {
  r <- x$regions
  if (is.null(r)) return(sprintf("%s (bp)", x$region$chrom))
  if (nrow(r) <= 2) {
    return(sprintf("%s (bp, end to end)", paste(loci_of(x), collapse = ", ")))
  }
  sprintf("%s, %d regions (bp, end to end)",
    paste(unique(r$chrom), collapse = ", "), nrow(r))
}

# The x scale: positions written with commas, and on an axis laying several
# regions end to end round genomic positions in each region, about as many as
# its share of the axis earns, each labelled with its position in its own
# region, and none on an edge two regions face across their gap, or within a
# label's width of the last one kept, where two labels would collide. A bar
# in the view's separator grey marks each seam across the band, over what the
# band draws, as the view masks its tracks.
genomic_x <- function(x, ...) {
  r <- x$regions
  if (is.null(r)) {
    return(scale_with(ggplot2::scale_x_continuous,
      list(labels = scales::label_comma()), ...))
  }
  at <- region_breaks(r)
  breaks <- unlist(lapply(seq_len(nrow(r)), function(i)
    at[[i]] - r$start[i] + r$offset[i]))
  labels <- format(unlist(at), big.mark = ",", trim = TRUE,
    scientific = FALSE)
  room <- 0.09 * (max(r$cum_end) - min(r$cum_start))
  kept <- Reduce(function(k, i)
    if (!length(k) || breaks[i] - breaks[k[length(k)]] >= room) c(k, i)
    else k, seq_along(breaks), integer(0))
  seams <- (utils::head(r$cum_end, -1) + r$offset[-1]) / 2
  list(scale_with(ggplot2::scale_x_continuous, list(breaks = breaks[kept],
      labels = labels[kept]), ...),
    ggplot2::geom_vline(xintercept = seams, colour = REGION_SEPARATOR,
      linewidth = 0.9))
}

# The view's seam between two regions, its palette's regionSeparator
REGION_SEPARATOR <- "#757575"
