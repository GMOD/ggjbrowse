# A multi-genome alignment file's adapter, read off its extension: an
# all-vs-all PAF of PanSN-named sequences (plain or gzipped), or the same
# indexed as a PIF
multiway_adapter <- function(file, index) {
  name <- tolower(sub("\\?.*$", "", file))
  if (grepl("\\.paf(\\.gz)?$", name)) {
    return(list(type = "MultiGenomePAFAdapter", pafLocation = location(file)))
  }
  if (grepl("\\.pif\\.gz$", name)) {
    index <- index %||% paste0(file, ".tbi")
    return(list(type = "MultiGenomeIndexedPAFAdapter",
      pifGzLocation = location(file),
      index = list(location = location(index),
        indexType = if (grepl("\\.csi$", index)) "CSI" else "TBI")))
  }
  stop("cannot tell the format of ", file,
    ": expected an all-vs-all .paf (optionally .gz) or a .pif.gz with its index",
    call. = FALSE)
}

# `lanes` as the display's `rows` and the lane picker's pins: unnamed, the
# lanes drawn in order; named, each value a contig, or a position on it, the
# lane is pinned to ("" pins nothing)
lane_choice <- function(lanes) {
  if (is.null(lanes)) return(list())
  names_of <- names(lanes) %||% as.character(unlist(lanes))
  pins <- if (is.null(names(lanes))) list() else
    Filter(Negate(is.null), stats::setNames(lapply(as.character(lanes),
      lane_pin), names_of))
  list(rows = list(domain = I(names_of), kept = I(names_of)), pins = pins)
}

lane_pin <- function(locus) {
  if (is.na(locus) || !nzchar(locus)) return(NULL)
  if (!grepl(":", locus)) return(list(refName = locus))
  m <- regmatches(locus, regexec("^(.+):([0-9,]+)(?:(?:-|\\.\\.)([0-9,]+))?$",
    locus, perl = TRUE))[[1]]
  if (!length(m)) {
    stop("cannot read \"", locus, "\" as a lane's contig or position: write ",
      "\"Chr2\", \"Chr2:1,000,000\" or \"Chr2:1,000,001-2,000,000\"",
      call. = FALSE)
  }
  number <- function(x) as.numeric(gsub(",", "", x))
  from <- number(m[3]) - 1
  to <- if (nzchar(m[4])) number(m[4]) else from + 1
  list(refName = m[2], at = (from + to) / 2)
}

#' Read a multi-genome alignment as JBrowse's multi-way synteny display stacks
#' it
#'
#' Runs JBrowse's own multi-way synteny display in V8 under a linear genome
#' view of `locus`: the anchor lane on top, then one lane per other genome,
#' each in its own coordinates. Each lane picks its contig, its block, its
#' zoom and its orientation from the alignments as the display does, and a
#' lane opens a hole where its sequence lacks a stretch of the anchor. Ribbons
#' join each alignment's placements on adjacent lanes. Coordinates are
#' 0-based, half-open, on each lane's assembly; `x1` and `x2` columns are the
#' view's pixels, where each lane's frame places them.
#'
#' @param file an all-vs-all PAF (plain or gzipped) whose sequences are named
#'   `sample#haplotype#contig` (PanSN), or the same as an indexed `.pif.gz`,
#'   a path or a URL; or a JBrowse multi-genome synteny track's `adapter`
#'   config, as a list, such as a `MultiGenomePAFAdapter` declaring its
#'   `attributeColumns`, an `MCScanBlocksAdapter` or a
#'   `MultiPairwiseSyntenyAdapter`
#' @param locus the anchor lane's regions: a locus as JBrowse's location box
#'   writes it (`"ctgA:1-50,001"`, 1-based and inclusive), several of them,
#'   drawn end to end, or a data frame of `chrom`, `start` and `end`
#' @param assembly the anchor's assembly, the PanSN sample name for a PAF;
#'   for an adapter config naming its assemblies, unset takes its first
#' @param lanes the other lanes: unset, every genome the alignments place,
#'   densest first; a character vector of assemblies, drawn in that order;
#'   or a named vector of assemblies whose values pin a lane, as the lane
#'   picker does, to a contig (`"Chr2"`) or a position on one
#'   (`"Chr2:1,000,000"`), `""` pinning nothing
#' @param display the display's settings, as a `MultiWaySyntenyDisplay`
#'   config writes them: a list with any of `ribbonColor` (a CSS color, or
#'   the color object's slots, such as `list(field = "strand")`,
#'   `list(field = "identity")` or a declared column), `hideUnlabelled`,
#'   `drawCurves`, `bridgeSkippedLanes` and `rows`, which `lanes` writes
#' @param width the view's pixels, which set each lane's zoom and which
#'   ribbons are wide enough to draw
#' @param index the tabix or CSI index of a `.pif.gz`, if not beside it
#' @return a `jbrowse_multiway` list:
#'   * `lanes`: one row per lane, top to bottom: `lane`, `assembly`,
#'     `label`, `anchor`, `chrom`, `flipped` (the lane runs right to left),
#'     `bp_per_px`, `also_on` (other contigs placing a comparable share)
#'     and `locus`, the lane's window, to read its tracks over
#'   * `segments`: the stretches of each lane's sequence across the view,
#'     cut where it opens a hole: `lane`, `chrom`, `start`, `end`, and `x1`,
#'     `x2`, the pixels `start` and `end` sit at
#'   * `blocks`: each alignment or gene's placement on each lane: `lane`,
#'     `group`, `name`, `chrom`, `start`, `end`, `strand` against the anchor,
#'     `x1`, `x2`
#'   * `ribbons`: what each gutter draws, in paint order: `alignment` (a row
#'     of `alignments`), `group`, `kind` (`"match"`, or `"mismatch"` for a
#'     base tile), `gutter` (the lane it hangs under), `upper_lane`,
#'     `lower_lane`, each foot's `*_chrom`, `*_start` and `*_end`, and its
#'     pixels `upper_x1`, `upper_x2`, `lower_x1`, `lower_x2`, where the ribbon
#'     joins `upper_x1` to `lower_x1`, so a reversed pair crosses; `color`,
#'     `alpha`, `thin` and `bridge`, a ribbon that skips a lane placing
#'     nothing for its group
#'   * `alignments`: the records the ribbons draw: `id`, `group`, `name`,
#'     `assembly`, `chrom`, `start`, `end`, `strand`, `mate_assembly`,
#'     `mate_chrom`, `mate_start`, `mate_end`, and the measurements and
#'     columns the file carries
#'   * `legend`: `label` and `color`, the ribbons' key
#'
#'   `width` is the view's pixels, `tick_interval` the bp between each lane's
#'   ticks, and `rows_vs_anchor` whether each lane is stated against the
#'   anchor alone, as a star of pairwise files states it.
#' @export
#' @examples
#' paf <- system.file("extdata", "volvox_all_vs_all.paf", package = "ggjbrowse")
#' m <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox")
#' m$lanes
#' head(m$ribbons)
jbrowse_multiway <- function(file, locus, assembly = NULL, lanes = NULL,
                             display = NULL, width = 800, index = NULL) {
  if (is.null(assembly) && !is.list(file)) {
    stop("name the anchor lane's assembly, the PanSN sample its sequences ",
      "are named after", call. = FALSE)
  }
  regions <- regions_of(locus, NULL, NULL)
  choice <- lane_choice(lanes)
  display <- display %||% list()
  if (!is.null(choice$rows)) display$rows <- choice$rows
  display$ribbonColor <- config_color(display$ribbonColor)
  args <- list(adapterConfig = adapter_of(file, function()
    multiway_adapter(file, index)),
    anchor = axis_args(regions, assembly, width), display = display)
  if (length(choice$pins)) args$pins <- choice$pins
  r <- call_engine("multiway", args)
  one_based <- function(d, columns) {
    for (column in intersect(columns, names(d))) d[[column]] <- d[[column]] + 1
    d
  }
  lane_frame <- one_based(snake(as_frame(r$lanes)), "lane")
  segments <- one_based(snake(as_frame(r$segments)), "lane")
  lane_frame$locus <- lane_loci(lane_frame, segments)
  ribbons <- one_based(snake(as_frame(r$ribbons)),
    c("alignment", "gutter", "upper_lane", "lower_lane"))
  ribbons$color <- jbrowse_color(ribbons$color)
  key <- key_frame(unlist(r$legend$label), unlist(r$legend$color),
    r$legendTitle)
  attr(key, "override") <- list(alpha = 0.6)
  attr(key, "complete") <- TRUE
  structure(
    list(lanes = lane_frame, segments = segments,
      blocks = one_based(snake(as_frame(r$blocks)), "lane"),
      ribbons = ribbons, alignments = snake(as_frame(r$alignments)),
      legend = key, keys = list(ribbons = key), ramps = ramp_frames(r$ramps),
      field = r$field, tick_interval = r$tickInterval,
      rows_vs_anchor = r$rowsVsAnchor, curves = r$curves, width = width,
      assemblies = unlist(r$assemblyNames), anchor = region_layout(regions)),
    primary = "ribbons",
    class = c("jbrowse_multiway", "jbrowse_result"))
}

# Each lane's window as a locus: the anchor's regions, or the stretch of its
# contig a mate lane shows across the view
lane_loci <- function(lanes, segments) {
  vapply(lanes$lane, function(i) {
    s <- segments[segments$lane == i, , drop = FALSE]
    if (!nrow(s)) return(NA_character_)
    if (lanes$anchor[i]) return(paste(loci_of_layout(s), collapse = " "))
    loci_of_layout(data.frame(chrom = s$chrom[1], start = min(s$start),
      end = max(s$end)))
  }, "")
}

#' @export
print.jbrowse_multiway <- function(x, ...) {
  cat(sprintf("<jbrowse_multiway> %d lanes under %s, %d ribbons\n",
    nrow(x$lanes), paste(loci_of_layout(x$anchor), collapse = ", "),
    nrow(x$ribbons)))
  for (i in seq_len(nrow(x$lanes))) {
    cat(sprintf("  %d %s %s%s\n", i, x$lanes$assembly[i],
      x$lanes$locus[i] %|NA|% "places nothing in view",
      if (x$lanes$flipped[i]) ", reversed" else ""))
  }
  invisible(x)
}

`%|NA|%` <- function(a, b) if (is.na(a)) b else a

# The stack's heights in the figure's pixels: a gutter between lanes, and a
# lane's ruler labels above its line
MULTIWAY_GUTTER_PX <- 100
LANE_RULER_PX <- 22
LANE_BOX_PX <- 5

# The stack laid out top to bottom, a lane's line at y = its pixels above
# the bottom lane's: each ribbon from under its upper lane to over its lower
# lane's ruler, each lane's ticks at the shared interval
multiway_on_axis <- function(x) {
  if (!is.null(x$ribbons$.on_axis)) return(x)
  n <- nrow(x$lanes)
  lane_y <- (n - seq_len(n)) * MULTIWAY_GUTTER_PX
  r <- x$ribbons
  r$y_top <- lane_y[r$gutter] - LANE_BOX_PX
  r$y_bottom <- lane_y[r$lower_lane] + LANE_RULER_PX
  r$.on_axis <- rep(TRUE, nrow(r))
  x$ribbons <- r
  x$segments$y <- lane_y[x$segments$lane]
  x$blocks$y <- lane_y[x$blocks$lane]
  x$ticks <- lane_ticks(x$segments, x$tick_interval)
  x$lane_y <- lane_y
  x
}

# Every multiple of `interval` on each segment, at its pixel on the lane
lane_ticks <- function(segments, interval) {
  ticks <- lapply(seq_len(nrow(segments)), function(i) {
    s <- segments[i, ]
    bp <- seq(ceiling(s$start / interval), floor(s$end / interval)) * interval
    bp <- bp[bp >= s$start & bp <= s$end]
    if (!length(bp) || s$end == s$start) return(NULL)
    data.frame(lane = s$lane, bp = bp, y = s$y,
      x = s$x1 + (bp - s$start) / (s$end - s$start) * (s$x2 - s$x1))
  })
  d <- do.call(rbind, ticks)
  if (is.null(d)) {
    return(data.frame(lane = numeric(0), bp = numeric(0), y = numeric(0),
      x = numeric(0)))
  }
  d[!duplicated(d[c("lane", "bp")]), , drop = FALSE]
}

multiway_feet <- function(r) {
  list(top1 = r$upper_x1, top2 = r$upper_x2, bottom1 = r$lower_x1,
    bottom2 = r$lower_x2, y_top = r$y_top, y_bottom = r$y_bottom)
}

#' A multi-way stack's ribbons and lanes as ggplot layers
#'
#' `geom_jbrowse_lane_ribbons()` draws each ribbon from under its upper lane
#' to over its lower lane as the display draws it: a filled trapezoid in its
#' color at its opacity, or a one-pixel line where it is thinner than a
#' pixel. `geom_jbrowse_lanes()` draws each lane's line, its ruler and, for a
#' gene table, the boxes the ribbons join. Lanes stack top to bottom
#' `r MULTIWAY_GUTTER_PX` pixels apart on y; x is the view's pixels.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_multiway()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param curves draw each ribbon edge as the display's `drawCurves` does,
#'   unset as the result's settings say
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' paf <- system.file("extdata", "volvox_all_vs_all.paf", package = "ggjbrowse")
#' m <- jbrowse_multiway(paf, "ctgA:1-50,001", "volvox")
#' library(ggplot2)
#' ggplot(m) + geom_jbrowse_lane_ribbons() + geom_jbrowse_lanes()
geom_jbrowse_lane_ribbons <- function(mapping = NULL, data = NULL,
                                      curves = NULL, ...) {
  data <- replacing(data, "ribbons", "geom_jbrowse_lane_ribbons")
  curved <- function(x) curves %||% isTRUE(x$curves)
  filled <- function(x) {
    r <- x$ribbons[!x$ribbons$thin, , drop = FALSE]
    ribbon_polygons(r, curved(x), multiway_feet(r))
  }
  lines <- function(x) {
    r <- x$ribbons[x$ribbons$thin, , drop = FALSE]
    ribbon_lines(r, curved(x), multiway_feet(r))
  }
  list(
    jbrowse_layer(ggplot2::GeomPolygon, data, filled,
      for_geom(mapping, ggplot2::GeomPolygon),
      ggplot2::aes(x = x, y = y, group = ribbon, fill = color,
        alpha = I(alpha)), list(colour = NA, ...), keyed = "complete"),
    jbrowse_layer(ggplot2::GeomPath, data, lines,
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = y, group = ribbon, colour = color,
        alpha = I(alpha)), list(linewidth = 0.3, ...), keyed = FALSE))
}

#' @rdname geom_jbrowse_lane_ribbons
#' @export
geom_jbrowse_lanes <- function(data = NULL, ...) {
  data <- replacing(data, NULL, "geom_jbrowse_lanes")
  named <- function(x) x$blocks[!is.na(x$blocks$name), , drop = FALSE]
  tick_px <- 4
  list(
    jbrowse_layer(ggplot2::GeomRect, data, named, NULL,
      ggplot2::aes(xmin = pmin(x1, x2), xmax = pmax(x1, x2),
        ymin = y - LANE_BOX_PX, ymax = y + LANE_BOX_PX),
      list(fill = "grey60", colour = "grey30", linewidth = 0.2, ...),
      keyed = FALSE),
    jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$segments, NULL,
      ggplot2::aes(x = x1, xend = x2, y = y, yend = y),
      list(colour = "grey35", linewidth = 0.4, ...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$ticks, NULL,
      ggplot2::aes(x = x, xend = x, y = y, yend = y + tick_px),
      list(colour = "grey35", linewidth = 0.3, ...), keyed = FALSE),
    jbrowse_layer(ggplot2::GeomText, data, function(x) x$ticks, NULL,
      ggplot2::aes(x = x, y = y + tick_px + 2,
        label = scales::label_number(scale_cut = scales::cut_si("bp"))(bp)),
      list(size = 2.2, vjust = 0, colour = "grey30", ...), keyed = FALSE))
}

# Each lane's name on the y axis: its assembly, and the contig it shows
lane_names <- function(x) {
  l <- x$lanes
  where <- ifelse(is.na(l$chrom), "nothing in view",
    paste0(l$chrom, ifelse(l$flipped, ", reversed", "")))
  paste0(l$label, "\n", where)
}

multiway_stack <- function(x, curves) {
  x <- on_axis(x)
  p <- ggplot2::ggplot(x) +
    geom_jbrowse_lane_ribbons(curves = curves) +
    geom_jbrowse_lanes() +
    scale_fill_jbrowse(x) +
    ggplot2::scale_y_continuous(breaks = x$lane_y, labels = lane_names(x)) +
    ggplot2::coord_cartesian(xlim = c(0, x$width),
      ylim = c(-LANE_BOX_PX, max(x$lane_y) + LANE_RULER_PX), expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 1, lineheight = 0.9))
  p
}

# A track over a lane's window drawn on the stack's pixels: each genomic
# position moved to its pixel on the lane, through the lane's segments
lane_band <- function(track, x, lane, where) {
  s <- x$segments[x$segments$lane == lane, , drop = FALSE]
  read <- regions_read(track)
  if (!nrow(s) || !all(read$chrom %in% s$chrom) ||
    min(read$start) > min(s$start) || max(read$end) < max(s$end)) {
    stop("a track ", where, " the stack draws over ", x$lanes$locus[lane],
      ", but one reads ", paste(loci_of(track), collapse = ", "),
      call. = FALSE)
  }
  bp <- c(rbind(s$start, s$end))
  px <- c(rbind(s$x1, s$x2))
  if (lane == 1 && nrow(x$anchor) > 1) {
    a <- x$anchor
    bp <- c(rbind(a$cum_start, a$cum_end))
  }
  if (is.unsorted(bp)) {
    stop("a track ", where, " the stack draws over one region only",
      call. = FALSE)
  }
  to_px <- function(v) stats::approx(bp, px, v, rule = 2, ties = "ordered")$y
  from_px <- function(v) stats::approx(px, bp, v, rule = 2,
    ties = "ordered")$y
  band <- track_bands(track)
  band$plots <- lapply(band$plots, function(p) suppressMessages(p +
    ggplot2::scale_x_continuous(transform = scales::new_transform("lane_px",
      to_px, from_px)) +
    ggplot2::coord_cartesian(xlim = from_px(c(0, x$width)), expand = FALSE)))
  band
}

#' Draw a multi-genome alignment as JBrowse's multi-way synteny display
#'
#' The lanes top to bottom, each with its own ruler and its genome named on
#' the left, the ribbons between them, and any tracks over the anchor lane
#' above the stack and over the bottom lane below it: a
#' [jbrowse_features()], [jbrowse_wiggle()], [jbrowse_variants()] or
#' [jbrowse_alignments()] result read over that lane's window, as
#' `x$lanes$locus` writes it.
#'
#' @param object a [jbrowse_multiway()] result
#' @param above,below tracks over the anchor lane and over the bottom lane, a
#'   result or a list of them
#' @param curves draw each ribbon edge as a curve, as the display's
#'   `drawCurves` does; unset, as the result's settings say
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot of the stack alone, or with tracks a patchwork, whose
#'   height [jbrowse_height()] measures
#' @exportS3Method ggplot2::autoplot
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' m <- jbrowse_multiway(e("volvox_all_vs_all.paf"), "ctgA:1-50,001", "volvox",
#'   display = list(ribbonColor = list(field = "strand")))
#' autoplot(m)
#' genes <- jbrowse_features(e("volvox.sort.gff3.gz"), m$lanes$locus[1],
#'   display = list(showOnlyGenes = TRUE))
#' autoplot(m, above = genes)
autoplot.jbrowse_multiway <- function(object, above = NULL, below = NULL,
                                      curves = NULL, ...) {
  rlang::check_dots_empty()
  stack <- multiway_stack(object, curves)
  top <- lapply(listed_tracks(above), lane_band, object, 1, "above")
  bottom <- lapply(listed_tracks(below), lane_band, object, nrow(object$lanes),
    "below")
  if (!length(top) && !length(bottom)) return(stack)
  n <- nrow(object$lanes)
  bands <- c(top,
    list(list(plots = list(stack),
      heights = (n - 1) * MULTIWAY_GUTTER_PX + LANE_RULER_PX + LANE_BOX_PX)),
    bottom)
  plots <- do.call(c, lapply(bands, `[[`, "plots"))
  heights <- unlist(lapply(bands, `[[`, "heights"))
  plots <- lapply(plots, function(p) p + no_x +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 14, 1, 5.5)))
  patchwork::wrap_plots(plots, ncol = 1, heights = heights,
    guides = "collect")
}
