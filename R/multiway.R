# A multi-genome alignment file's adapter, read off its extension: an
# all-vs-all PAF of PanSN-named sequences (plain or gzipped), or the same
# indexed as a PIF
multiway_adapter <- function(file, index) {
  guess_adapter(file, index, list(
    list(pattern = "\\.paf(\\.gz)?$", type = "MultiGenomePAFAdapter",
      slot = "pafLocation"),
    list(pattern = "\\.pif\\.gz$", type = "MultiGenomeIndexedPAFAdapter",
      slot = "pifGzLocation", index = "TBI")),
    "an all-vs-all .paf (optionally .gz) or a .pif.gz with its index")
}

# `lanes` as the display's `rows` and the lane picker's pins, records of
# `lane`, `refName` and `at`: unnamed, the lanes drawn in order; named, each
# value a contig, or a position on it, the lane is pinned to ("" pins
# nothing). An unnamed entry among named ones is a lane, pinned nowhere.
lane_choice <- function(lanes) {
  if (is.null(lanes)) return(list())
  values <- vapply(as.list(lanes), function(v)
    if (length(v)) as.character(v)[1] else "", "")
  named <- nzchar(names(lanes) %||% rep("", length(values)))
  names_of <- ifelse(named, names(lanes) %||% "", values)
  pins <- Map(function(lane, locus) {
    pin <- lane_pin(locus)
    if (!is.null(pin)) c(list(lane = lane), pin)
  }, names_of, ifelse(named, values, ""), USE.NAMES = FALSE)
  list(rows = list(sort = names_of, kept = names_of),
    pins = Filter(Negate(is.null), pins))
}

lane_pin <- function(locus) {
  if (is.na(locus) || !nzchar(locus)) return(NULL)
  if (!grepl(":", locus)) return(list(refName = locus))
  m <- regmatches(locus, regexec("^(.+):([0-9,]+)(?:(?:-|\\.\\.)([0-9,]+))?$",
    locus, perl = TRUE))[[1]]
  if (!length(m)) {
    rlang::abort(paste0("cannot read \"", locus, "\" as a lane's contig or ",
      "position: write \"Chr2\", \"Chr2:1,000,000\" or ",
      "\"Chr2:1,000,001-2,000,000\""), class = "ggjbrowse_bad_region")
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
#' view of `region`: the anchor lane on top, then one lane per other genome,
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
#' @param region the anchor lane's regions, drawn end to end: loci as
#'   JBrowse's location box writes them (`"ctgA:1-50,001"`, 1-based and
#'   inclusive), a data frame of `chrom`, `start` and `end` (0-based,
#'   half-open), or a GRanges. A multi-genome file states no contig's
#'   length, so a contig's name alone or `NULL`, which a linear reader reads
#'   whole, raises `ggjbrowse_bad_region`
#' @param anchor the anchor lane's assembly name, the PanSN sample name for a
#'   PAF; for an adapter config naming its assemblies, unset takes its first
#' @param lanes the other lanes: unset, every genome the alignments place,
#'   densest first; a character vector of assemblies, drawn in that order;
#'   or a named vector of assemblies whose values pin a lane, as the lane
#'   picker does, to a contig (`"Chr2"`) or a position on one
#'   (`"Chr2:1,000,000"`), `""` pinning nothing; an unnamed entry among
#'   named ones is a lane pinned nowhere
#' @param display the display's settings, as a `MultiWaySyntenyDisplay`
#'   config writes them: a list with any of `ribbons` (its `encoding$color`,
#'   `list(value = )` for a CSS color, or a field such as
#'   `list(field = "strand")`, `list(field = "identity")` or a declared
#'   column, its `encoding$opacity`, its `filter`, a list of predicates
#'   such as `list(field = "length", gte = 10000)`, and its
#'   `mark$interpolate`, `"bump"` for curved sides), and `rows`, which
#'   `lanes` writes
#' @param width the view's pixels, 1,500 as for the linear readers, which set
#'   each lane's zoom and which ribbons are wide enough to draw
#' @param index the tabix or CSI index of a `.pif.gz`, if not beside it
#' @return a `jbrowse_multiway` list, beside the slots [jbrowse_result] lists:
#'   * `lanes`: one row per lane, top to bottom: `lane`, `assembly`,
#'     `label`, `anchor`, `chrom`, `flipped` (the lane runs right to left),
#'     `bp_per_px`, `also_on` (other contigs placing a comparable share),
#'     `y_top` and `y_bottom`, the band its gene boxes fill, and `locus`, the
#'     lane's window, to read its tracks over
#'   * `segments`: the stretches of each lane's sequence across the view,
#'     cut where it opens a hole: `lane`, `chrom`, `start`, `end`, and `x1`,
#'     `x2`, the pixels `start` and `end` sit at, and `region`, the row of
#'     `regions` it lies in
#'   * `blocks`: each alignment or gene's placement on each lane: `lane`,
#'     `group`, `name`, `chrom`, `start`, `end`, `strand` against the anchor,
#'     `x1`, `x2`, `region`
#'   * `ribbons`: what each gutter draws, in paint order: `alignment` (the
#'     record's key in `alignments`), `group`, `kind` (`"match"`, or `"mismatch"` for a
#'     base tile), `gutter` (the lane it hangs under), `upper_lane`,
#'     `lower_lane`, each foot's `*_region`, `*_chrom`, `*_start` and
#'     `*_end`, and its
#'     pixels `upper_x1`, `upper_x2`, `lower_x1`, `lower_x2`, where the ribbon
#'     joins `upper_x1` to `lower_x1`, so a reversed pair crosses; `y_top`
#'     and `y_bottom`, its gutter; `color`, `alpha`, `thin` and `bridge`, a
#'     ribbon that skips a lane placing nothing for its group
#'   * `alignments`: the records the ribbons draw: `alignment`, the key
#'     assigning `alignments` drops the ribbons of, `id`, `group`, `name`,
#'     `assembly`, `chrom`, `start`, `end`, `strand`, `mate_assembly`,
#'     `mate_chrom`, `mate_start`, `mate_end`, and the measurements and
#'     columns the file carries
#'
#'   `keys$ribbons` is the ribbons' key. `regions` lays out each lane, its `axis` the lane's assembly: the
#'   anchor's regions, spaced as a linear track over them spaces them, and
#'   the window each other lane shows. `display` is the display's config as
#'   the engine resolves it, `rows` as `lanes` writes it. `width` is the view's pixels, `height` the stack's, each `y_top`
#'   and `y_bottom` pixels down from its top, `tick_interval` the bp between
#'   each lane's ticks, and `rows_vs_anchor` whether each lane is stated
#'   against the anchor alone, as a star of pairwise files states it.
#' @export
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' liftover <- function(genome, assembly) list(type = "PairwiseIndexedPAFAdapter",
#'   pifGzLocation = list(uri = e(sprintf("hg38To%s.over.pif.gz", genome))),
#'   index = list(location = list(uri = e(sprintf("hg38To%s.over.pif.gz.csi",
#'     genome))), indexType = "CSI"), assemblyNames = c(assembly, "hg38"))
#' primates <- list(type = "MultiPairwiseSyntenyAdapter", adapters = list(
#'   liftover("PanTro6", "panTro6"), liftover("GorGor6", "gorGor6"),
#'   liftover("PonAbe3", "ponAbe3")))
#' m <- jbrowse_multiway(primates, "chr17:7,400,001-7,700,000", "hg38",
#'   lanes = c("panTro6", "gorGor6", "ponAbe3"))
#' m$lanes
#' head(m$ribbons)
jbrowse_multiway <- function(file, region, anchor = NULL, lanes = NULL,
                             display = NULL, width = DEFAULT_WIDTH,
                             index = NULL) {
  if (is.null(anchor) && !is.list(file)) {
    stop("name the anchor lane's assembly, the PanSN sample its sequences ",
      "are named after", call. = FALSE)
  }
  regions <- anchor_regions(region)
  choice <- lane_choice(lanes)
  config <- display %||% list()
  config$rows <- choice$rows %||% config$rows
  args <- list(adapterConfig = adapter_of(file, function()
    multiway_adapter(file, index)),
    anchor = axis_args(regions, anchor, width),
    display = display_config(config) %||% list())
  if (length(choice$pins)) args$pins <- choice$pins
  read <- warned_call("multiway", args)
  r <- read$value
  missing <- check_missing(as.list(unlist(r$missingContigs)), regions$chrom)
  regions$missing <- as.vector(missing)
  one_based <- function(d, columns) {
    for (column in intersect(columns, names(d))) d[[column]] <- d[[column]] + 1
    d
  }
  lane_frame <- one_based(snake_frame(r$lanes), "lane")
  segments <- one_based(snake_frame(r$segments), "lane")
  layout <- lane_regions(lane_frame, segments, regions)
  lane_frame$locus <- lane_loci(lane_frame, layout)
  axis <- function(lane) lane_frame$assembly[lane]
  mid <- function(d, prefix = "") (d[[paste0(prefix, "start")]] +
    d[[paste0(prefix, "end")]]) / 2
  segments$region <- region_at(layout, axis(segments$lane), segments$chrom,
    mid(segments))
  blocks <- one_based(snake_frame(r$blocks), "lane")
  blocks$region <- region_at(layout, axis(blocks$lane), blocks$chrom,
    mid(blocks))
  ribbons <- one_based(snake_frame(r$ribbons),
    c("alignment", "gutter", "upper_lane", "lower_lane"))
  ribbons$alignment <- as.integer(ribbons$alignment)
  ribbons$upper_region <- region_at(layout, axis(ribbons$upper_lane),
    ribbons$upper_chrom, mid(ribbons, "upper_"))
  ribbons$lower_region <- region_at(layout, axis(ribbons$lower_lane),
    ribbons$lower_chrom, mid(ribbons, "lower_"))
  comparative_result("jbrowse_multiway", "ribbons",
    list(lanes = lane_frame, segments = segments, blocks = blocks,
      ribbons = ribbons, alignments = keyed_alignments(r$alignments)),
    layout, width, r$config, list(ribbons = legend_key(r, 0.6)),
    ramp_frames(r$ramps),
    unique(c(attr(missing, "warning"), read$warnings)),
    list(height = r$height, tick_interval = r$tickInterval,
      rows_vs_anchor = r$rowsVsAnchor))
}

# The anchor lane's regions: loci, a data frame or a GRanges. A
# multi-genome file states no contig's length, so a whole contig, or the
# whole genome, has no extent to lay out.
anchor_regions <- function(region) {
  if (is.null(region) || (is.character(region) && !all(grepl(LOCUS, region)))) {
    rlang::abort(paste("the anchor lane reads loci such as",
      "\"ctgA:1-50,001\", a data frame of chrom, start and end, or a GRanges:",
      "a multi-genome file states no contig's length, so a whole contig has",
      "no extent"), class = "ggjbrowse_bad_region")
  }
  regions_of(region)
}

# Each lane's regions, its `axis` its assembly: the anchor's regions, laid
# out as a linear track over them lays them, and for each other lane the
# stretch of its contig it shows across the view. A lane placing nothing has
# none.
lane_regions <- function(lanes, segments, anchor) {
  do.call(rbind, lapply(lanes$lane, function(i) {
    name <- lanes$assembly[i]
    if (lanes$anchor[i]) return(axis_rows(anchor, name, name, 0.012))
    s <- segments[segments$lane == i, , drop = FALSE]
    if (!nrow(s)) return(NULL)
    axis_rows(data.frame(chrom = s$chrom[1], start = min(s$start),
      end = max(s$end)), name, name)
  }))
}

# The row of `regions` on each `axis` and contig holding `at`, or its first
# on that contig where none holds it
region_at <- function(regions, axis, chrom, at) {
  vapply(seq_along(axis), function(i) {
    on <- which(regions$axis == axis[i] & regions$chrom == chrom[i])
    inside <- on[at[i] >= regions$start[on] & at[i] <= regions$end[on]]
    c(inside, on, NA_integer_)[1]
  }, 1L)
}

# Each lane's regions as loci, NA for a lane placing nothing
lane_loci <- function(lanes, regions) {
  vapply(lanes$assembly, function(name) {
    r <- on_axis_of(regions, name)
    if (!nrow(r)) NA_character_ else paste(loci_of_layout(r), collapse = " ")
  }, "", USE.NAMES = FALSE)
}

#' @export
print.jbrowse_multiway <- function(x, ...) {
  cat(sprintf("<jbrowse_multiway> %d lanes under %s, %d ribbons\n",
    nrow(x$lanes), x$lanes$locus[1],
    nrow(x$ribbons)))
  for (i in seq_len(nrow(x$lanes))) {
    cat(sprintf("  %d %s %s%s\n", i, x$lanes$assembly[i],
      x$lanes$locus[i] %|NA|% "places nothing in view",
      if (x$lanes$flipped[i]) ", reversed" else ""))
  }
  invisible(x)
}

`%|NA|%` <- function(a, b) if (is.na(a)) b else a

# The stack as the display lays it out, y the px up from its bottom: each
# lane's line through the middle of its gene boxes, each ribbon across its
# gutter, each lane's ticks at the shared interval
multiway_on_axis <- function(x) {
  if (is_on_axis(x$ribbons)) return(x)
  up <- function(y) x$height - y
  l <- x$lanes
  l$box_top <- up(l$y_top)
  l$box_bottom <- up(l$y_bottom)
  x$lanes <- l
  lane_y <- (l$box_top + l$box_bottom) / 2
  r <- x$ribbons
  r$y_top <- up(r$y_top)
  r$y_bottom <- up(r$y_bottom)
  x$ribbons <- structure(r, on_axis = TRUE)
  x$segments$y <- lane_y[x$segments$lane]
  x$blocks$y <- lane_y[x$blocks$lane]
  x$blocks$box_top <- l$box_top[x$blocks$lane]
  x$blocks$box_bottom <- l$box_bottom[x$blocks$lane]
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
#' gene table, the boxes the ribbons join. Lanes stack top to bottom as the
#' display lays them out over its `height`, y the pixels up from the stack's
#' bottom; x is the view's pixels.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_multiway()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param curves draw each ribbon edge as the display's
#'   `ribbons = list(mark = list(interpolate = "bump"))` does, unset as the
#'   result's settings say
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' liftover <- function(genome, assembly) list(type = "PairwiseIndexedPAFAdapter",
#'   pifGzLocation = list(uri = e(sprintf("hg38To%s.over.pif.gz", genome))),
#'   index = list(location = list(uri = e(sprintf("hg38To%s.over.pif.gz.csi",
#'     genome))), indexType = "CSI"), assemblyNames = c(assembly, "hg38"))
#' primates <- list(type = "MultiPairwiseSyntenyAdapter", adapters = list(
#'   liftover("PanTro6", "panTro6"), liftover("GorGor6", "gorGor6"),
#'   liftover("PonAbe3", "ponAbe3")))
#' m <- jbrowse_multiway(primates, "chr17:7,400,001-7,700,000", "hg38",
#'   lanes = c("panTro6", "gorGor6", "ponAbe3"))
#' library(ggplot2)
#' ggplot(m) + geom_jbrowse_lane_ribbons() + geom_jbrowse_lanes()
geom_jbrowse_lane_ribbons <- function(mapping = NULL, data = NULL,
                                      curves = NULL, ...) {
  ribbon_layers(mapping,
    replacing(data, "ribbons", "geom_jbrowse_lane_ribbons"),
    function(x) curves %||% bump(x$display$ribbons$mark), multiway_feet, ...)
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
        ymin = box_bottom, ymax = box_top),
      list(fill = "grey60", colour = "grey30", linewidth = 0.2, ...)),
    jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$segments, NULL,
      ggplot2::aes(x = x1, xend = x2, y = y, yend = y),
      list(colour = "grey35", linewidth = 0.4, ...)),
    jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$ticks, NULL,
      ggplot2::aes(x = x, xend = x, y = y, yend = y + tick_px),
      list(colour = "grey35", linewidth = 0.3, ...)),
    jbrowse_layer(ggplot2::GeomText, data, function(x) x$ticks, NULL,
      ggplot2::aes(x = x, y = y + tick_px + 2,
        label = scales::label_number(scale_cut = scales::cut_si("bp"))(bp)),
      list(size = 2.2, vjust = 0, colour = "grey30", ...)))
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
    guides_jbrowse(x) +
    ggplot2::scale_y_continuous(breaks = x$lane_y, labels = lane_names(x)) +
    ggplot2::coord_cartesian(xlim = c(0, x$width),
      ylim = c(0, x$height), expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 1, lineheight = 0.9))
  p
}

# A track over a lane's window drawn on the stack's pixels: each genomic
# position moved to its pixel on the lane, through the lane's segments, and
# past the lane's ends at its bp per px, where the lane leaves the view
# empty
lane_band <- function(track, x, lane, where) {
  s <- x$segments[x$segments$lane == lane, , drop = FALSE]
  read <- track$regions
  a <- on_axis_of(x$regions, x$lanes$assembly[lane])
  anchor <- lane == 1 && nrow(a) > 1
  if (!nrow(s) || !all(read$chrom %in% s$chrom) ||
    min(read$start) > min(s$start) || max(read$end) < max(s$end) ||
    (anchor && !identical(paste(read$start, read$end),
      paste(a$start, a$end)))) {
    stop("a track ", where, " the stack draws over ", x$lanes$locus[lane],
      ", but one reads ", paste(loci_of(track), collapse = ", "),
      call. = FALSE)
  }
  bp <- c(rbind(s$start, s$end))
  px <- c(rbind(s$x1, s$x2))
  if (anchor) bp <- c(rbind(a$cum_start, a$cum_end))
  if (is.unsorted(bp)) {
    stop("a track ", where, " the stack draws over one region only",
      call. = FALSE)
  }
  if (is.unsorted(px)) {
    stop("a track ", where, " the stack draws over a lane laid out left to ",
      "right, but ", x$lanes$locus[lane], " runs flipped", call. = FALSE)
  }
  per_px <- x$lanes$bp_per_px[lane]
  bp <- c(bp[1] - px[1] * per_px, bp,
    bp[length(bp)] + (x$width - px[length(px)]) * per_px)
  px <- c(0, px, x$width)
  to_px <- function(v) stats::approx(bp, px, v, rule = 2, ties = "ordered")$y
  from_px <- function(v) stats::approx(px, bp, v, rule = 2,
    ties = "ordered")$y
  band <- track_bands(track)
  band$plots <- lapply(band$plots, function(p) suppressMessages(p +
    ggplot2::scale_x_continuous(transform = scales::new_transform("lane_px",
      to_px, from_px)) +
    ggplot2::coord_cartesian(xlim = from_px(c(0, x$width)),
      ylim = p$coordinates$limits$y, expand = FALSE)))
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
#'   result or a list of them; a lane drawn flipped takes none
#' @param curves draw each ribbon edge as a curve, as the display's
#'   `ribbons = list(mark = list(interpolate = "bump"))` does; unset, as the
#'   result's settings say
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot of the stack alone, or with tracks a patchwork, whose
#'   height and width [jbrowse_height()] and [jbrowse_width()] measure
#' @exportS3Method ggplot2::autoplot
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' liftover <- function(genome, assembly) list(type = "PairwiseIndexedPAFAdapter",
#'   pifGzLocation = list(uri = e(sprintf("hg38To%s.over.pif.gz", genome))),
#'   index = list(location = list(uri = e(sprintf("hg38To%s.over.pif.gz.csi",
#'     genome))), indexType = "CSI"), assemblyNames = c(assembly, "hg38"))
#' primates <- list(type = "MultiPairwiseSyntenyAdapter", adapters = list(
#'   liftover("PanTro6", "panTro6"), liftover("GorGor6", "gorGor6"),
#'   liftover("PonAbe3", "ponAbe3")))
#' m <- jbrowse_multiway(primates, "chr17:7,400,001-7,700,000", "hg38",
#'   lanes = c("panTro6", "gorGor6", "ponAbe3"),
#'   display = list(ribbons = list(encoding = list(color = list(field = "strand")))))
#' autoplot(m)
autoplot.jbrowse_multiway <- function(object, above = NULL, below = NULL,
                                      curves = NULL, ...) {
  rlang::check_dots_empty()
  stack <- multiway_stack(object, curves)
  top <- lapply(listed_tracks(above), lane_band, object, 1, "above")
  bottom <- lapply(listed_tracks(below), lane_band, object, nrow(object$lanes),
    "below")
  if (!length(top) && !length(bottom)) return(stack)
  bands <- c(top, list(list(plots = list(stack), heights = object$height)),
    bottom)
  stack_bands(bands, axis_at = 0)
}
