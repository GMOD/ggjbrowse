# A synteny file's adapter, read off its extension: PAF (plain or gzipped),
# indexed PIF, MUMmer delta, UCSC chain or MashMap output. The adapter names
# its sides `query` and `target`, as the file's columns do, unless `names`
# gives two others, in query, target order.
synteny_adapter <- function(file, index, names = NULL) {
  a <- guess_adapter(file, index, SYNTENY_FORMATS,
    ".paf (optionally .gz), .pif.gz with its index, .delta, .chain or .out")
  sides <- c("query", "target")
  named <- length(names) == 2 && !all(names %in% sides)
  c(a, list(assemblyNames = if (named) names else sides))
}

SYNTENY_FORMATS <- list(
  list(pattern = "\\.paf(\\.gz)?$", type = "PAFAdapter", slot = "pafLocation"),
  list(pattern = "\\.delta$", type = "DeltaAdapter", slot = "deltaLocation"),
  list(pattern = "\\.chain$", type = "ChainAdapter", slot = "chainLocation"),
  list(pattern = "\\.out$", type = "MashMapAdapter", slot = "outLocation"),
  list(pattern = "\\.pif\\.gz$", type = "PairwiseIndexedPAFAdapter",
    slot = "pifGzLocation", index = "TBI"))

# One axis of a comparative view for the engine: its regions, the assembly
# they are on where the caller names it, the pixels it spans, and its side's
# contigs in the order a whole-genome axis lays them, which hand out the
# palette a color by contig paints where no assembly `config` deals it in the
# assembly's own order
axis_args <- function(regions, assembly, width, config = NULL) {
  a <- list(regions = lapply(seq_len(nrow(regions)), function(i)
    list(refName = regions$chrom[i], start = regions$start[i],
      end = regions$end[i])), width = width,
    refNames = I(unique(c(attr(regions, "contigs"), regions$chrom))))
  if (!is.null(assembly) && !is.na(assembly)) a$assemblyName <- assembly
  a$assembly <- config
  a
}

# The engine's alignments, each keyed by `alignment`, the number the tables
# drawn from it name it by
keyed_alignments <- function(columns) {
  d <- snake_frame(columns)
  cbind(alignment = seq_len(nrow(d)), d)
}

snake_frame <- function(columns) snake(as_frame(columns))

# A comparative view's key as the engine titles and colors it, its chips
# drawn at `alpha`
legend_key <- function(r, alpha) {
  structure(engine_key(r$legend, r$legendTitle), alpha = alpha)
}

# The engine's camelCase columns as this package's snake_case
snake <- function(d) {
  names(d) <- gsub("([A-Z])", "_\\L\\1", names(d), perl = TRUE)
  d
}


# The assemblies a pairwise adapter config names, in query, target order
adapter_assemblies <- function(config) {
  unlist(config$assemblyNames) %||%
    unlist(c(config$queryAssembly, config$targetAssembly))
}

comparative_args <- function(file, index, display, names = NULL) {
  args <- list(adapterConfig = adapter_of(file, function()
    synteny_adapter(file, index, names)))
  args$view <- display_config(display)
  args
}

# The engine's tables and the warnings it raised reading them, each also
# raised as an R warning
warned_call <- function(table, args) {
  said <- character(0)
  value <- withCallingHandlers(call_engine(table, args), warning = function(w)
    said <<- c(said, conditionMessage(w)))
  list(value = value, warnings = unique(said))
}

# A comparative result: its tables, then the slots every result holds, then
# its own drawing state
comparative_result <- function(cls, primary, tables, regions, width, display,
                               keys, ramps, warnings, state = list()) {
  structure(c(tables, list(regions = regions, width = width,
      display = display, keys = keys, ramps = ramps, y_scale = NULL,
      fields = list(), warnings = warnings), state),
    primary = primary, class = c(cls, "jbrowse_result"))
}

# The names an axis title shows: each assembly name `labels` renames, as the
# named vector `c(query = "peach", target = "grape")` gives it, or an unnamed
# pair in the adapter's query, target order
display_names <- function(names, labels, config) {
  if (is.null(labels)) return(names)
  if (is.null(names(labels))) {
    if (length(labels) != 2) {
      stop("`labels` is a named vector, or a pair in query, target order",
        call. = FALSE)
    }
    names(labels) <- adapter_assemblies(config)
  }
  hit <- match(names, names(labels))
  ifelse(is.na(hit), names, unname(labels)[hit])
}

# `assemblies` as names, and as the engine's assembly configs where it gives
# them, NULL for an axis it names alone: a pair of names, or a list of two,
# each a name, a jbrowse_assembly() or an assembly config
given_assemblies <- function(assemblies) {
  if (is.null(assemblies) || is.character(assemblies)) {
    return(list(names = assemblies, configs = list(NULL, NULL)))
  }
  if (!is.list(assemblies) || length(assemblies) != 2 ||
    inherits(assemblies, "jbrowse_assembly") || !is.null(assemblies$name)) {
    stop("`assemblies` is two assembly names, or a list of two assemblies",
      call. = FALSE)
  }
  configs <- lapply(unname(assemblies), function(a)
    if (!is.character(a)) assembly_config(a))
  names <- vapply(seq_len(2), function(i) {
    name <- if (is.null(configs[[i]])) assemblies[[i]] else configs[[i]]$name
    if (!is.character(name) || length(name) != 1) {
      stop("each assembly config in `assemblies` has a `name`", call. = FALSE)
    }
    name
  }, "")
  list(names = names, configs = configs)
}

# The assembly each axis is on, first and second, its name and config:
# `given`'s where it names them, else the adapter's own names, reversed by
# `swap`
axis_assemblies <- function(config, given, swap) {
  names <- given$names %||% adapter_assemblies(config)
  configs <- given$configs
  if (!swap) return(list(names = names, configs = configs))
  if (length(names) < 2) {
    stop("swap = TRUE needs two assembly names, the adapter's or `assemblies`",
      call. = FALSE)
  }
  list(names = rev(names[1:2]), configs = rev(configs))
}

# Each side's contigs in the order the file first names them: `side`
# (`"query"` or `"target"`), `assembly`, `chrom` and `length`, `NA` where the
# file does not say
contigs_of <- function(config) {
  r <- call_engine("contigs", list(adapterConfig = config))
  side <- function(name) {
    chrom <- as.character(unlist(r[[name]]$refName))
    data.frame(side = rep(name, length(chrom)),
      assembly = rep(r[[name]]$assemblyName, length(chrom)), chrom = chrom,
      length = as.numeric(unlist(r[[name]]$length)))
  }
  rbind(side("query"), side("target"))
}

memo <- function(f) {
  value <- NULL
  function() {
    if (is.null(value)) value <<- f()
    value
  }
}

# The file's contigs `chrom` in the order a whole-genome axis lays them: the
# assembly `config` lists them in, by their canonical names, then any it does
# not list, in natural order, as all are with no assembly
axis_order <- function(chrom, config) {
  natural <- integer(length(chrom))
  natural[natural_order(chrom)] <- seq_along(chrom)
  if (is.null(config)) return(chrom[order(natural)])
  listed <- match(canonical_names(config, chrom),
    assembly_contigs(config)$canonical)
  chrom[order(listed, natural)]
}

# The regions an axis shows: any region a linear reader takes, a data frame's
# rows narrowed to the axis's assembly where it names an `assembly`, and with
# none given every contig of the axis's assembly in `axis_order()`, which the
# `contigs` attribute holds either way. `side` picks the file's side where the
# axis names no assembly the file has.
axis_regions <- function(given, contigs, assembly, side, config = NULL) {
  all <- contigs()
  named <- length(assembly) && !is.na(assembly) && assembly %in% all$assembly
  name <- if (named) assembly else side
  own <- all[if (named) all$assembly == assembly else all$side == side, ]
  own <- own[!duplicated(own$chrom), ]
  where <- if (name %in% c("query", "target")) paste(name, "side")
    else paste("assembly", name)
  bad <- function(...) rlang::abort(paste0(...), class = "ggjbrowse_bad_region",
      call = NULL)
  contigs <- axis_order(own$chrom, config)
  if (is.null(given)) {
    if (!nrow(own)) bad("the file names no contig on its ", where)
    given <- contigs
  }
  if (is.data.frame(given) && "assembly" %in% names(given)) {
    mine <- given[as.character(given[["assembly"]]) %in% name, , drop = FALSE]
    if (!nrow(mine)) {
      bad("the region frame names assembly ", paste(unique(
        given[["assembly"]]), collapse = ", "), ", none of them this axis's ",
        name)
    }
    given <- mine
  }
  regions <- regions_of(given, function()
    data.frame(chrom = own$chrom, start = 0, end = own$length))
  if (anyNA(regions$end)) {
    bad("the file does not state the length of ", regions$chrom[
      is.na(regions$end)][1], ": give the axis its regions, as loci")
  }
  attr(regions, "contigs") <- contigs
  regions
}

# Both axes of a comparative view: the assembly each is on, its config where
# `given` has one, and its regions
comparative_axes <- function(config, first, second, given, swap) {
  on <- axis_assemblies(config, given, swap)
  contigs <- memo(function() contigs_of(config))
  list(names = on$names, configs = on$configs,
    first = axis_regions(first, contigs, on$names[1], "query", on$configs[[1]]),
    second = axis_regions(second, contigs, on$names[2], "target",
      on$configs[[2]]))
}

# An axis's rows of a result's `regions`: each region named by the `axis` it
# lays out on and the `assembly` it is on, laid end to end from the axis's
# first region as JBrowse's synteny and dotplot workers lay them, with no gap
# unless `gap_frac` gives one
axis_rows <- function(regions, axis, assembly, gap_frac = 0) {
  l <- region_layout(regions[c("chrom", "start", "end")], gap_frac)
  n <- nrow(l)
  data.frame(axis = rep(axis, n),
    assembly = rep(if (length(assembly)) assembly else NA_character_, n), l,
    missing = regions$missing %||% rep("", n))
}

# The rows of `regions` on `axis`
on_axis_of <- function(regions, axis) {
  r <- regions[regions$axis == axis, , drop = FALSE]
  rownames(r) <- NULL
  r
}

#' Read a pairwise alignment as JBrowse's linear synteny view draws it
#'
#' Runs JBrowse's own synteny display in V8 over two rows, the query assembly
#' above the target: the alignment blocks the two rows share, the ribbons the
#' band draws between them, with each CIGAR insertion, deletion and skip wide
#' enough to show as its own colored piece, and the color and opacity the
#' view's settings give each. Coordinates are 0-based, half-open, on each
#' assembly.
#'
#' @param file a PAF (plain or gzipped), an indexed PIF (`.pif.gz` with its
#'   tabix index beside it), a MUMmer `.delta`, a UCSC `.chain` or a MashMap
#'   `.out` file, a path or a URL; or a JBrowse synteny track's `adapter`
#'   config, as a list, such as an `MCScanAnchorsAdapter`'s
#' @param query,target the regions each row shows, drawn end to end as a
#'   view lays them out, in any form a linear reader's `region` takes: loci
#'   as JBrowse's location box writes them (`"ctgA:1-50,001"`, 1-based and
#'   inclusive), contig names alone for whole contigs, a data frame of
#'   `chrom`, `start` and `end` (0-based, half-open), or a GRanges. A frame
#'   with an `assembly` column, as [jbrowse_contigs()] and a result's
#'   `regions` write it, gives the row the regions on its assembly. With
#'   `NULL`, the default, every contig of that row's assembly, in the order
#'   the assembly lists them where `assemblies` gives it, else in natural
#'   order, `chr2` before `chr10`, as JBrowse shows a whole genome. A region
#'   on a contig the file lacks reads empty and warns, and rows with none of
#'   the file's contigs between them stop
#' @param display the view's settings, as a `LinearSyntenyView` snapshot
#'   writes them: a list with any of `encoding` (its `color`, a CSS color or
#'   the color object's slots, such as `list(field = "strand")`,
#'   `list(field = "identity")`, `list(field = "query")` or
#'   `list(field = "mapq")`, and its `opacity`, a number or an opacity object
#'   with a `field`, so `list(encoding = list(color = list(field =
#'   "strand"), opacity = 0.5))`), `filter` (a list of field predicates, such
#'   as `list(list(field = "length", gte = 10000))`), `indels` (unset for the
#'   indels in their colors, `list(encoding = list(opacity = 0))` for the
#'   matches alone with the indels left clear, `list(show = FALSE)` for the
#'   block alone), `mark` (`list(interpolate = "bump")` curves each ribbon's
#'   sides), `axis` (`list(grid = TRUE)` draws the location markers), `render`
#'   (`list(antialias = FALSE)` stops the fade of sub-pixel ribbons; unset,
#'   it fades once they dominate) and `height`, the band's pixels. Any other
#'   key is an error naming these
#' @param width the width in pixels both rows draw their regions across,
#'   1,500 as for the linear readers, which decides which CIGAR ops are wide
#'   enough to draw and which ribbons are thin enough to fade
#' @param assemblies which assembly each row is on, `c(query, target)`: two
#'   names, or a list of two, each a name or an assembly, a
#'   [jbrowse_assembly()] or an assembly config as a list, named by its
#'   `name`. An assembly lays a whole-genome row in its own contig order and
#'   deals a color by contig its palette in that order, aliases included, as
#'   the view does. For an adapter config naming its assemblies, the names
#'   pick its sides, so `c("volvox", "volvox_del")` on a PAF of `volvox_del`
#'   against `volvox` puts `volvox` on top; a file guessed from its
#'   extension takes them as its query's and target's names, unless they are
#'   `"query"` and `"target"` themselves. Unset, the adapter's own order
#' @param swap put the file's target side on the top row and its query side
#'   on the bottom: the adapter's two assembly names, or `assemblies`,
#'   reversed. A file guessed from its extension names its sides `"query"`
#'   and `"target"`, so `swap = TRUE` is `assemblies = c("target",
#'   "query")`
#' @param labels what each axis title calls its assembly, named by the
#'   assembly: `c(query = "peach", target = "grape")` for a file guessed from
#'   its extension, or the adapter's own assembly names. An unnamed pair is
#'   the adapter's query and target in that order. Unset, a title shows the
#'   assembly name itself, `"query"` and `"target"` for a guessed file
#' @param index the tabix or CSI index of a `.pif.gz`, if not beside it
#' @return a `jbrowse_synteny` list, beside the slots [jbrowse_result] lists:
#'   * `alignments`: each block, in the band's paint order: `alignment`,
#'     its key, `id`, `name`,
#'     `query_chrom`, `query_start`, `query_end`, `target_chrom`,
#'     `target_start`, `target_end`, `strand`, `length`, the measurements the
#'     file carries (`identity`, `mean_identity`, `mapq`, and any column the
#'     adapter's `attributeColumns` declares), `color` and `alpha`
#'   * `ribbons`: what the band draws, in paint order, each a trapezoid from
#'     a span on the query row to a span on the target row: `alignment` (the
#'     block's key in `alignments`), `kind` (`"match"`, `"insertion"`,
#'     `"deletion"`, `"skip"`, or `"marker"` for a location tick),
#'     `query_region`, `query_chrom`, `query_start`, `query_end`,
#'     `target_region`, `target_chrom`, `target_start`, `target_end`, each
#'     `*_region` the row of `regions` its foot lies in, where
#'     the ribbon joins `query_start` to `target_start` and `query_end` to
#'     `target_end`, so a reverse-strand block crosses; `query_x1`,
#'     `query_x2`, `target_x1` and `target_x2`, the pixel across `width` of
#'     each foot at `*_start` and `*_end`, within its region; `color`,
#'     `alpha`, and `thin`, whether the band draws it as a one-pixel line
#'   * `offscreen`: the blocks on one `axis` whose mate the other row is not
#'     showing, on a contig it does not show or past its window and the
#'     band's overdraw, the `region` they lie in, their span and their
#'     mate's, the `color` and `alpha` the strip along that row's edge paints
#'     them, and `drawn`, whether they clear the strip's length and contig
#'     floors
#'   * `offscreen_contigs`: each contig those mates lie on, by `axis`:
#'     `mate_chrom`, `count`, `aligned_bp`, the `label` the strip's hover
#'     prints, `color` and `alpha`
#'   * `offscreen_strip`: what the strips paint, in pixels from the band's
#'     top left: `marks` (`axis`, `x`, `y`, `width`, `height`, `color`,
#'     `alpha`) and `labels` (`axis`, `label`, `x`, `y` at the baseline,
#'     `color`, `alpha`, `halo`), in a `font_px` font
#'
#'   Assigning `alignments` drops the ribbons of the blocks it leaves out;
#'   the off-screen tables summarize the band's edges as read, and stay.
#'   `keys$ribbons` is the key the view shows, `ramps` a continuous field's
#'   color bar, as [guides_jbrowse()] draws it. `regions` lays out both rows, its `axis` `"query"` for the top row and
#'   `"target"` for the bottom. `display` holds the view's settings as the
#'   engine resolves them, `encoding$color` and `encoding$opacity` as config
#'   objects and `height` the band's pixels.
#' @export
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' s <- jbrowse_synteny(chain,
#'   "chr8_PATERNAL:5,250,001-14,250,000",
#'   "chr8_MATERNAL:5,250,001-14,250,000")
#' head(s$alignments)
jbrowse_synteny <- function(file, query = NULL, target = NULL,
                            display = NULL, width = DEFAULT_WIDTH,
                            assemblies = NULL, swap = FALSE, labels = NULL,
                            index = NULL) {
  p <- pairwise_read("synteny", file, query, target,
    c(query = width, target = width), display, assemblies, swap, labels, index)
  r <- p$r
  offscreen <- snake_frame(r$offscreen)
  offscreen$region <- offscreen$region + 1 +
    ifelse(offscreen$axis == "target", p$first_regions, 0)
  contigs <- snake_frame(r$offscreenContigs)
  comparative_result("jbrowse_synteny", "ribbons",
    list(alignments = p$alignments, ribbons = p$children(r$ribbons),
      offscreen = offscreen, offscreen_contigs = contigs),
    p$regions, width, r$config,
    c(list(ribbons = p$key), offscreen_keys(offscreen, contigs)),
    ramp_frames(r$ramps), p$warnings,
    list(offscreen_strip = list(marks = snake_frame(r$offscreenStrip$marks),
        labels = snake_frame(r$offscreenStrip$labels),
        font_px = r$offscreenStrip$fontPx),
      assemblies = p$assemblies))
}

# A pairwise view's read: the engine's `view` tables of `file` over the
# `first` and `second` axes, named and sized by `px`, with what both views
# make of them. `children()` keys a table's rows to their alignment and each
# axis's regions, the second's after the first's. Axes none of whose contigs
# the file has stop; a region on a contig the file lacks beside others warns,
# its warning in `warnings`. A row with no contig of the file still draws the
# strip of the other row's off-screen mates.
pairwise_read <- function(view, file, first, second, px, display, assemblies,
                          swap, labels, index) {
  given <- given_assemblies(assemblies)
  args <- comparative_args(file, index, display, given$names)
  axes <- comparative_axes(args$adapterConfig, first, second, given, swap)
  sides <- names(px)
  args[[sides[1]]] <- axis_args(axes$first, axes$names[1], px[[1]],
    axes$configs[[1]])
  args[[sides[2]]] <- axis_args(axes$second, axes$names[2], px[[2]],
    axes$configs[[2]])
  read <- warned_call(view, args)
  r <- read$value
  said <- function(side) as.character(unlist(r$missingContigs[[side]]))
  missing <- check_missing(c(said(sides[1]), said(sides[2])),
    c(axes$first$chrom, axes$second$chrom))
  axes$first$missing <- said(sides[1])
  axes$second$missing <- said(sides[2])
  first_regions <- nrow(axes$first)
  children <- function(table) {
    d <- snake_frame(table)
    d$alignment <- as.integer(d$alignment) + 1L
    on_first <- paste0(sides[1], "_region")
    on_second <- paste0(sides[2], "_region")
    d[[on_first]] <- d[[on_first]] + 1
    d[[on_second]] <- d[[on_second]] + 1 + first_regions
    d
  }
  list(r = r, children = children, first_regions = first_regions,
    regions = rbind(axis_rows(axes$first, sides[1], axes$names[1]),
      axis_rows(axes$second, sides[2], axes$names[2])),
    alignments = keyed_alignments(r$alignments),
    key = legend_key(r, r$legendAlpha),
    warnings = unique(c(attr(missing, "warning"), read$warnings)),
    assemblies = display_names(unlist(r$assemblyNames), labels,
      args$adapterConfig))
}

#' @export
print.jbrowse_synteny <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_synteny> %s over %s, %d alignments, %d ribbons, %d off screen\n",
    loci_text(on_axis_of(x$regions, "query")),
    loci_text(on_axis_of(x$regions, "target")),
    nrow(x$alignments), nrow(x$ribbons), nrow(x$offscreen)))
  invisible(x)
}

# A position in a row of `regions` as a pixel of the band: its place on its
# axis's regions laid end to end, as a share of the axis's span, across the
# view's width
band_px <- function(regions, axis, region, coord, width) {
  lo <- tapply(regions$cum_start, regions$axis, min)[axis]
  hi <- tapply(regions$cum_end, regions$axis, max)[axis]
  unname((on_offsets(regions, region, coord) - lo) / (hi - lo) * width)
}

# The result with its ribbons and off-screen blocks marked as on the band's
# axis, where the engine placed each foot
synteny_on_axis <- function(x) {
  r <- x$ribbons
  if (is_on_axis(r)) return(x)
  x$ribbons <- structure(r, on_axis = TRUE)
  x$offscreen <- structure(x$offscreen, on_axis = TRUE)
  x
}

# The points of each ribbon edge from (x0, y0) to (x1, y1), a row per edge:
# a straight line, or the band's cubic Bezier with both control points at mid
# height
edge_points <- function(x0, x1, y0, y1, curves, n = 16) {
  at <- if (curves) seq(0, 1, length.out = n) else c(0, 1)
  ease_x <- if (curves) at * at * (3 - 2 * at) else at
  ease_y <- if (curves) 1.5 * at * (1 - at) + at^3 else at
  list(x = x0 + outer(x1 - x0, ease_x), y = y0 + outer(y1 - y0, ease_y))
}

# Where each ribbon's feet sit: px along its upper and lower edges, and the
# heights of those edges. The band's query row is at y = 1, its target at 0.
band_feet <- function(r) {
  list(top1 = r$query_x1, top2 = r$query_x2,
    bottom1 = r$target_x1, bottom2 = r$target_x2,
    y_top = rep(1, nrow(r)), y_bottom = rep(0, nrow(r)))
}

# Ribbons `r` as one frame of points, ribbon by ribbon along the rows of `x`
# and `y`, each point carrying every column of its ribbon
ribbon_marks <- function(r, x, y) {
  d <- data.frame(x = as.vector(t(x)), y = as.vector(t(y)),
    ribbon = rep(seq_len(nrow(r)), each = ncol(x)))
  r <- r[setdiff(names(r), names(d))]
  cbind(d, r[d$ribbon, , drop = FALSE], row.names = NULL)
}

# Each filled ribbon's outline, top to bottom on its first edge and back up
# its second
ribbon_polygons <- function(r, curves, feet = band_feet(r)) {
  down <- edge_points(feet$top1, feet$bottom1, feet$y_top, feet$y_bottom,
    curves)
  up <- edge_points(feet$bottom2, feet$top2, feet$y_bottom, feet$y_top,
    curves)
  ribbon_marks(r, cbind(down$x, up$x), cbind(down$y, up$y))
}

# Each thin ribbon's centerline, from the middle of its top edge to the
# middle of its bottom edge
ribbon_lines <- function(r, curves, feet = band_feet(r)) {
  e <- edge_points((feet$top1 + feet$top2) / 2,
    (feet$bottom1 + feet$bottom2) / 2, feet$y_top, feet$y_bottom, curves)
  ribbon_marks(r, e$x, e$y)
}

#' A synteny band's ribbons as ggplot layers
#'
#' `geom_jbrowse_ribbons()` draws each ribbon between the query row at the
#' top of the panel (`y = 1`) and the target row at the bottom (`y = 0`) as
#' the band draws it: a filled trapezoid in its color at its opacity, or a
#' one-pixel line where it is thinner than a pixel (`linewidth` in `...`
#' widens it), with the band's location
#' ticks where the view draws them. `geom_jbrowse_offscreen_mates()` draws
#' the strip the band hangs off a row's edge for the blocks whose mate the
#' other row is not showing: the display's own marks, 6 px tall, grey or in
#' the mate contig's ribbon color where the ribbons are colored by that
#' row's contigs, and the contig names it has room for.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_synteny()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param curves draw each ribbon edge as the view's
#'   `mark = list(interpolate = "bump")` does, a cubic curve vertical at both
#'   rows; unset, as the result's display says
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' s <- jbrowse_synteny(chain,
#'   "chr8_PATERNAL:5,250,001-14,250,000",
#'   "chr8_MATERNAL:5,250,001-14,250,000")
#' library(ggplot2)
#' ggplot(s) + geom_jbrowse_ribbons(curves = TRUE) + guides_jbrowse(s)
geom_jbrowse_ribbons <- function(mapping = NULL, data = NULL, curves = NULL,
                                 ...) {
  ribbon_layers(mapping, replacing(data, "ribbons", "geom_jbrowse_ribbons"),
    function(x) curves %||% bump(x$display$mark), band_feet, ...)
}

# Whether a ribbons `mark` config curves each edge
bump <- function(mark) identical(mark$interpolate, "bump")

# A ribbon geom's two layers: the ribbons the engine flags `thin` as lines,
# the rest as polygons, each between `feet(ribbons)` and curved where
# `curved(result)` says
ribbon_layers <- function(mapping, data, curved, feet, ...) {
  pick <- function(x, thin) x$ribbons[x$ribbons$thin == thin, , drop = FALSE]
  filled <- function(x) {
    r <- pick(x, FALSE)
    ribbon_polygons(r, curved(x), feet(r))
  }
  lines <- function(x) {
    r <- pick(x, TRUE)
    ribbon_lines(r, curved(x), feet(r))
  }
  list(
    jbrowse_layer(ggplot2::GeomPolygon, data, filled,
      for_geom(mapping, ggplot2::GeomPolygon),
      ggplot2::aes(x = x, y = y, group = ribbon, fill = I(color),
        alpha = I(alpha)), list(colour = NA, ...)),
    jbrowse_layer(ggplot2::GeomPath, data, lines,
      for_geom(mapping, ggplot2::GeomPath),
      ggplot2::aes(x = x, y = y, group = ribbon, colour = I(color),
        alpha = I(alpha)), with_defaults(list(...),
        list(linewidth = px_to_mm(1)))))
}

#' @rdname geom_jbrowse_ribbons
#' @param labels name each stretch of marks with its contig where the strip
#'   has room, as the band does
#' @export
geom_jbrowse_offscreen_mates <- function(mapping = NULL, data = NULL,
                                         labels = TRUE, ...) {
  data <- replacing(data, NULL, "geom_jbrowse_offscreen_mates")
  c(list(jbrowse_layer(ggplot2::GeomRect, data, strip_marks,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = I(color), alpha = I(alpha)),
      list(colour = NA, ...))),
    if (labels) list(jbrowse_layer(GeomDisplayText, data, strip_labels, NULL,
      ggplot2::aes(x = x, y = y, label = label, colour = I(color),
        alpha = I(alpha), size = I(size), px_wide = I(px_wide),
        px_tall = I(px_tall), halo = I(halo), halo_px = I(HALO_PX)),
      list(hjust = 0, vjust = 0))))
}

# The strip's marks on the band, which runs from the target row (y = 0) up
# to the query row (y = 1) over `height` pixels
strip_marks <- function(x) {
  m <- x$offscreen_strip$marks
  m$xmin <- m$x
  m$xmax <- m$x + m$width
  px <- x$display$height
  m$ymax <- 1 - m$y / px
  m$ymin <- 1 - (m$y + m$height) / px
  m
}

HALO_PX <- 3

# The strip's labels at their baselines, each haloed in its `halo` color
strip_labels <- function(x) {
  l <- x$offscreen_strip$labels
  n <- nrow(l)
  px <- x$display$height
  l$y <- 1 - l$y / px
  l$size <- rep(x$offscreen_strip$font_px, n)
  l$px_wide <- rep(x$width, n)
  l$px_tall <- rep(px, n)
  l
}

OFFSCREEN_KEY_ROWS <- 10

# A key per row whose strip draws marks: a chip of each contig in the color
# the strip paints it, beside the line the strip's hover prints for it, the
# most sequence first
offscreen_keys <- function(offscreen, contigs) {
  drawn <- offscreen[offscreen$drawn, , drop = FALSE]
  axes <- intersect(c("query", "target"), drawn$axis)
  keys <- lapply(axes, function(ax) {
    on <- contigs[contigs$axis == ax &
      contigs$mate_chrom %in% drawn$mate_chrom[drawn$axis == ax], ]
    on <- on[order(-on$aligned_bp), ]
    labels <- on$label
    fill <- scales::alpha(on$color, on$alpha)
    rest <- length(labels) - OFFSCREEN_KEY_ROWS
    if (rest > 0) {
      labels <- c(labels[seq_len(OFFSCREEN_KEY_ROWS)], sprintf("+%d more", rest))
      fill <- c(fill[seq_len(OFFSCREEN_KEY_ROWS)], NA)
    }
    key_frame(labels, fill, sprintf("Off-screen mates, %s row", ax))
  })
  stats::setNames(keys, sprintf("offscreen_%s", axes))
}

# Whether an axis names its contigs rather than numbering positions, as
# JBrowse's dotplot labels each region by name: an axis of more than four
# regions, or of several where one spans under a tenth of it
names_contigs <- function(layout) {
  w <- layout$end - layout$start
  nrow(layout) > 4 || (nrow(layout) > 1 && min(w) < sum(w) / 10)
}

# The pixels apart two contig names stand at least, the dotplot's LABEL_PX
LABEL_PX <- 12

# Which labels at `at` px to keep: each its region's, widest region first,
# where it stands LABEL_PX clear of every label kept before it
thin_labels <- function(at, widths) {
  kept <- numeric(0)
  keep <- logical(length(at))
  for (i in order(widths, decreasing = TRUE)) {
    if (all(abs(kept - at[i]) >= LABEL_PX)) {
      kept <- c(kept, at[i])
      keep[i] <- TRUE
    }
  }
  keep
}

# The labels on an axis of `width` px, `px(layout, region, coord)` placing a
# position on it: round positions in each region, labelled in the region's
# own coordinates, or where names_contigs() says, each contig's name at its
# middle
layout_breaks <- function(layout, px, width) {
  r <- layout
  if (names_contigs(r)) {
    n <- nrow(r)
    mid <- px(r, seq_len(n), (r$start + r$end) / 2)
    edges <- px(r, c(1, n), c(r$start[1], r$end[n]))
    keep <- thin_labels((mid - edges[1]) / diff(edges) * width,
      r$end - r$start)
    list(breaks = mid[keep], labels = r$chrom[keep])
  } else region_ticks(r, px)
}

layout_label <- function(layout, assembly = "") {
  if (names_contigs(layout)) {
    return(sprintf("%s%d regions, end to end", if (nzchar(assembly))
      paste0(assembly, ", ") else "", nrow(layout)))
  }
  what <- if (nrow(layout) == 1) layout$chrom
    else loci_text(layout)
  name <- if (nzchar(assembly)) paste0(assembly, " ") else ""
  sprintf("%s%s (bp%s)", name, what,
    if (nrow(layout) == 1) "" else ", end to end")
}

# The band's two rulers: the target row's positions under the panel and the
# query row's above it, each at its pixel
synteny_x <- function(x) {
  on_row <- function(axis) layout_breaks(on_axis_of(x$regions, axis),
    function(r, i, at) band_px(r, axis, i, at, x$width), x$width)
  q <- on_row("query")
  t <- on_row("target")
  ggplot2::scale_x_continuous(breaks = t$breaks, labels = t$labels,
    name = layout_label(on_axis_of(x$regions, "target"), x$assemblies[2]),
    sec.axis = ggplot2::sec_axis(~., breaks = q$breaks, labels = q$labels,
      name = layout_label(on_axis_of(x$regions, "query"), x$assemblies[1])))
}

synteny_band <- function(x, curves) {
  x <- on_axis(x)
  ggplot2::ggplot(x) +
    geom_jbrowse_ribbons(curves = curves) +
    geom_jbrowse_offscreen_mates() +
    guides_jbrowse(x) +
    ggplot2::coord_cartesian(xlim = c(0, x$width), ylim = c(0, 1),
      expand = FALSE) +
    synteny_x(x) +
    ggplot2::labs(y = NULL) +
    ggplot2::theme_minimal() +
    no_y +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    upright_names(on_axis_of(x$regions, "target"), "x.bottom") +
    upright_names(on_axis_of(x$regions, "query"), "x.top")
}

# An axis's contig names turned to run along the axis's normal, as the
# dotplot's horizontal axis draws them, so close names do not overlap
upright_names <- function(layout, axis) {
  if (!names_contigs(layout)) return(NULL)
  hjust <- if (axis == "x.top") 0 else 1
  do.call(ggplot2::theme, stats::setNames(list(ggplot2::element_text(
    angle = 90, hjust = hjust, vjust = 0.5)), paste0("axis.text.", axis)))
}

# The tracks a row shows, each read over the row's regions
row_tracks <- function(tracks, layout, row) {
  tracks <- named_tracks(tracks, "")
  for (t in tracks) {
    if (!same_loci(t$regions, layout)) {
      stop("a ", row, " track draws over the ", row, " row's regions, ",
        loci_text(layout), ", but one reads ",
        loci_text(t$regions), call. = FALSE)
    }
  }
  lapply(tracks, function(t) {
    axis <- c("offset", "cum_start", "cum_end")
    t$regions[axis] <- layout[axis]
    t
  })
}

#' Draw a pairwise alignment as JBrowse's linear synteny view
#'
#' The band of ribbons between the two rows, the query row's ruler above it
#' and the target row's below, and any tracks over each row: a
#' [jbrowse_features()], [jbrowse_wiggle()], [jbrowse_variants()] or
#' [jbrowse_alignments()] result over that row's regions, drawn as
#' [jbrowse_stack()] draws it, the top row's tracks above the band and the
#' bottom row's below.
#'
#' @param object a [jbrowse_synteny()] result
#' @param above,below tracks over the top row's regions, drawn above the
#'   band, and over the bottom row's, drawn below it: each a result or a
#'   list of them, a list's names titling its tracks as [jbrowse_stack()]
#'   titles them
#' @param curves draw each ribbon edge as a curve, as the view's
#'   `mark = list(interpolate = "bump")` does; unset, as the result's display
#'   says
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot of the band alone, or with tracks a patchwork, whose
#'   height and width [jbrowse_height()] and [jbrowse_width()] measure
#' @exportS3Method ggplot2::autoplot
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' s <- jbrowse_synteny(e("hg002_chr8.chain"), "chr8_PATERNAL:5,250,001-14,250,000",
#'   "chr8_MATERNAL:5,250,001-14,250,000",
#'   display = list(encoding = list(color = list(field = "strand"))),
#'   labels = c(query = "paternal", target = "maternal"))
#' autoplot(s)
#' genes <- jbrowse_features(e("hg002_mat_genes.gff3.gz"), "chr8_MATERNAL:5,250,001-14,250,000",
#'   display = list(filter = "geneLike"))
#' autoplot(s, below = genes)
autoplot.jbrowse_synteny <- function(object, above = NULL, below = NULL,
                                     curves = NULL, ...) {
  rlang::check_dots_empty()
  band <- synteny_band(object, curves)
  row_bands <- function(tracks, axis, row) titled_bands(
    row_tracks(tracks, on_axis_of(object$regions, axis), row), track_bands)
  around_band(band, object$display$height, row_bands(above, "query", "top"),
    row_bands(below, "target", "bottom"), main_axis = TRUE)
}

#' Read a pairwise alignment as JBrowse's dotplot view draws it
#'
#' Runs JBrowse's own dotplot display in V8: the alignment blocks between the
#' two axes, and the lines the plot draws for them, one per block or, where a
#' block is wide enough on the plot for its CIGAR to show, one per CIGAR
#' step, each in the color the view's settings paint it. Coordinates are
#' 0-based, half-open, on each assembly.
#'
#' @inheritParams jbrowse_synteny
#' @param x,y the regions each axis shows, written as [jbrowse_synteny()]'s
#'   `query` and `target` are: the horizontal axis is the adapter's query
#'   side and the vertical its target, unless `assemblies` or `swap` says
#'   otherwise. Unset, every contig of that axis's assembly, so
#'   `jbrowse_dotplot(paf)` is the whole genome against the whole genome, in
#'   the order [jbrowse_synteny()]'s `query` describes
#' @param display the view's settings, as a `DotplotView` snapshot writes
#'   them: a list with any of `encoding` (its `color` and `opacity`, as for
#'   [jbrowse_synteny()]), `filter` (a list of field predicates, such as
#'   `list(list(field = "identity", gte = 0.9))`), `indels` (as for
#'   [jbrowse_synteny()]) and `lineWidth`. Any other key is an error naming
#'   these
#' @param width,height the plot's pixels across and up, which decide which
#'   blocks are wide enough for their CIGAR to show: 1,500 across as for the
#'   linear readers, and up as many as across
#' @param assemblies which assembly each axis is on, `c(x, y)`, in any form
#'   [jbrowse_synteny()]'s `assemblies` takes
#' @param swap put the file's target side on the horizontal axis and its
#'   query on the vertical, as [jbrowse_synteny()]'s `swap` does the rows
#' @param labels what each axis title calls its assembly, as for
#'   [jbrowse_synteny()]: `c(query = "peach", target = "grape")`
#' @return a `jbrowse_dotplot` list, beside the slots [jbrowse_result] lists:
#'   * `alignments`: each block: `alignment`, its key, `name`, `x_chrom`, `x_start`, `x_end`,
#'     `y_chrom`, `y_start`, `y_end`, `strand`, `length`, the measurements
#'     the file carries, `color` and `alpha`
#'   * `segments`: the lines the plot draws, in paint order: `alignment` (the
#'     block's key in `alignments`), `op` (`"match"`, `"insertion"`,
#'     `"deletion"` or `"skip"`), `x_region`, `x_chrom`, `x_start`, `x_end`,
#'     `y_region`, `y_chrom`, `y_start`, `y_end`, from (`x_start`, `y_start`)
#'     to (`x_end`, `y_end`), each `*_region` a row of `regions`, `color` and
#'     `alpha`
#'
#'   Assigning `alignments` drops the segments of the blocks it leaves out.
#'   `keys$segments` is the key the view shows, `ramps` a continuous field's
#'   color bar. `regions` lays out both axes, its `axis` `"x"` or `"y"`. `display`
#'   holds the view's settings as the engine resolves them, `encoding$color`
#'   and `encoding$opacity` as config objects and `lineWidth` a line's pixels. `skipped`
#'   is the blocks neither axis places: their `count` of the `total` read,
#'   and the contigs they lie on, `x_chroms` and `y_chroms`.
#' @export
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' d <- jbrowse_dotplot(chain,
#'   "chr8_PATERNAL:5,250,001-14,250,000",
#'   "chr8_MATERNAL:5,250,001-14,250,000")
#' head(d$segments)
jbrowse_dotplot <- function(file, x = NULL, y = NULL, display = NULL,
                            width = DEFAULT_WIDTH, height = width,
                            assemblies = NULL, swap = FALSE, labels = NULL,
                            index = NULL) {
  p <- pairwise_read("dotplot", file, x, y, c(x = width, y = height), display,
    assemblies, swap, labels, index)
  r <- p$r
  comparative_result("jbrowse_dotplot", "segments",
    list(alignments = p$alignments, segments = p$children(r$segments)),
    p$regions, width, r$config, list(segments = p$key), ramp_frames(r$ramps),
    p$warnings,
    list(height = height, assemblies = p$assemblies, skipped = list(
      count = r$skipped$count, total = r$skipped$total,
      x_chroms = as.character(unlist(r$skipped$xChroms)),
      y_chroms = as.character(unlist(r$skipped$yChroms)))))
}

#' @export
print.jbrowse_dotplot <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_dotplot> %s against %s, %d alignments, %d segments\n",
    loci_text(on_axis_of(x$regions, "x")),
    loci_text(on_axis_of(x$regions, "y")),
    nrow(x$alignments), nrow(x$segments)))
  invisible(x)
}

# Each segment's ends on the axes that lay each side's regions end to end
dotplot_on_axis <- function(x) {
  s <- x$segments
  if (is_on_axis(s)) return(x)
  r <- x$regions
  s$x_start <- on_offsets(r, s$x_region, s$x_start)
  s$x_end <- on_offsets(r, s$x_region, s$x_end)
  s$y_start <- on_offsets(r, s$y_region, s$y_start)
  s$y_end <- on_offsets(r, s$y_region, s$y_end)
  x$segments <- structure(s, on_axis = TRUE)
  x
}

#' A dotplot's lines as a ggplot layer
#'
#' Draws each segment from its start to its end in its color at its
#' opacity, with round ends, as the view does.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_dotplot()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param ... passed to [ggplot2::layer()] params, such as `linewidth`
#' @return a ggplot layer
#' @export
#' @examples
#' paf <- system.file("extdata", "hprc_amylase.paf", package = "ggjbrowse")
#' d <- jbrowse_dotplot(paf, "GRCh38#0#chr1:103,520,894-103,832,637",
#'   "HG00133#1#CM090045.1:103,669,666-103,981,330")
#' library(ggplot2)
#' ggplot(d) + geom_jbrowse_dotplot(aes(colour = op))
geom_jbrowse_dotplot <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "segments", "geom_jbrowse_dotplot")
  jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$segments,
    for_geom(mapping, ggplot2::GeomSegment),
    ggplot2::aes(x = x_start, xend = x_end, y = y_start, yend = y_end,
      colour = I(color), alpha = I(alpha)),
    with_defaults(list(...), list(linewidth = px_to_mm(2.5),
      lineend = "round")))
}

# An axis of regions laid end to end: breaks at round positions in each,
# labelled in the region's own coordinates
layout_scale <- function(scale, layout, width) {
  b <- layout_breaks(layout, on_offsets, width)
  scale(breaks = b$breaks, labels = b$labels)
}

# The seams between an axis's regions
layout_seams <- function(layout) {
  if (nrow(layout) < 2) return(numeric(0))
  (layout$cum_end[-nrow(layout)] + layout$cum_start[-1]) / 2
}

#' Draw a pairwise alignment as JBrowse's dotplot
#'
#' [geom_jbrowse_dotplot()] at the view's line width, each axis labelled in
#' its regions' positions, with a line at each seam where an axis lays
#' several regions end to end.
#' @param object a [jbrowse_dotplot()] result
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot
#' @exportS3Method ggplot2::autoplot
#' @examples
#' chain <- system.file("extdata", "hg002_chr8.chain", package = "ggjbrowse")
#' autoplot(jbrowse_dotplot(chain, display = list(encoding = list(color = list(field = "strand"))),
#'   labels = c(query = "paternal", target = "maternal")))
autoplot.jbrowse_dotplot <- function(object, ...) {
  rlang::check_dots_empty()
  x <- on_axis(object)
  h <- on_axis_of(x$regions, "x")
  v <- on_axis_of(x$regions, "y")
  seams <- data.frame(x = layout_seams(h))
  rows <- data.frame(y = layout_seams(v))
  ggplot2::ggplot(x) +
    ggplot2::geom_vline(data = seams, ggplot2::aes(xintercept = x),
      colour = "grey70") +
    ggplot2::geom_hline(data = rows, ggplot2::aes(yintercept = y),
      colour = "grey70") +
    geom_jbrowse_dotplot(linewidth = px_to_mm(x$display$lineWidth)) +
    guides_jbrowse(x) +
    layout_scale(ggplot2::scale_x_continuous, h, x$width) +
    layout_scale(ggplot2::scale_y_continuous, v, x$height) +
    ggplot2::coord_cartesian(xlim = c(min(h$cum_start), max(h$cum_end)),
      ylim = c(min(v$cum_start), max(v$cum_end)), expand = FALSE) +
    ggplot2::labs(x = layout_label(h, x$assemblies[1]),
      y = layout_label(v, x$assemblies[2])) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey70")) +
    upright_names(h, "x.bottom")
}

# Names in the order a reader expects them: chr2 before chr10, each run of
# digits compared as a number
natural_order <- function(x) {
  at <- gregexpr("[0-9]+", x)
  runs <- regmatches(x, at)
  width <- max(0, nchar(unlist(runs)))
  key <- x
  regmatches(key, at) <- lapply(runs, function(r)
    paste0(strrep("0", width - nchar(r)), r))
  order(key, x)
}
