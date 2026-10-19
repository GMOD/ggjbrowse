VCF_FORMATS <- list(
  list(pattern = "\\.vcf$", type = "VcfAdapter", slot = "vcfLocation"),
  list(pattern = "\\.vcf\\.b?gz$", type = "VcfTabixAdapter",
    slot = "vcfGzLocation", index = "TBI"))

# A VCF's adapter, read off its extension, any other name read as bgzipped
variants_adapter <- function(file, index, samples) {
  adapter <- guess_adapter(file, index, c(VCF_FORMATS, list(list(pattern = "",
    type = "VcfTabixAdapter", slot = "vcfGzLocation", index = "TBI"))))
  if (!is.null(samples)) adapter$samplesTsvLocation <- location(samples)
  adapter
}

#' Read a region of a multi-sample VCF as JBrowse paints it
#'
#' Runs JBrowse's own multi-sample variant display in V8: the records its
#' filters keep, each sample's genotype at each, and the fill the display
#' paints each cell, its hue shaded by the alt dosage. Coordinates are 0-based,
#' half-open, on the reference.
#'
#' The display keeps a record only where some sample calls an allele, so a
#' calls VCF whose FORMAT has no GT, such as a somatic SV caller's, keeps
#' none. [jbrowse_features()] reads its records as JBrowse's plain variant
#' track draws them, with their FILTER, INFO and FORMAT columns.
#'
#' @param file a VCF, bgzipped with a tabix or CSI index beside it, or plain;
#'   a path or a URL; or a JBrowse track's `adapter` config, as a list
#' @param region the regions to read, read as one view and drawn end to end
#'   on one axis as a JBrowse view lays them out: loci as JBrowse's location
#'   box writes them (`"ctgA:1,001-9,000"`, 1-based and inclusive), contig
#'   names alone for whole contigs, a data frame of `chrom`, `start` and `end`
#'   (0-based, half-open, as BED writes them), or a GRanges. With `NULL`, the
#'   default, the reader reads the whole file: each contig it lists but those
#'   its index shows empty, to the length the file or `assembly` states. A whole
#'   file or contig is refused past 5 million rows in any one table
#' @param display the display's settings, as a track config's `displays`
#'   entry writes them for a `LinearMultiSampleVariantDisplay`: a list with
#'   any of `unit` (`"sample"` for a row per sample, `"haplotype"` for a row
#'   per haplotype), `color` (a field such as `"svType"` or `"impact"`, or a
#'   list of the color object's slots), `filter` (a list of any of
#'   `expressions`, jexl expressions, `minMinorAlleleFrequency` and
#'   `maxMissingness`),
#'   `referenceDrawingMode` (`"skip"` or `"draw"`), `shadeByDosage`, `rows`
#'   (the arrangement a session saves: `domain`, the rows listed first in
#'   order, `labels`, a label by row name, and `kept`, the rows shown) and
#'   `facet` (a samples file attribute, or `list(field, domain)`, whose values
#'   each take a band of rows) and `rowHeight` (pixels a row, where
#'   [jbrowse_stack()] otherwise fits the rows into the display's height)
#' @param index the tabix or CSI index, if not beside `file` under its usual
#'   name
#' @param samples a tab-separated file of the samples' attributes, a `name`
#'   column and one column per attribute, which `display$rowColor` colors
#'   the rows by, e.g. `list(rowColor = list(field = "group"))`
#' @param byte_limit the most compressed bytes the region may need, refused
#'   off the index before any record is fetched
#' @param assembly a [jbrowse_assembly()], an assembly as a JBrowse config's
#'   `assemblies` entry writes it, or its FASTA or 2bit alone: `region` may
#'   name a contig by any alias it lists, and the reader names the contig in
#'   its tables as `region` does
#' @param width the pixels the figure spans across the regions, at which
#'   records sharing a pixel fade as the display fades them
#' @return a `jbrowse_variants` list of data frames, beside the slots
#'   [jbrowse_result] lists:
#'   * `variants`: each record the display keeps: `variant` (its key, from
#'     1, which `cells` and `annotations` name it by), `start`, `end`, `id`,
#'     `name`, `ref`, `alt`, `type`, `color`, its lane color, `filter`, the
#'     FILTER column as the file writes it (`"PASS"`, or the failed filters
#'     joined by `;`, `NA` for `.`), and a column per INFO key the header
#'     declares or a record uses, named as the key (`AC`, `AN`, ...), or
#'     `info_<key>` where the key names a column above. A flag is logical,
#'     `FALSE` where the record lacks it; a `Number=1` key is integer,
#'     numeric or character as its `Type` says; any other `Number` (`A`,
#'     `R`, `G`, `.`, 2) is a list column of vectors, so `AC` at a site with
#'     two alts is `c(33L, 1L)`. A record lacking a key has `NA`. Over
#'     several regions a record is listed once in each region it falls in,
#'     its copies sharing `variant`. `opacity` is the fade the display paints
#'     the record's cells at, where records too short to see share a pixel at
#'     `width`
#'   * `annotations`: each record's SnpEff `ANN` or VEP `CSQ` consequences, a
#'     row per record and consequence: `variant`, `key` (`"ANN"` or `"CSQ"`) and a column per field the
#'     key's header line names (`Consequence`, `IMPACT`, `SYMBOL`, ...), each
#'     character, as the header types none. The raw string leaves the INFO
#'     columns
#'   * `cells`: each row's cell at a record: `variant`, `row` (the row's
#'     `name`, a factor in `rows`' order), `genotype`, `dosage` (the share of the called alleles that
#'     are an ALT; under `unit = "haplotype"` the row's allele, 1 for an ALT;
#'     `NA` for a no-call or a haplotype of an unphased call) and `color`, the
#'     fill the display paints; under `referenceDrawingMode = "skip"`, the
#'     default, a reference call has no cell. A cell takes its record's
#'     `start` and `end` and its row's `sample` when drawn. [as.matrix()]
#'     gives the genotypes as the display's clustering reads them
#'   * `rows`: `name`, `label`, `sample`, `haplotype`, `ploidy`, `band`,
#'     its `facet` value, and `color`, the bar the display marks the row with
#'     under `rowColor`, `""` where none. The rows come in the display's
#'     arrangement, and the matrix draws them in this table's order, and only
#'     these rows: reorder or filter it to rearrange them. Assigning `rows`
#'     drops the cells of the rows it leaves out, as assigning `variants`
#'     drops those of its records. Set `tree` to an
#'     [hclust()] result, or anything [as.hclust()] takes, labelled by
#'     `name`, such as `hclust(dist(as.matrix(v)), "average")`, and the rows
#'     take its leaf order and [autoplot()] draws it beside them. A list of
#'     trees is a forest, one per band, as the display clusters each band
#'     apart
#'   * `bands`: each band's `key` and `label`, in the order they stack. Banded
#'     rows stay in their band whatever order `rows` or `tree` gives them
#'
#'   `fields$info` declares each INFO column's `key`, `number`, `type` and
#'   `description` off the header, and the `column` it takes in `variants`;
#'   `fields$annotations` each consequence key's `key` and `field`, in the
#'   header's order. `keys` holds `genotypes`, the display's key, and `rows`,
#'   the row colors' key.
#' @export
#' @examples
#' vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "1:155,000,001-155,050,000")
#' head(v$cells)
jbrowse_variants <- function(file, region = NULL, display = NULL,
                             width = DEFAULT_WIDTH, assembly = NULL,
                             samples = NULL, index = NULL, byte_limit = 2e8) {
  adapter <- adapter_of(file, function() variants_adapter(file, index, samples))
  assembly <- assembly_config(assembly)
  regions <- reader_regions(adapter, region, assembly)
  args <- list(adapterConfig = adapter, regions = engine_regions(regions),
    byteLimit = byte_limit)
  args$assembly <- assembly
  if (isTRUE(attr(regions, "whole"))) args$rowLimit <- WHOLE_FILE_ROWS
  args$display <- display_config(display)
  t <- call_engine("variants", args)
  x <- in_one_call(variants_result(t), regions, t, width)
  x$variants$opacity <- density_alpha(x$variants, diff(xlim_of(x)) / width)
  x
}

variants_result <- function(t) {
  variants <- as_frame(t$variants)
  variants <- cbind(variant = variants$key + 1L,
    variants[names(variants) != "key"])
  variants$filter <- as.character(variants$filter)
  info <- typed_columns(t$info, t$infoFields, nrow(variants), names(variants),
    "info")
  variants <- cbind(variants, info$frame)
  info_fields <- info$fields
  annotations <- annotation_frame(t$annotations, "variant", variants$variant)
  rows <- as_frame(t$rows)
  bands <- as_frame(t$bands)[c("key", "label")]
  cells <- as_frame(t$cells)
  cells$variant <- variants$variant[cells$variant + 1L]
  cells$row <- factor(cells$row, levels = rows$name)
  genotypes <- as_frame(t$legend)
  attr(genotypes, "title") <- c(genotypes$section[genotypes$section != ""],
    "Genotypes")[1]
  structure(
    list(variants = variants, cells = cells, rows = rows, bands = bands,
      annotations = annotations, display = t$config,
      keys = list(genotypes = genotypes, rows = key_frame(
        unlist(t$rowLegend$label) %||% character(0),
        unlist(t$rowLegend$color) %||% character(0), t$rowLegendTitle)),
      ramps = ramp_frames(t$ramps), y_scale = NULL,
      fields = Filter(Negate(is.null), list(info = info_fields,
        annotations = annotation_fields(t$annotationFields)))),
    primary = "cells",
    class = c("jbrowse_variants", "jbrowse_result"))
}

# The columns R adds to a table after the engine's, which a file's own key
# never takes
ADDED_COLUMNS <- c("chrom", "region", "opacity", "feature", "subfeature",
  "variant")

# A VCF's INFO or FORMAT values as a data frame of `n` rows, a column per
# key typed by `info_column()`, each named as its key, or `<prefix>_<key>`
# where the key names one of `taken` or ADDED_COLUMNS; and the keys' `fields`
# with the `column` each took
typed_columns <- function(values, fields, n, taken, prefix) {
  fields <- as_frame(fields)
  fields$column <- ifelse(fields$key %in% c(taken, ADDED_COLUMNS),
    paste0(prefix, "_", fields$key), fields$key)
  frame <- data.frame(row.names = seq_len(n))
  for (i in seq_len(nrow(fields))) {
    frame[[fields$column[i]]] <- info_column(
      values[[fields$key[i]]] %||% rep(NA, n), fields$number[i],
      fields$type[i])
  }
  list(frame = frame, fields = fields)
}

# A VCF read as features' own columns, aligned to the features `keys` name:
# each record's `ref`, `alt`, `filter` and INFO in `records`, and each
# sample's FORMAT values in `samples`, `feature` its record's key; NULL for a
# file that is not a VCF
vcf_columns <- function(vcf, keys, taken) {
  if (is.null(vcf)) return(NULL)
  n <- length(keys)
  records <- data.frame(ref = as.character(vcf$records$ref),
    alt = as.character(vcf$records$alt),
    filter = as.character(vcf$records$filter), row.names = seq_len(n))
  info <- typed_columns(vcf$info, vcf$infoFields, n,
    c(taken, names(records)), "info")
  samples <- data.frame(feature = keys[as.integer(vcf$samples$feature) + 1L],
    sample = as.character(vcf$samples$sample))
  format <- typed_columns(vcf$format, vcf$formatFields, nrow(samples),
    c(names(samples), "chrom"), "format")
  list(records = cbind(records, info$frame), info_fields = info$fields,
    samples = cbind(samples, format$frame), format_fields = format$fields,
    annotations = annotation_frame(vcf$annotations, "feature", keys),
    annotation_fields = annotation_fields(vcf$annotationFields))
}

# Each record's ANN or CSQ consequences, a row each: `index`, the key of the
# record, from its row's in `keys`, `key`, the INFO key, and a column per
# field its header line names, each character, as the header types none
annotation_frame <- function(columns, index, keys) {
  d <- as_frame(columns)
  d[[index]] <- keys[as.integer(d[[index]]) + 1L]
  for (column in setdiff(names(d), c(index, "key"))) {
    d[[column]] <- as.character(d[[column]])
  }
  d
}

annotation_fields <- function(fields) {
  if (length(unlist(fields$key))) as_frame(fields)
}

# An INFO key's values as a column, typed by its header line: a flag as
# logical, a Number=1 key as a vector, any other Number as a list of vectors,
# one per record; NA where a record lacks the key. A key no record in the
# region gives arrives as NAs, which a list column takes one per record.
info_column <- function(x, number, type) {
  if (type == "Flag" || number == "0") return(as.logical(unlist(x)))
  cast <- switch(type, Integer = as.integer, Float = as.numeric,
    as.character)
  if (number == "1") return(cast(unlist(x)))
  lapply(as.list(x), cast)
}

#' @export
print.jbrowse_variants <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_variants> %s, %d records, %d rows, %d cells\n",
    print_loci(x), nrow(x$variants),
    nrow(x$rows), nrow(x$cells)))
  invisible(x)
}

#' A multi-sample VCF's genotype cells as a ggplot layer
#'
#' Each row's cell at each record, filled as the display paints it. At
#' `"genomic"` positions a cell spans its record, widened to two pixels of a
#' figure `width` pixels across, as the display widens a record narrower than
#' a cell; in `"columns"` each record is a column of the genotype matrix, x
#' being its row number in `variants`. Map `genotype`, `dosage`, `sample` or
#' any other column of the cells table to draw it your own way.
#'
#' Zoomed out, records too short to see share pixels, and the display fades
#' each by how many share its pixels, so a row shades by the share of its
#' records that are alt rather than filling wherever any is. The reader
#' gives each record that fade as its `opacity`, at the `width` it read at,
#' which the cells take at genomic positions; set `v$variants$opacity <- 1`
#' to draw every cell opaque.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_variants()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param layout `"genomic"` or `"columns"`
#' @param width the pixels the figure spans across its region at genomic
#'   positions; `NULL` takes the result's `width`
#' @param ... passed to [ggplot2::layer()] params
#' @return a ggplot layer
#' @export
#' @examples
#' vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "1:155,000,001-155,050,000")
#' library(ggplot2)
#' ggplot(v) + geom_jbrowse_genotypes(aes(fill = genotype), layout = "columns")
geom_jbrowse_genotypes <- function(mapping = NULL, data = NULL,
                                   layout = c("genomic", "columns"),
                                   width = NULL, ...) {
  data <- replacing(data, "cells", "geom_jbrowse_genotypes")
  layout <- match.arg(layout)
  cells <- function(x) {
    d <- x$cells
    d$alpha <- rep(1, nrow(d))
    if (layout == "columns") {
      column <- match(d$variant, unique(x$variants$variant))
      d$xmin <- column - 0.5
      d$xmax <- column + 0.5
    } else {
      floor_bp <- 2 / px_of(x, width)
      mid <- (d$start + d$end) / 2
      half <- pmax(d$end - d$start, floor_bp) / 2
      d$xmin <- mid - half
      d$xmax <- mid + half
      d$alpha <- d$opacity %||% d$alpha
    }
    d$y <- match(as.character(d$row), x$rows$name)
    d[!is.na(d$y), , drop = FALSE]
  }
  jbrowse_layer(ggplot2::GeomRect, data, cells, mapping,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - 0.5, ymax = y + 0.5,
      fill = I(color), alpha = I(alpha)), list(...))
}

# Each record's fade where records share pixels at `bp_per_px`, region by
# region, as the display fades each displayed region's cells
density_alpha <- function(variants, bp_per_px) {
  alpha <- rep(1, nrow(variants))
  region <- variants$region
  for (i in split(seq_len(nrow(variants)), region)) {
    if (length(i) < 2) next
    alpha[i] <- unlist(call_engine("variantDensityAlpha", list(
      start = variants$start[i], end = variants$end[i],
      bpPerPx = bp_per_px)))
  }
  alpha
}

# Every row's background in the reference call's grey, the key's first entry,
# as the display paints a row behind the cells it draws
reference_rows <- function(x) {
  background <- data.frame(xmin = -Inf, xmax = Inf, ymin = 0.5,
    ymax = nrow(x$rows) + 0.5, color = x$keys$genotypes$color[1])
  ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin,
    ymax = ymax, fill = I(color)), data = background, inherit.aes = FALSE)
}

#' Draw a multi-sample VCF as JBrowse's variant matrix
#'
#' [geom_jbrowse_genotypes()] with a row per sample (or haplotype) named on
#' the y axis, keyed by the display's genotype colors.
#' @param object a [jbrowse_variants()] result
#' @param layout `"genomic"` for cells at their records' positions,
#'   `"columns"` for the genotype matrix, a column per record
#' @param width the pixels the figure draws the regions across, which widen
#'   a record to two of them at genomic positions; `NULL` takes the `width`
#'   the variants were read at, which their fade keeps
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot
#' @exportS3Method ggplot2::autoplot
#' @examples
#' vcf <- system.file("extdata", "1000g_dcst2.vcf.gz", package = "ggjbrowse")
#' autoplot(jbrowse_variants(vcf, "1:155,000,001-155,050,000"))
autoplot.jbrowse_variants <- function(object,
                                      layout = c("genomic", "columns"),
                                      width = NULL, ...) {
  rlang::check_dots_empty()
  layout <- match.arg(layout)
  object <- on_axis(drawn_at(object, width))
  rows <- as.character(object$rows$label %||% object$rows$name)
  genomic <- layout == "genomic"
  xlim <- if (genomic) xlim_of(object)
    else c(0.5, length(unique(object$variants$variant)) + 0.5)
  inner <- if (any(object$rows$color != "")) 9 else 2.2
  tree <- object$tree
  strip <- if (is.null(tree)) inner else inner + TREE_PT + 3
  bands <- row_bands(object$rows, object$bands)
  outer <- if (is.null(bands)) strip else strip + BAND_PT + 2
  named <- length(rows) <= MAX_NAMED_ROWS
  ggplot2::ggplot(object) +
    reference_rows(object) +
    geom_jbrowse_genotypes(layout = layout) +
    guides_jbrowse(object) +
    row_bars(object, xlim) +
    (if (!is.null(tree)) geom_jbrowse_tree(gap = inner + 1, span = TREE_PT)) +
    band_labels(bands, strip + 1) +
    (if (length(rows)) ggplot2::scale_y_reverse(breaks = seq_along(rows),
      labels = if (named) rows else rep("", length(rows)))
    else ggplot2::scale_y_reverse(breaks = NULL)) +
    ggplot2::coord_cartesian(xlim = xlim, expand = FALSE, clip = "off") +
    (if (genomic) genomic_x(object)) +
    ggplot2::labs(x = if (genomic) axis_label(object)
      else "record", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    if (named) row_names(rows, seq_along(rows), gap_pt = outer)
    else if (outer > inner) {
      ggplot2::theme(axis.text.y = ggplot2::element_blank()) + axis_gap(outer)
    } else no_y
}

# The points a band's name takes beside the rows
BAND_PT <- 10

# Each band's name written up its rows, `gap` points left of the panel, as
# the display writes its bands in a strip beside the tree; none where its
# rows are too short for it
band_labels <- function(bands, gap) {
  if (is.null(bands)) return(NULL)
  ggplot2::layer(geom = GeomBandLabels, stat = "identity",
    position = "identity", data = data.frame(ymin = bands$first - 0.5,
      ymax = bands$last + 0.5, label = bands$label),
    mapping = ggplot2::aes(ymin = ymin, ymax = ymax, label = label),
    inherit.aes = FALSE, show.legend = FALSE, params = list(gap = gap))
}

GeomBandLabels <- ggplot2::ggproto("GeomBandLabels", ggplot2::Geom,
  required_aes = c("ymin", "ymax", "label"),
  default_aes = ggplot2::aes(colour = "grey20"),
  extra_params = c("na.rm", "gap"),
  draw_panel = function(data, panel_params, coord, gap) {
    d <- coord$transform(data, panel_params)
    grid::gTree(children = grid::gList(grid::textGrob(d$label,
      x = grid::unit(0, "npc") - grid::unit(gap + BAND_PT / 2, "pt"),
      y = grid::unit((d$ymin + d$ymax) / 2, "npc"), rot = 90,
      gp = grid::gpar(col = d$colour, fontsize = BAND_PT - 3))),
      span = abs(d$ymax - d$ymin), cl = "ggjbrowse_band_labels")
  })

#' @exportS3Method grid::makeContent
makeContent.ggjbrowse_band_labels <- function(x) {
  text <- x$children[[1]]
  room <- x$span * grid::convertHeight(grid::unit(1, "npc"), "pt",
    valueOnly = TRUE)
  fits <- vapply(seq_along(text$label), function(i) grid::convertWidth(
    grid::grobWidth(grid::textGrob(text$label[i], gp = text$gp)), "pt",
    valueOnly = TRUE) <= room[i] - 6, NA)
  text$label[!fits] <- ""
  grid::setChildren(x, grid::gList(text))
}

# The most rows a matrix names. Past it, rows carry only their group's color,
# as the display's rows do once they are too short for text.
MAX_NAMED_ROWS <- 150

# Each row's color as a bar beside its label, outside the matrix, as the
# display marks a row by its rowColor, keyed under the attribute
row_bars <- function(x, xlim) {
  colored <- x$rows$color != ""
  if (!any(colored)) return(NULL)
  bars <- data.frame(y = which(colored), color = x$rows$color[colored],
    x = xlim[1] - 0.006 * diff(xlim))
  list(ggplot2::geom_segment(ggplot2::aes(x = x, xend = x, y = y - 0.5,
      yend = y + 0.5, colour = I(color)), data = bars, inherit.aes = FALSE,
      linewidth = 2.5, lineend = "butt"))
}

# A matrix row's height in display pixels: the display's rowHeight, or for
# 0, the rows fitted to the display's height, as the display fits them
row_px <- function(x) {
  if (isTRUE(x$display$rowHeight > 0)) x$display$rowHeight
  else x$display$height / max(1, nrow(x$rows))
}
