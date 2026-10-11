variants_adapter <- function(file, index) {
  if (grepl("\\.vcf$", sub("\\?.*$", "", tolower(file)))) {
    return(list(type = "VcfAdapter", vcfLocation = location(file)))
  }
  index <- index %||% paste0(file, ".tbi")
  list(type = "VcfTabixAdapter", vcfGzLocation = location(file),
    index = list(location = location(index),
      indexType = if (grepl("\\.csi$", index)) "CSI" else "TBI"))
}

#' Read a region of a multi-sample VCF as JBrowse paints it
#'
#' Runs JBrowse's own multi-sample variant display in V8: the records its
#' filters keep, each sample's genotype at each, and the fill the display
#' paints each cell, its hue shaded by the alt dosage. Coordinates are 0-based,
#' half-open, on the reference.
#'
#' @param file a VCF, bgzipped with a tabix or CSI index beside it, or plain;
#'   a path or a URL; or a JBrowse track's `adapter` config, as a list
#' @param chrom,start,end the region, 0-based half-open
#' @param display the display's settings, as a track config's `displays`
#'   entry writes them for a `LinearMultiSampleVariantDisplay`: a list with
#'   any of `unit` (`"allele"` for a row per sample, `"haplotype"` for a row
#'   per haplotype), `color` (a field such as `"svType"` or `"impact"`, or a
#'   list of the color object's slots), `filter` (jexl expressions),
#'   `minorAlleleFrequencyFilter`, `maxMissingnessFilter`,
#'   `referenceDrawingMode` (`"skip"` or `"draw"`) and `shadeByDosage`
#' @param index the tabix or CSI index, if not beside `file` under its usual
#'   name
#' @param byte_limit the most compressed bytes the region may need, refused
#'   off the index before any record is fetched
#' @return a `jbrowse_variants` list of data frames:
#'   * `variants`: each record the display keeps: `start`, `end`, `id`,
#'     `name`, `ref`, `alt`, `type` and `color`, its lane color
#'   * `cells`: each row's cell at a record: `variant` (a row number of
#'     `variants`), `start`, `end`, `row` (a factor in the display's row
#'     order), `sample`, `genotype`, `dosage` (the alt share of the called
#'     alleles) and `color`, the fill the display paints; under
#'     `referenceDrawingMode = "skip"`, the default, a reference call has no
#'     cell
#'   * `rows`: `name`, `sample`, `haplotype` and `ploidy`
#'   * `legend`: the display's key, its title as the frame's `"title"`
#'     attribute
#' @export
#' @examples
#' vcf <- system.file("extdata", "dog10k_cyp1a2_snvs.vcf.gz",
#'   package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "chr30", 38258000, 38265000)
#' head(v$cells)
jbrowse_variants <- function(file, chrom, start, end, display = NULL,
                             index = NULL, byte_limit = 2e8) {
  args <- list(
    adapterConfig = adapter_of(file, function()
      variants_adapter(file, index)),
    region = list(refName = chrom, start = start, end = end),
    byteLimit = byte_limit)
  if (!is.null(display)) {
    display$color <- config_color(display$color)
    if (!is.null(display$filter)) display$filter <- I(display$filter)
    args$display <- display
  }
  t <- call_engine("variants", args)
  variants <- as_frame(t$variants)
  variants$color <- jbrowse_color(variants$color)
  rows <- as_frame(t$rows)
  cells <- as_frame(t$cells)
  cells$variant <- cells$variant + 1
  cells$start <- variants$start[cells$variant]
  cells$end <- variants$end[cells$variant]
  cells$sample <- rows$sample[match(cells$row, rows$name)]
  cells$row <- factor(cells$row, levels = rows$name)
  cells$color <- jbrowse_color(cells$color)
  legend <- as_frame(t$legend)
  legend$color <- jbrowse_color(legend$color)
  attr(legend, "title") <- c(legend$section[legend$section != ""],
    "Genotypes")[1]
  structure(
    list(variants = variants, cells = cells, rows = rows, legend = legend,
      keys = list(genotypes = legend),
      region = list(chrom = chrom, start = start, end = end)),
    primary = "cells",
    class = c("jbrowse_variants", "jbrowse_result"))
}

#' @export
print.jbrowse_variants <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_variants> %s:%.0f-%.0f, %d records, %d rows, %d cells\n",
    x$region$chrom, x$region$start, x$region$end, nrow(x$variants),
    nrow(x$rows), nrow(x$cells)))
  invisible(x)
}

#' A multi-sample VCF's genotype cells as a ggplot layer
#'
#' Each row's cell at each record, filled as the display paints it. At
#' `"genomic"` positions a cell spans its record, widened to two pixels of a
#' figure `px_per_bp` describes, as the display widens a record narrower than
#' a cell; in `"columns"` each record is a column of the genotype matrix, x
#' being its row number in `variants`. Map `genotype`, `dosage`, `sample` or
#' any other column of the cells table to draw it your own way.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_variants()] result, or `NULL` for the plot's
#' @param layout `"genomic"` or `"columns"`
#' @param px_per_bp how many pixels a base spans at genomic positions; `NULL`
#'   takes a figure about 1200 pixels wide
#' @param ... passed to [ggplot2::layer()] params
#' @return a ggplot layer
#' @export
#' @examples
#' vcf <- system.file("extdata", "dog10k_cyp1a2_snvs.vcf.gz",
#'   package = "ggjbrowse")
#' v <- jbrowse_variants(vcf, "chr30", 38258000, 38265000)
#' library(ggplot2)
#' ggplot(v) + geom_jbrowse_genotypes(aes(fill = genotype), layout = "columns")
geom_jbrowse_genotypes <- function(mapping = NULL, data = NULL,
                                   layout = c("genomic", "columns"),
                                   px_per_bp = NULL, ...) {
  layout <- match.arg(layout)
  cells <- function(x) {
    d <- x$cells
    if (layout == "columns") {
      d$xmin <- d$variant - 0.5
      d$xmax <- d$variant + 0.5
    } else {
      floor_bp <- 2 / px_of(x, px_per_bp)
      mid <- (d$start + d$end) / 2
      half <- pmax(d$end - d$start, floor_bp) / 2
      d$xmin <- mid - half
      d$xmax <- mid + half
    }
    d$y <- as.integer(d$row)
    d
  }
  jbrowse_layer(ggplot2::GeomRect, data, cells, mapping,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - 0.5, ymax = y + 0.5,
      fill = color), list(...))
}

#' Draw a multi-sample VCF as JBrowse's variant matrix
#'
#' [geom_jbrowse_genotypes()] with a row per sample (or haplotype) named on
#' the y axis, keyed by the display's genotype colors.
#' @param object a [jbrowse_variants()] result
#' @param layout `"genomic"` for cells at their records' positions,
#'   `"columns"` for the genotype matrix, a column per record
#' @param ... unused
#' @return a ggplot
#' @exportS3Method ggplot2::autoplot
#' @examples
#' vcf <- system.file("extdata", "dog10k_cyp1a2_snvs.vcf.gz",
#'   package = "ggjbrowse")
#' autoplot(jbrowse_variants(vcf, "chr30", 38258000, 38265000))
autoplot.jbrowse_variants <- function(object,
                                      layout = c("genomic", "columns"), ...) {
  layout <- match.arg(layout)
  rows <- object$rows$name
  genomic <- layout == "genomic"
  ggplot2::ggplot(object) +
    geom_jbrowse_genotypes(layout = layout) +
    scale_fill_jbrowse(object) +
    ggplot2::scale_y_reverse(breaks = seq_along(rows), labels = rows) +
    ggplot2::coord_cartesian(
      xlim = if (genomic) c(object$region$start, object$region$end)
        else c(0.5, nrow(object$variants) + 0.5),
      expand = FALSE) +
    ggplot2::labs(x = if (genomic) sprintf("%s (bp)", object$region$chrom)
      else "record", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6))
}
