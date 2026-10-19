# A feature file's adapter, read off its extension: GFF3, GTF and BED plain or
# bgzipped with a tabix (or CSI) index beside it, or BigBed.
features_adapter <- function(file, index) {
  entry <- function(pattern, type, slot, index = NULL) {
    list(pattern = pattern, type = type, slot = slot, index = index)
  }
  guess_adapter(file, index, list(
    entry("\\.(bb|bigbed)$", "BigBedAdapter", "bigBedLocation"),
    entry("\\.gff3?\\.gz$", "Gff3TabixAdapter", "gffGzLocation", "TBI"),
    entry("\\.gff3?$", "Gff3Adapter", "gffLocation"),
    entry("\\.gtf\\.gz$", "GtfTabixAdapter", "gtfGzLocation", "TBI"),
    entry("\\.gtf$", "GtfAdapter", "gtfLocation"),
    entry("\\.bed\\.gz$", "BedTabixAdapter", "bedGzLocation", "TBI"),
    entry("\\.vcf\\.gz$", "VcfTabixAdapter", "vcfGzLocation", "TBI"),
    entry("\\.vcf$", "VcfAdapter", "vcfLocation"),
    entry("\\.bed$", "BedAdapter", "bedLocation")),
    ".gff3, .gtf, .bed or .vcf (each optionally .gz with an index) or .bb")
}

#' Read a region of a feature file as JBrowse draws it
#'
#' Runs JBrowse's own feature display in V8: the features it fetches, the gene
#' glyph it builds from each (transcripts, exons, CDS, UTRs, intron lines,
#' strand chevrons and arrows), the rows its layout packs them into, the
#' labels it shows at this zoom and the colors its `color` object paints.
#' Coordinates are 0-based, half-open, on the reference; vertical positions
#' are the display's pixels, down from the track's top.
#'
#' @param file a GFF3, GTF, BED or VCF file (plain, or bgzipped with a tabix
#'   or CSI index beside it) or a BigBed, a path or a URL; or a JBrowse track's
#'   `adapter` config, as a list, for an adapter's own settings, or a
#'   `FromConfigAdapter` whose `features` list them inline
#' @param region the regions to read, read as one view and drawn end to end
#'   on one axis as a JBrowse view lays them out: loci as JBrowse's location
#'   box writes them (`"ctgA:1,001-9,000"`, 1-based and inclusive), contig
#'   names alone for whole contigs, a data frame of `chrom`, `start` and `end`
#'   (0-based, half-open, as BED writes them), or a GRanges. With `NULL`, the
#'   default, the reader reads the whole file: each contig it lists but those
#'   its index shows empty, to the length the file or `assembly` states. A whole
#'   file or contig is refused past 5 million rows in any one table
#' @param display the display's settings, as a track config's `displays`
#'   entry writes them for a `LinearBasicDisplay`: a list with any of
#'   `color` (a constant such as `"purple"`, or the color object's slots,
#'   such as `list(field = "type")`), `displayMode` (`"normal"`, `"compact"`,
#'   `"superCompact"`, `"collapsed"`), `showLabels` (`"auto"`, `"name"`,
#'   `"nameAndDescription"`, `"description"`, `"none"`), `geneGlyphMode`
#'   (`"auto"`, `"all"`, `"longestCoding"`), `subfeatureLabels` (`"none"`,
#'   `"below"`, `"overlay"`), `featureHeight`, `filter` (jexl
#'   expressions under `expressions`, `genes = "only"` or `"exclude"`),
#'   `showChevrons`, `facet` (a field such as
#'   `"strand"` or `"type"`, or `list(field = , domain = , hidden = )`, which
#'   packs a section per value, as alignments' `facet` does), and `heightMode`
#'   with `height`. A figure has no track height to fit, so `heightMode`
#'   defaults to `"grow"`, every feature at full size; `list(heightMode =
#'   "fit", height = 100)` gives up the labels, isoforms and row height the
#'   browser gives up to fit a 100 px track
#' @param width the width in pixels the display lays the region out across:
#'   whether a name fits beside its feature, which chevrons and arrows draw and
#'   how the `auto` modes resolve all follow from it. [autoplot()] scales those
#'   pixels to the panel, so a smaller `width` gives larger labels in a figure
#'   of the same size
#' @param assembly the reference: a [jbrowse_assembly()], an assembly as a
#'   JBrowse config's `assemblies` entry writes it, or an indexed FASTA or a
#'   2bit alone, path or URL. `color_by_cds` and `amino_acids` read its
#'   sequence. `region` may name a contig by any alias it lists, the engine
#'   reads the file under the file's own, and the reader names the contig in
#'   its tables as `region` does
#' @param color_by_cds paint each CDS by its reading frame, as the view's
#'   "Color by CDS" does
#' @param amino_acids with `assembly`, translate each CDS once the layout
#'   reaches a base a pixel, as the view's "Show amino acids" does: codon
#'   stripes in `rects`, and each residue in `amino_acids`
#' @param index the tabix or CSI index, if not beside `file` under its usual
#'   name
#' @param byte_limit the most compressed bytes the region may need, refused
#'   off the index before any feature is fetched
#' @return a `jbrowse_features` list of data frames, beside the slots
#'   [jbrowse_result] lists:
#'   * `features`: each top-level feature on its row: `feature` (its key,
#'     from 1, which the tables below name it by), `id` (the engine's id for
#'     it), `name`, `type`, `start`, `end`, `strand`, `top` and `bottom` (its
#'     row's pixels, labels included), `gene` and `section`, then a column
#'     per field the file states (see below)
#'   * `subfeatures`: a gene's transcripts: `subfeature` (its key, from 1),
#'     `id`, `name`, `type`, `start`, `end`, `top`, `bottom`, `feature`, and
#'     the file's fields as `features` has them
#'   * `rects`: the boxes the glyph draws (exons, CDS, UTRs, a plain feature):
#'     `start`, `end`, `top`, `height`, `color`, `feature`, `subfeature`
#'   * `lines`: the intron and connector lines: `start`, `end`, `y`, `color`,
#'     `feature`, `subfeature`
#'   * `chevrons`: the strand chevrons on those lines where the display draws
#'     them: `position`, `y`, `direction`, `color`, `feature`, `subfeature`
#'   * `arrows`: the strand arrow off a feature's end where it is wide enough
#'     on screen: `position`, `y`, `direction`, `head_half_height`, `color`,
#'     `feature`, `subfeature`
#'   * `labels`: each name, description and transcript label the display
#'     shows: `text`, `kind`, `position` (its left edge, in bp), `y` (its top),
#'     `color`, `feature`, and `subfeature`, the transcript a transcript label
#'     names (`NA` for a feature's own label)
#'   * `sections`: under `display$facet`, each section's `key`, `label`, `top`
#'     and `height` in stacking order, the chip row the display names it in
#'     included; ungrouped, one section labelled `""`
#'   * `amino_acids`: each codon over a CDS box in the region, a codon split
#'     by an intron in a piece on each side: `start`, `end`, `amino_acid`,
#'     `protein_index` (0-based), `top`, `height`, `stop_or_non_triplet`,
#'     `transl_except`, `feature`, and the letter the display writes over it,
#'     `text` (empty under 8 pixels a base or in a box too short), `font_size`
#'     and `y` (its baseline) in pixels, and `color`
#'
#'   **The file's fields.** `features` and `subfeatures` carry a column per
#'   field the file states, keyed as the file writes it. A GFF3 gives its
#'   column-9 attributes (`ID`, `Name`, `Parent`, `gene_name`, ...) as
#'   character, beside `source`, `score` and `phase`. `id` and `name` are the
#'   engine's own; a GFF3's `ID` and `Name` sit beside them. A key naming one
#'   of the table's own columns takes `attr_<key>`, as an INFO key takes
#'   `info_<key>`. A key the GFF3 spec lets list several values (`Parent`,
#'   `Alias`, `Note`, `Dbxref`, `Ontology_term`), or that a record lists
#'   several of, is a list column. A BED or BigBed field takes the type its
#'   autoSql declares. A GFF3 has no header, so its columns are the keys the
#'   region's records use, and a region lacking one reads `NA` there when
#'   regions bind. `fields$attributes` declares each field: its `key`,
#'   `type`, `number` (`"1"`, or `"."` for a list column), `description` and
#'   the `column` it takes. `keys$features` is the key the display shows for
#'   a color field.
#'
#'   A VCF reads as JBrowse's plain variant track draws it, every record
#'   whether or not it calls a genotype, so a calls VCF with no GT reads
#'   here where [jbrowse_variants()] keeps none. `features` then also holds
#'   each record's `ref`, `alt`, `filter` and a column per INFO key, typed as
#'   [jbrowse_variants()] types them, and the result adds:
#'   * `samples`: a row per record and sample, `feature`, `sample` and a
#'     column per FORMAT key, typed as INFO is (`AD`, `DP`, ...), named
#'     `format_<key>` where the key names a column above
#'   * `annotations`: its SnpEff `ANN` or VEP `CSQ` consequences, as
#'     [jbrowse_variants()] gives them, keyed by `feature`
#'
#'   `fields$info` and `fields$format` declare each INFO and FORMAT column's
#'   `key`, `number`, `type`, `description` and `column`, and
#'   `fields$annotations` each consequence key's fields.
#'
#'   Each table below `features` names its feature by `feature`, and a
#'   gene's mark its transcript by `subfeature` (`NA` for a feature's own
#'   mark), both keys the same in every region a record is listed in; a
#'   mark's `region` picks the copy. A mark takes its feature's `name`,
#'   `type`, `strand` and `section` (a factor in stacking order) when drawn,
#'   and a chevron, arrow or label its `start` and `end` too. Assigning
#'   `features` or `subfeatures` drops the rows that name a record it leaves
#'   out. `glyph` gives the pixel
#'   sizes the display draws text, arrows, chevrons and section chips
#'   (`sectionLabelHeight`) at, and `height` the track's laid-out height in
#'   pixels.
#' @export
#' @examples
#' gff <- system.file("extdata", "gapdh_gencode.gff3.gz", package = "ggjbrowse")
#' g <- jbrowse_features(gff, "chr12:6,534,201-6,538,700")
#' head(g$features)
jbrowse_features <- function(file, region = NULL, display = NULL,
                             width = DEFAULT_WIDTH, assembly = NULL,
                             color_by_cds = FALSE, amino_acids = TRUE,
                             index = NULL, byte_limit = 2e8) {
  adapter <- adapter_of(file, function() features_adapter(file, index))
  assembly <- assembly_config(assembly)
  regions <- reader_regions(adapter, region, assembly)
  args <- list(
    adapterConfig = adapter,
    regions = engine_regions(regions),
    bpPerPx = sum(regions$end - regions$start) / width,
    colorByCDS = color_by_cds,
    aminoAcids = amino_acids,
    byteLimit = byte_limit)
  args$assembly <- assembly
  display <- display %||% list()
  display$heightMode <- display$heightMode %||% "grow"
  args$display <- display_config(display)
  if (isTRUE(attr(regions, "whole"))) args$rowLimit <- WHOLE_FILE_ROWS
  t <- call_engine("features", args)
  key <- structure(as_frame(t$legend), title = t$legendTitle)
  features <- as_frame(t$features)
  features <- cbind(feature = match(features$id, unique(features$id)),
    features)
  subfeatures <- as_frame(t$subfeatures)
  subfeature_key <- match(subfeatures$id, unique(subfeatures$id))
  sections <- as_frame(t$sections)
  vcf <- vcf_columns(t$vcf, features$feature, names(features))
  taken <- c(names(features), names(subfeatures), "strand", "section")
  attributes <- typed_columns(t$attributes, t$attributeFields,
    nrow(features), taken, "attr")
  subfeature_attributes <- typed_columns(t$subfeatureAttributes,
    t$attributeFields, nrow(subfeatures), taken, "attr")$frame
  features$section <- factor(features$section, levels = sections$label)
  features <- cbind(features, attributes$frame)
  if (!is.null(vcf)) features <- cbind(features, vcf$records)
  on_feature <- function(table) {
    d <- as_frame(table)
    d$feature <- features$feature[d$feature + 1L]
    if (!is.null(d$subfeature)) {
      at <- as.integer(d$subfeature) + 1L
      at[at < 1L] <- NA
      d$subfeature <- subfeature_key[at]
    }
    d
  }
  subfeatures <- cbind(subfeature = subfeature_key,
    on_feature(t$subfeatures), subfeature_attributes)
  arrows <- on_feature(t$arrows)
  names(arrows)[names(arrows) == "headHalfHeight"] <- "head_half_height"
  tables <- list(features = features, subfeatures = subfeatures,
    rects = on_feature(t$rects), lines = on_feature(t$lines),
    chevrons = on_feature(t$chevrons), arrows = arrows,
    amino_acids = snake(on_feature(t$aminoAcids)),
    labels = on_feature(t$labels), sections = sections,
    samples = vcf$samples, annotations = vcf$annotations)
  fields <- list(
    attributes = if (nrow(attributes$fields)) attributes$fields,
    info = vcf$info_fields, format = vcf$format_fields,
    annotations = vcf$annotation_fields)
  in_one_call(structure(
    c(Filter(Negate(is.null), tables),
      list(display = t$config, keys = list(features = key),
        ramps = ramp_frames(t$ramps), y_scale = NULL,
        fields = Filter(Negate(is.null), fields), glyph = t$glyph,
        height = t$height, bp_per_px = t$bpPerPx)),
    primary = "features",
    class = c("jbrowse_features", "jbrowse_result")),
    regions, t, width)
}

# Text sized in the display's pixels, scaled to the panel when it draws:
# `px_wide` and `px_tall` are the display's extent in pixels, and a `halo`
# color strokes each letter a display pixel wide, as canvas `strokeText` does
GeomDisplayText <- ggplot2::ggproto("GeomDisplayText", ggplot2::GeomText,
  default_aes = ggplot2::aes(!!!ggplot2::GeomText$default_aes,
    px_wide = NA, px_tall = NA, halo = NA),
  draw_panel = function(self, data, panel_params, coord, ...) {
    text <- ggplot2::ggproto_parent(ggplot2::GeomText, self)$draw_panel(data,
      panel_params, coord, ...)
    grid::gTree(children = grid::gList(text), px_wide = data$px_wide[1],
      px_tall = data$px_tall[1], halo = data$halo[1],
      cl = "ggjbrowse_display_text")
  })

#' @exportS3Method grid::makeContent
makeContent.ggjbrowse_display_text <- function(x) {
  mm_per_px <- min(
    grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE) /
      x$px_wide,
    grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE) /
      x$px_tall)
  text <- x$children[[1]]
  text$gp$fontsize <- text$gp$fontsize * mm_per_px
  if (is.na(x$halo)) return(grid::setChildren(x, grid::gList(text)))
  angle <- seq(0, 2 * pi, length.out = 9)[-9]
  r <- grid::unit(mm_per_px / 2, "mm")
  halos <- lapply(seq_along(angle), function(i) {
    grid::editGrob(text, x = text$x + r * cos(angle[i]),
      y = text$y + r * sin(angle[i]), gp = grid::gpar(col = x$halo),
      name = paste0(text$name, ".halo", i))
  })
  grid::setChildren(x, do.call(grid::gList, c(halos, list(text))))
}

#' @export
print.jbrowse_features <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_features> %s, %d features, %d transcripts, %d boxes, %d labels, %.0f px tall\n",
    print_loci(x), nrow(x$features),
    nrow(x$subfeatures), nrow(x$rects), nrow(x$labels), x$height))
  invisible(x)
}

# A chevron's two arms, each a segment, at the display's px sizes
chevron_arms <- function(x) {
  chev <- x$chevrons
  half_w <- x$glyph$chevronWidth / 2 * x$bp_per_px
  half_h <- x$glyph$chevronHeight / 2
  rbind(
    transform(chev, x = position - half_w * direction, y0 = y - half_h,
      xend = position + half_w * direction, yend = y),
    transform(chev, x = position + half_w * direction, y0 = y,
      xend = position - half_w * direction, yend = y + half_h))
}

# The glyph's boxes, each at least a display pixel wide, as the display draws
# a 1 bp feature at any zoom
boxes <- function(x) {
  d <- x$rects
  d$end <- pmax(d$end, d$start + x$bp_per_px)
  d
}

# A strand arrow's stem, as a box in data units
arrow_stems <- function(x) {
  a <- x$arrows
  stem_end <- a$position + x$glyph$arrowLength / 2 * x$bp_per_px * a$direction
  half <- x$glyph$arrowStemHalfHeight
  transform(a, xmin = pmin(position, stem_end), xmax = pmax(position, stem_end),
    ymin = y - half, ymax = y + half)
}

# A strand arrow's head, as a triangle's three corners
arrow_heads <- function(x) {
  a <- x$arrows
  length <- x$glyph$arrowLength * x$bp_per_px
  rows <- rep(seq_len(nrow(a)), each = 3)
  h <- a[rows, , drop = FALSE]
  corner <- rep(1:3, nrow(a))
  h$id <- rows
  h$x <- h$position + length * h$direction * ifelse(corner == 3, 1, 0.5)
  h$y <- h$y + c(-1, 1, 0)[corner] * h$head_half_height
  h
}

# The entries of `mapping` a geom draws
for_geom <- function(mapping, geom) {
  if (is.null(mapping)) return(NULL)
  mapping[intersect(names(mapping), geom$aesthetics())]
}

#' A feature track's marks as ggplot layers
#'
#' `geom_jbrowse_glyphs()` draws each feature's glyph as the display does:
#' the intron and connector lines with their strand chevrons, the boxes, and
#' the strand arrows. `geom_jbrowse_amino_acids()` writes the residues the
#' display letters over a zoomed-in CDS, in a monospace face haloed in
#' white, red at a stop or a partial codon. `geom_jbrowse_labels()` writes
#' the names and descriptions the display shows at this zoom. Both take the
#' display's text size scaled to the panel. Each mark is drawn with its
#' feature's `name`, `type` and `strand` beside the display's `color`, and a
#' mark with no span of its own (a chevron, an arrow, a label) its feature's
#' `start` and `end`, so `aes(fill = type)` colors the boxes by feature type
#' through ggplot's own scales, beside the display's colors on the rest. A
#' `fill` mapping reaches the boxes and arrowheads, and a `colour` mapping
#' every part that has a stroke, the boxes' outlines included.
#' @param mapping aesthetics that override the defaults
#' @param data a [jbrowse_features()] result, `NULL` for the plot's,
#'   or a data frame in place of the table the geom draws from the plot's
#'   result, such as a filtered copy of it
#' @param ... passed to [ggplot2::layer()] params
#' @return a list of ggplot layers
#' @export
#' @examples
#' gff <- system.file("extdata", "gapdh_gencode.gff3.gz", package = "ggjbrowse")
#' g <- jbrowse_features(gff, "chr12:6,534,201-6,538,700")
#' library(ggplot2)
#' ggplot(g) + geom_jbrowse_glyphs(aes(fill = type)) + scale_y_reverse()
geom_jbrowse_glyphs <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, NULL, "geom_jbrowse_glyphs")
  list(
    jbrowse_layer(ggplot2::GeomSegment, data, function(x) x$lines,
      for_geom(mapping, ggplot2::GeomSegment),
      ggplot2::aes(x = start, xend = end, y = y, yend = y, colour = I(color)),
      list(linewidth = 0.3, ...)),
    jbrowse_layer(ggplot2::GeomSegment, data, chevron_arms,
      for_geom(mapping, ggplot2::GeomSegment),
      ggplot2::aes(x = x, xend = xend, y = y0, yend = yend, colour = I(color)),
      list(linewidth = 0.3, ...)),
    jbrowse_layer(ggplot2::GeomRect, data, boxes,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = start, xmax = end, ymin = top, ymax = top + height,
        fill = I(color)), list(...)),
    jbrowse_layer(ggplot2::GeomRect, data, arrow_stems,
      for_geom(mapping, ggplot2::GeomRect),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = I(color)), list(...)),
    jbrowse_layer(ggplot2::GeomPolygon, data, arrow_heads,
      for_geom(mapping, ggplot2::GeomPolygon),
      ggplot2::aes(x = x, y = y, group = id, fill = I(color)), list(...)))
}

#' @rdname geom_jbrowse_glyphs
#' @export
geom_jbrowse_amino_acids <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "amino_acids", "geom_jbrowse_amino_acids")
  lettered <- function(x) {
    d <- x$amino_acids[nzchar(x$amino_acids$text), , drop = FALSE]
    n <- nrow(d)
    d$position <- (d$start + d$end) / 2
    d$px_wide <- rep(diff(xlim_of(x)) / x$bp_per_px, n)
    d$px_tall <- if (is.null(d$section)) rep(x$height, n)
      else section_px(x)[as.integer(d$section)]
    d$halo <- rep(x$glyph$aminoAcidHalo, n)
    d
  }
  jbrowse_layer(GeomDisplayText, data, lettered, mapping,
    ggplot2::aes(x = position, y = y, label = text, colour = I(color),
      size = I(font_size), px_wide = I(px_wide), px_tall = I(px_tall),
      halo = I(halo)),
    list(hjust = 0.5, vjust = 0, family = "mono", ...))
}

#' @rdname geom_jbrowse_glyphs
#' @export
geom_jbrowse_labels <- function(mapping = NULL, data = NULL, ...) {
  data <- replacing(data, "labels", "geom_jbrowse_labels")
  labels <- function(x) {
    d <- x$labels
    n <- nrow(d)
    d$size <- rep(x$glyph$fontSize, n)
    d$px_wide <- rep(diff(xlim_of(x)) / x$bp_per_px, n)
    d$px_tall <- if (is.null(d$section)) rep(x$height, n)
      else section_px(x)[as.integer(d$section)]
    d
  }
  jbrowse_layer(GeomDisplayText, data, labels, mapping,
    ggplot2::aes(x = position, y = y, label = text, colour = I(color),
      size = I(size), px_wide = I(px_wide), px_tall = I(px_tall)),
    list(hjust = 0, vjust = 1, ...))
}

#' Draw features as JBrowse's feature track
#'
#' Each feature's glyph on its row ([geom_jbrowse_glyphs()]), the residues
#' over a zoomed-in CDS ([geom_jbrowse_amino_acids()]) and the labels
#' the display shows ([geom_jbrowse_labels()]), keyed by the display's
#' `color` object. The y axis is the display's pixels. Labels scale with the
#' panel: at draw time they take the display's pixel size across the panel's
#' width or down its height, whichever is smaller, so they keep inside the
#' room the layout left them at any output size. Under the display's `facet`,
#' each section draws in its own panel over its rows, named by its strip.
#'
#' @param object a [jbrowse_features()] result
#' @param width the pixels the figure draws the regions across, which
#'   [jbrowse_width()] sizes a stack by; `NULL` takes the `width` the track
#'   was read at. Narrower than that warns, since the labels were laid out
#'   for more pixels
#' @param ... must be empty, so a misspelt argument errors
#' @return a ggplot, to extend or to stack over an alignments track with
#'   patchwork
#' @exportS3Method ggplot2::autoplot
#' @examples
#' gff <- system.file("extdata", "gapdh_gencode.gff3.gz", package = "ggjbrowse")
#' g <- jbrowse_features(gff, "chr12:6,534,201-6,538,700")
#' autoplot(g)
autoplot.jbrowse_features <- function(object, width = NULL, ...) {
  rlang::check_dots_empty()
  object <- on_axis(drawn_at(object, width))
  faceted <- !identical(object$sections$label %||% "", "")
  p <- ggplot2::ggplot(object)
  if (faceted) p <- p + section_extents(object)
  p +
    geom_jbrowse_glyphs() +
    geom_jbrowse_amino_acids() +
    guides_jbrowse(object) +
    geom_jbrowse_labels() +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(
      xlim = xlim_of(object),
      ylim = if (!faceted) c(0, object$height), expand = FALSE) +
    by_section(object, rows = TRUE) +
    genomic_x(object) +
    ggplot2::labs(x = axis_label(object), y = NULL) +
    ggplot2::theme_minimal() +
    across +
    no_y
}

# Each section's rows in display pixels, below the chip row a facet strip
# stands in for
section_rows <- function(x) {
  s <- x$sections
  chip <- x$glyph$sectionLabelHeight %||% 0
  data.frame(section = factor(s$label, levels = s$label),
    top = s$top + chip, bottom = s$top + s$height)
}

section_px <- function(x) {
  r <- section_rows(x)
  r$bottom - r$top
}

# Each panel spans its section's rows, whatever marks it holds
section_extents <- function(x) {
  r <- section_rows(x)
  ggplot2::geom_blank(ggplot2::aes(y = y),
    data = data.frame(section = rep(r$section, 2), y = c(r$top, r$bottom)),
    inherit.aes = FALSE)
}
