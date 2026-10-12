alignments_adapter <- function(file, index) {
  if (grepl("\\.cram$", file, ignore.case = TRUE)) {
    list(type = "CramAdapter", cramLocation = location(file),
      craiLocation = location(index %||% paste0(file, ".crai")))
  } else {
    index <- index %||% paste0(file, ".bai")
    list(type = "BamAdapter", bamLocation = location(file),
      index = list(location = location(index),
        indexType = if (grepl("\\.csi$", index)) "CSI" else "BAI"))
  }
}

sequence_adapter <- function(reference) {
  if (grepl("\\.2bit$", reference, ignore.case = TRUE)) {
    list(type = "TwoBitAdapter", twoBitLocation = location(reference))
  } else if (grepl("\\.gz$", reference)) {
    list(type = "BgzipFastaAdapter", fastaLocation = location(reference),
      faiLocation = location(paste0(reference, ".fai")),
      gziLocation = location(paste0(reference, ".gzi")))
  } else {
    list(type = "IndexedFastaAdapter", fastaLocation = location(reference),
      faiLocation = location(paste0(reference, ".fai")))
  }
}

#' Read a region of a BAM or CRAM file as JBrowse draws it
#'
#' Runs JBrowse's own alignments engine in V8: the reads it fetches, the
#' pileup rows its display packs them into, their mismatches against the
#' reference, their deletions and splice gaps, and the coverage. Coordinates
#' are 0-based, half-open, on the reference.
#'
#' @param file a BAM or CRAM file, a path or a URL; or a JBrowse track's
#'   `adapter` config, as a list
#' @param region the region: a locus as JBrowse's location box writes it
#'   (`"ctgA:1,001-9,000"`, 1-based and inclusive), several of them, which
#'   are read in turn and drawn end to end on one axis as a JBrowse view lays
#'   them out, a data frame of `chrom`, `start` and `end`, or a contig name
#'   with `start` and `end`
#' @param start,end with a contig name as `region`, the region's bounds,
#'   0-based half-open
#' @param reference the reference: an indexed FASTA, plain or bgzipped, or a
#'   2bit, path or URL; a CRAM needs it, and a BAM read without an MD tag is
#'   compared against it
#' @param reference_chrom the contig's name in `reference`, where the file
#'   spells it another way (`1` in the BAM, `chr1` in the FASTA)
#' @param index the BAI, CSI or CRAI, if not beside `file` under its usual name
#' @param display the display's config, as a track's `displays` entry
#'   writes it for a `LinearAlignmentsDisplay`, which is what JBrowse's
#'   "Export R script" passes. Its main slots:
#'   * `color`: what fills each read: `list(field = "strand")` (or `"mapq"`,
#'     `"insertSize"`, `"pairOrientation"`, `"firstOfPairStrand"`,
#'     `"tags.XS"`), with optional `domain` and `range`, e.g.
#'     `list(field = "insertSize", domain = c(200, 600))`; a bare string is
#'     one constant color, as in a track config; unset fills every read plain
#'   * `filter`: which reads are kept: any of `flagInclude`, `flagExclude`
#'     (1540 by default), `readName`, `tagFilters` (a list of
#'     `list(tag = , value = )`), `spliced`, `properPairs`, `singletons` and
#'     `split` (each `"only"` or `"exclude"`)
#'   * `baseColor`: the per-base layer painted over the reads:
#'     `list(field = "modifications")` for the reads' MM/ML base
#'     modifications, `list(field = "baseQuality")` for each base's Phred
#'     score; unset paints none. Under `"modifications"` or `"bisulfite"` the
#'     display greys every mismatched base so the calls stand out
#'   * `modifications`: what that layer reads: any of `threshold` (the
#'     percent probability under which a call is hidden, 10), `unmodified`
#'     (`"hidden"`, `"calls"` or `"all"`), `cytosineContext` (`"CG"`,
#'     `"CHG"`, `"CHH"` or `"all"`) and `shownModifications` (the type codes
#'     drawn, every one when empty)
#'   * `sort`: how the rows are laid out: `list(type = "position")` (or
#'     `"length"`, `"spliced"`, `"split"`), or a column sort, `type`
#'     (`"strand"`, `"basePair"`, `"tag"`, `"insertion"`, `"softclip"`,
#'     `"hardclip"`) with `pos`, the 0-based column ranked over, and `tag`
#'   * `facet`: what splits the reads into stacked sections, each laid out on
#'     its own: `list(field = "strand")` (or `"pairOrientation"`,
#'     `"tags.HP"`, ...), with optional `domain` (section keys in stacking
#'     order) and `hidden` (keys left out)
#'   * `readConnections`: the read-connection band, `"arc"` for each mate
#'     pair's and split read's arc or `"cloud"` for the read cloud; unset
#'     reads none
#'   * `showCoverage`, `showPileup`, `showSashimiArcs`, and any other slot
#'     such as `showOrdinaryPairs`, `arcColor` or `minSashimiScore`
#' @param byte_limit the most compressed bytes the region may need, refused off
#'   the index before any read is fetched; the whole table is held in memory
#' @return a `jbrowse_alignments` list of data frames:
#'   * `reads`: `read` (the id the other tables' `read` column names it by,
#'     kept through filtering), `start`, `end`, `row`, `strand`, `mapq`,
#'     `flags`, `name`,
#'     `insert_size`, `pair_orientation`, `color` (the fill the display paints,
#'     as `#rrggbb`) and `category` (the bucket the display files the read
#'     under: `fwdStrand`, `mapq`, `unmappedMate`, ...)
#'   * `mismatches`: `position`, `base`, `read`, `row`, `frequency` (how many
#'     reads carry that base there, over the depth), `quality`, `fade`
#'   * `gaps`: `start`, `end`, `type` (`"deletion"` or `"skip"`), `read`, `row`,
#'     `fade`
#'   * `insertions`: `position`, `length`, `sequence`, `read`, `row`, `fade`
#'   * `clips`: `position`, `length`, `type` (`"soft"` or `"hard"`), `read`,
#'     `row`, `fade`
#'   * `modifications`: each base modification call the layer paints:
#'     `position`, `read`, `row`, `type` (the code, `m` for 5mC), `name`,
#'     `probability`, `unmodified` (whether the call is that the base is not
#'     modified), and the `color` and `opacity` the display paints it at;
#'     empty without `display$baseColor`
#'   * `qualities`: under `baseColor` `"baseQuality"`, each aligned base:
#'     `position`, `read`, `row`, `score` (Phred) and the `color` the display
#'     paints it
#'   * `coverage`: `start`, `end`, `depth`, `forward`, `reverse`
#'   * `allele_coverage`: the coverage band's stack at each mismatch position,
#'     as the display stacks it: `position`, `base`, `from` and `to` (the
#'     allele's share of the bar, from the bottom) and the `depth` there
#'   * `modification_coverage`: the band's stack of modification states under
#'     `display$baseColor`: `position`, `from`, `to`, `depth`, `color`, `opacity`
#'   * `interbase_coverage`: the band's histogram of the reads' insertions,
#'     soft clips and hard clips at each position, stacked in that order:
#'     `position`, `type` (`"insertion"`, `"softclip"` or `"hardclip"`),
#'     `from` and `to` (in reads, from the band's top), `color`
#'   * `indicators`: the positions where those events pass the display's
#'     share of the local depth: `position`, `type` (the most common there),
#'     `color`
#'   * `junctions`: the reads' splice gaps merged, each once: `start`, `end`,
#'     `count` (the reads supporting it), `strand` (1, -1, or 0 where neither
#'     the reads' tags nor the splice motif says) and `motif` (`"GT-AG
#'     (canonical)"`, `"non-canonical"`, or `""` where the region's sequence
#'     does not cover both ends)
#'   * `legend`: `label` and `color`, the key the display shows for `color`,
#'     its title as the frame's `"title"` attribute
#'   * `ramps`: the color bars the display keys (a linear tag or MAPQ, base
#'     quality), as [guides_jbrowse()] draws them
#'   * `arcs`: under `display$readConnections`, each connection the band draws:
#'     `start` and `end` (its feet), `shape` (`"arc"`, or the read cloud's
#'     `"flat"`, `"flatSplit"` for a split read, and `"unplaced"` for one whose
#'     partner is outside the region), `category`, `color`, `opacity`, `y`
#'     (where the band plots it in bp: an arc's radius, or the read cloud's
#'     insert size, jittered as the display jitters it), `span` (the distance
#'     itself), `support` (the reads behind it) and `width` (the stroke in px)
#'   * `arc_ticks`: the interchromosomal ticks: `position`, `support`,
#'     `partners` (the contigs reached), `color`, `width`
#'   * `arc_legend`: the band's key, its title as the `"title"` attribute
#'   * `sections`: under `display$facet`, each section's `key` and `label` in stacking
#'     order; ungrouped, one section labelled `""`
#'
#'   Every table but the keys carries `section`, a factor in stacking order,
#'   and each section's rows start at 0. `read` is a row number of `reads`. `fade` is the share of the depth a
#'   mark's site holds, zeroed under the display's depth-dependent threshold,
#'   which the pileup fades a sub-pixel mark by.
#' @export
#' @examples
#' bam <- system.file("extdata", "volvox-sorted.bam", package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "ctgA", 1000, 3000,
#'   display = list(color = list(field = "strand")))
#' head(a$reads)
#' a$legend
jbrowse_alignments <- function(file, region, start = NULL, end = NULL,
                               reference = NULL, reference_chrom = NULL,
                               index = NULL, display = NULL,
                               byte_limit = 2e8) {
  adapter <- adapter_of(file, function() alignments_adapter(file, index))
  if (adapter$type == "CramAdapter" && is.null(reference)) {
    stop("a CRAM stores its reads against the reference it was encoded with: ",
      "pass that FASTA or 2bit as `reference =`", call. = FALSE)
  }
  read_regions(region, start, end, function(chrom, start, end, share, i)
    alignments_region(adapter, chrom, start, end, reference,
      reference_chrom[i] %||% chrom, display, byte_limit))
}

alignments_region <- function(adapter, chrom, start, end, reference,
                              reference_chrom, display, byte_limit) {
  display <- display %||% list()
  display$color <- config_color(display$color)
  if (!is.null(display$facet)) display$facet <- in_json(display$facet,
    c("domain", "hidden"))
  if (!is.null(display$modifications)) {
    display$modifications <- in_json(display$modifications,
      "shownModifications")
  }
  args <- list(
    adapterConfig = adapter,
    region = list(refName = chrom, originalRefName = reference_chrom,
      start = start, end = end),
    byteLimit = byte_limit)
  if (length(display)) args$display <- display
  if (!is.null(reference)) args$sequenceAdapter <- sequence_adapter(reference)
  t <- call_engine("alignments", args)
  sections <- as_frame(t$sections)
  frame <- function(table) {
    d <- as_frame(table)
    if (!is.null(d$section)) {
      d$section <- factor(d$section, levels = sections$label)
    }
    d
  }
  reads <- frame(t$reads)
  reads$read <- seq_len(nrow(reads))
  names(reads)[names(reads) == "insertSize"] <- "insert_size"
  names(reads)[names(reads) == "pairOrientation"] <- "pair_orientation"
  on_read <- function(table) {
    d <- frame(table)
    d$read <- d$read + 1
    d$row <- reads$row[d$read]
    d
  }
  colors <- jbrowse_colors()
  muted <- paints_modifications(display)
  if (muted) colors$bases[] <- MUTED_BASE
  reads$color <- jbrowse_color(reads$color)
  mismatches <- on_read(t$mismatches)
  mismatches$color <- jbrowse_color(colors$bases[mismatches$base])
  gaps <- on_read(t$gaps)
  gaps$color <- jbrowse_color(ifelse(gaps$type == "skip", colors$skip,
    colors$deletion))
  insertions <- on_read(t$insertions)
  insertions$color <- jbrowse_color(rep(colors$insertion, nrow(insertions)))
  clips <- on_read(t$clips)
  clips$color <- jbrowse_color(ifelse(clips$type == "soft", colors$softclip,
    colors$hardclip))
  modifications <- on_read(t$modifications)
  modifications$color <- jbrowse_color(modifications$color)
  qualities <- on_read(t$qualities)
  qualities$color <- jbrowse_color(qualities$color)
  coverage <- frame(t$coverage)
  coverage$color <- jbrowse_color(rep(colors$coverage, nrow(coverage)))
  alleles <- frame(t$alleleCoverage)
  alleles$color <- jbrowse_color(colors$bases[alleles$base])
  modification_coverage <- frame(t$modificationCoverage)
  modification_coverage$color <- jbrowse_color(modification_coverage$color)
  interbase_colors <- c(insertion = colors$insertion,
    softclip = colors$softclip, hardclip = colors$hardclip)
  interbase_coverage <- frame(t$interbaseCoverage)
  interbase_coverage$color <- jbrowse_color(
    unname(interbase_colors[interbase_coverage$type]))
  indicators <- frame(t$indicators)
  indicators$color <- jbrowse_color(unname(interbase_colors[indicators$type]))
  junctions <- frame(t$junctions)
  strand_colors <- c("1" = colors$forward, "-1" = colors$reverse,
    "0" = colors$read)
  junctions$color <- jbrowse_color(strand_colors[as.character(junctions$strand)])
  arcs <- frame(t$arcs)
  arcs$color <- jbrowse_color(arcs$color)
  arc_ticks <- frame(t$arcTicks)
  arc_ticks$color <- jbrowse_color(arc_ticks$color)
  legend <- frame(t$legend)
  attr(legend, "title") <- t$legendTitle
  arc_legend <- frame(t$arcLegend)
  attr(arc_legend, "title") <- t$arcLegendTitle
  base_key <- if (muted) key_frame(character(0), character(0), "Base")
    else key_frame(names(colors$bases), colors$bases, "Base")
  junction_key <- key_frame(c("Forward", "Reverse", "Unknown"), strand_colors,
    "Junction strand")
  junction_key <- junction_key[junction_key$color %in% junctions$color, ,
    drop = FALSE]
  attr(junction_key, "title") <- "Junction strand"
  structure(
    list(reads = reads, mismatches = mismatches, gaps = gaps,
      insertions = insertions, clips = clips, modifications = modifications,
      qualities = qualities,
      coverage = coverage, allele_coverage = alleles,
      modification_coverage = modification_coverage,
      interbase_coverage = interbase_coverage, indicators = indicators,
      junctions = junctions,
      legend = legend, arcs = arcs, arc_ticks = arc_ticks,
      arc_legend = arc_legend, arc_span_max = t$arcSpanMax,
      read_connections = display$readConnections %||% "off",
      sections = sections,
      keys = list(reads = legend, bases = base_key, arcs = arc_legend,
        junctions = junction_key), ramps = ramp_frames(t$ramps),
      region = list(chrom = chrom, start = start, end = end),
      reference = if (is.null(reference)) NULL else
        list(file = reference, chrom = reference_chrom),
      settings = t$settings, glyph = t$glyph,
      missing_contig = t$missingContig),
    primary = "reads",
    class = c("jbrowse_alignments", "jbrowse_result"))
}

#' Read a region of the reference
#'
#' The bases of a region, read through JBrowse's own FASTA or 2bit reader.
#' @param reference an indexed FASTA, plain or bgzipped, or a 2bit, path or
#'   URL
#' @param chrom,start,end the region, 0-based half-open
#' @return one string of `end - start` bases
#' @export
#' @examples
#' fa <- system.file("extdata", "volvox.fa", package = "ggjbrowse")
#' jbrowse_sequence(fa, "ctgA", 1000, 1020)
jbrowse_sequence <- function(reference, chrom, start, end) {
  args <- list(sequenceAdapter = sequence_adapter(reference),
    region = list(refName = chrom, start = start, end = end))
  in_engine(function(ctx) ctx$call("jbrowse.sequence",
    as.character(jsonlite::toJSON(args, auto_unbox = TRUE)), await = TRUE))
}

# The grey JBrowse paints every mismatched base in while the per-base layer
# paints modifications, its palette's mutedSnpBase
MUTED_BASE <- "#888"

paints_modifications <- function(display) {
  is.list(display$baseColor) &&
    isTRUE(display$baseColor$field %in% c("modifications", "bisulfite"))
}

# A config object whose vector slots stay JSON arrays however many values
# they hold
in_json <- function(object, slots) {
  for (slot in intersect(slots, names(object))) {
    object[[slot]] <- I(object[[slot]])
  }
  object
}

# A config's color as it writes it, a bare string its constant, and a list
# slot staying a list in JSON however many values it holds
config_color <- function(color) {
  if (!is.list(color)) return(color)
  for (slot in intersect(c("domain", "range", "labels"), names(color))) {
    color[[slot]] <- I(color[[slot]])
  }
  color
}

#' @export
print.jbrowse_alignments <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_alignments> %s:%.0f-%.0f, %d reads in %d sections, %d mismatches, %d gaps, %d insertions, %d clips, %d modifications, %d junctions, %d arcs\n",
    x$region$chrom, x$region$start, x$region$end, nrow(x$reads),
    nrow(x$sections), nrow(x$mismatches), nrow(x$gaps),
    nrow(x$insertions), nrow(x$clips), nrow(x$modifications),
    nrow(x$junctions), nrow(x$arcs)))
  invisible(x)
}
