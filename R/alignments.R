alignments_adapter <- function(file, index) {
  a <- guess_adapter(file, index, list(
    list(pattern = "\\.cram$", type = "CramAdapter", slot = "cramLocation"),
    list(pattern = "", type = "BamAdapter", slot = "bamLocation",
      index = "BAI")))
  if (a$type == "CramAdapter") {
    a$craiLocation <- location(index %||% beside(file, ".crai"))
  }
  a
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
#' @param region the regions to read, read as one view and drawn end to end
#'   on one axis as a JBrowse view lays them out: loci as JBrowse's location
#'   box writes them (`"ctgA:1,001-9,000"`, 1-based and inclusive), contig
#'   names alone for whole contigs, a data frame of `chrom`, `start` and `end`
#'   (0-based, half-open, as BED writes them), or a GRanges. With `NULL`, the
#'   default, the reader reads the whole file: each contig it lists but those
#'   its index shows empty, to the length the file or `assembly` states. A whole
#'   file or contig is refused past 5 million rows in any one table or past
#'   `span_limit`
#' @param assembly the reference: a [jbrowse_assembly()], an assembly as a
#'   JBrowse config's `assemblies` entry writes it, or an indexed FASTA,
#'   plain or bgzipped, or a 2bit alone, path or URL. A CRAM and the
#'   bisulfite layer need its sequence, and a BAM read without an MD tag is
#'   compared against it. `region` may name a contig by any alias it lists,
#'   the engine reads each file under the file's own, and the reader names
#'   the contig in its tables as `region` does
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
#'     (1540 by default), `readName`, `tags` (a list of
#'     `list(tag = , value = )`), `spliced`, `svSignal` (SV evidence: a split
#'     read, an insertion or deletion of 50 bp or more, or a discordant pair),
#'     `alone` and `split` (each `"only"` or `"exclude"`); and which
#'     junctions and connections draw:
#'     `junctions = list(minReads = , nonCanonical = "exclude")`,
#'     `connections = list(interchrom = , offscreenMates = , minInterchromReads = )`
#'   * `baseColor`: the per-base layer painted over the reads:
#'     `list(field = "modifications")` for the reads' MM/ML base
#'     modifications, `list(field = "bisulfite")` for the methylation a
#'     bisulfite or EM-seq read shows against `assembly` (a retained C is
#'     methylated, a C read as T is not), `list(field = "baseQuality")` for
#'     each base's Phred score; unset paints none. Under `"modifications"` or
#'     `"bisulfite"` the display greys every mismatched base so the calls
#'     stand out
#'   * `modifications`: what that layer reads: any of `threshold` (the
#'     percent probability under which a call is hidden, 10), `unmodified`
#'     (`"hidden"`, `"calls"` or `"all"`; bisulfite paints its unmethylated
#'     calls under either of the last two), `cytosineContext` (`"CG"`,
#'     `"CHG"`, `"CHH"` or `"all"`) and `shownModifications` (the type codes
#'     drawn, every one when empty)
#'   * `sort`: how the rows are laid out: `list(type = "position")` (or
#'     `"length"`, `"spliced"`, `"split"`, `"svSignal"`: SV evidence first,
#'     split reads above SV-size indels above discordant pairs), or a column sort, `type`
#'     (`"strand"`, `"basePair"`, `"tag"`, `"insertion"`, `"softclip"`,
#'     `"hardclip"`) with `pos`, the 0-based column ranked over, and `tag`
#'   * `facet`: what splits the reads into stacked sections, each laid out on
#'     its own: `list(field = "strand")` (or `"pairOrientation"`,
#'     `"tags.HP"`, ...), with optional `domain` (section keys in stacking
#'     order) and `hidden` (keys left out)
#'   * `unit`: what one row stands for, `"read"`, or `"chain"` for a read
#'     with its mate and split segments on one row, a hairline across it, and
#'     the strand colors framed against the chain's own
#'   * the four bands, each a list with its own `show`:
#'     `coverage` (`height`, `showInterbaseIndicators`, `minAlleleFraction`,
#'     the share of the depth an allele needs to be colored);
#'     `pileup` (`readHeight`, 7, or 3 and 1 for compact and super-compact,
#'     `showMismatches`, `showSoftClipping`, `fadeLowFrequency`, ...);
#'     `sashimi` (`height`, `placement`, `showCounts`);
#'     `readConnections` (`show = TRUE` for each mate pair's and split
#'     read's arc, `mark = "cloud"` for the read cloud, `height`,
#'     `showOrdinaryPairs`, `color = list(field = )`),
#'     e.g. `list(coverage = list(show = FALSE), pileup = list(readHeight = 3))`.
#'     A band or layer set not to show (`coverage`, `pileup` or `sashimi`'s
#'     `show`, `pileup$showMismatches`, `coverage$showInterbaseIndicators`)
#'     reads as empty tables
#' @param byte_limit the most compressed bytes the regions may need, refused
#'   off the index before any read is fetched; the whole table is held in
#'   memory
#' @param span_limit the most bp the regions may span in all, refused before
#'   anything is fetched, since the coverage band holds every base of them
#' @param width the pixels the figure spans across the regions, which the
#'   sub-pixel fade and the read-connection band's reach read
#' @return a `jbrowse_alignments` list of data frames, beside the slots
#'   [jbrowse_result] lists:
#'   * `reads`: `read` (the read's key, from 1, which the marks below name
#'     it by), `id` (the record's own id, the same in every call, where its
#'     mates and supplementary segments share its `name`), `start`, `end`,
#'     `row`, `strand`, `mapq`, `flags`, `name`, `insert_size`,
#'     `pair_orientation`, `color` (the fill the display paints, as
#'     `#rrggbb`), `category` (the bucket the display files the read under:
#'     `fwdStrand`, `mapq`, `unmappedMate`, ...) and `section`. Over several regions a read is listed once in each region it
#'     reaches, its copies sharing `read` and `row`
#'   * `mismatches`: `position`, `base`, `read`, `frequency` (how many reads
#'     carry that base there, over the depth), `quality`, `fade`
#'   * `gaps`: `start`, `end`, `type` (`"deletion"` or `"skip"`), `read`,
#'     `fade`
#'   * `chains`: under `unit = "chain"`, each chain of more than one read:
#'     `start`, `end`, `row`, and the `color` and `opacity` of the hairline
#'     across it; empty otherwise
#'   * `overlaps`: under `unit = "chain"`, each span where two reads of one row
#'     both align: `start`, `end`, `row` and the neutral `color` the display
#'     fills it with, which is neither read's color; empty otherwise
#'   * `insertions`: `position`, `length`, `sequence`, `read`, `fade`
#'   * `clips`: `position`, `length`, `type` (`"soft"` or `"hard"`), `read`,
#'     `fade`
#'   * `modifications`: each base modification call the layer paints:
#'     `position`, `read`, `type` (the code, `m` for 5mC), `name`,
#'     `probability` (1 for a bisulfite call), `unmodified` (whether the
#'     call is that the base is not modified), and the `color` and `opacity`
#'     the display paints it at; empty without `display$baseColor`
#'   * `qualities`: under `baseColor` `"baseQuality"`, each aligned base:
#'     `position`, `read`, `score` (Phred) and the `color` the display
#'     paints it
#'   * `coverage`: `start`, `end`, `depth`, `forward`, `reverse`, over the
#'     region alone, where the reads reach past it
#'   * `allele_coverage`: the coverage band's stack at each mismatch position,
#'     as the display stacks it: `position`, `base`, `from` and `to` (the
#'     allele's share of the bar, from the bottom) and the `depth` there
#'   * `modification_coverage`: the band's stack of modification states under
#'     `display$baseColor`: `position`, `from`, `to`, `depth`, `color`, `opacity`
#'   * `interbase_coverage`: the band's histogram of the reads' insertions,
#'     soft clips and hard clips at each position, stacked in that order:
#'     `position`, `type` (`"insertion"`, `"softclip"` or `"hardclip"`),
#'     `from` and `to` (in reads, from the band's top), `color`
#'   * `indicators`: the positions the display marks with a triangle, where
#'     those events number over 30% of a depth of 8 reads or more:
#'     `position`, `type` (the most common there), `color`
#'   * `junctions`: the reads' splice gaps merged, each once over every
#'     region: `start`, `end`,
#'     `count` (the reads supporting it), `strand` (1, -1, or 0 where neither
#'     the reads' tags nor the splice motif says) and `motif` (`"GT-AG
#'     (canonical)"`, `"non-canonical"`, or `""` where the region's sequence
#'     does not cover both ends)
#'   * `arcs`: under `display$readConnections`, each connection the band draws:
#'     `start` and `end` (its feet), `shape` (`"arc"`, or the read cloud's
#'     `"flat"`, `"flatSplit"` for a split read, and `"unplaced"` for one whose
#'     partner is outside the region), `category`, `color`, `opacity`, `y`
#'     (where the band plots it in bp: an arc's radius, or the read cloud's
#'     insert size, jittered as the display jitters it), `span` (the distance
#'     itself), `support` (the reads behind it) and `width` (the stroke in px)
#'   * `arc_ticks`: the interchromosomal ticks: `position`, `support`,
#'     `partners` (the contigs reached), `color`, `width`
#'   * `sections`: under `display$facet`, each section's `key` and `label` in stacking
#'     order; ungrouped, one section labelled `""`
#'
#'   `keys` holds `reads`, `bases` (empty while a modifications or bisulfite
#'   layer greys the bases), `junctions`, `interbase` and `arcs`, and
#'   `y_scale` the coverage band's depth axis, shared by every section.
#'   A read's marks (`mismatches`, `gaps`, `insertions`, `clips`,
#'   `modifications`, `qualities`) name it by `read` and their `region`, and
#'   take its `row` and `section` from `reads` when drawn, so reassigning
#'   `reads$row` moves them; assigning `reads` drops the marks of the reads it
#'   leaves out. The other tables have `section`, a factor in stacking
#'   order, each section's rows starting at 0. Every table has the 1-based
#'   `region` each row was read in, with its `chrom`, and every mark table
#'   `color`, the color the display paints the mark. The coverage band's
#'   tables, the junctions and the arcs count every read the engine read and
#'   stay as they were when `reads` is filtered: filter with `display =
#'   list(filter = ...)` to recompute them. `fade` is the share of the depth a mark's site holds, zeroed
#'   under the display's depth-dependent threshold, which the pileup fades a
#'   sub-pixel mark by.
#'   `read_pitch` gives the `featureHeight` a read draws at and the
#'   `featureSpacing` under it.
#' @export
#' @examples
#' bam <- system.file("extdata", "gapdh_direct_rna.bam",
#'   package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "chr12:6,534,201-6,538,700",
#'   display = list(color = list(field = "strand")))
#' head(a$reads)
#' a$keys$reads
jbrowse_alignments <- function(file, region = NULL, display = NULL,
                               width = DEFAULT_WIDTH, assembly = NULL,
                               index = NULL,
                               byte_limit = 2e8, span_limit = 1e7) {
  adapter <- adapter_of(file, function() alignments_adapter(file, index))
  assembly <- assembly_config(assembly)
  if (adapter$type == "CramAdapter" && !has_sequence(assembly)) {
    stop("a CRAM stores its reads against the reference it was encoded with: ",
      "pass that FASTA or 2bit as `assembly =`", call. = FALSE)
  }
  if (paints_bisulfite(display) && !has_sequence(assembly)) {
    stop("bisulfite calls compare each read's cytosines to the reference: ",
      "pass that FASTA or 2bit as `assembly =`", call. = FALSE)
  }
  regions <- reader_regions(adapter, region, assembly)
  span <- sum(regions$end - regions$start)
  if (is.null(region) && span > span_limit) {
    rlang::abort(sprintf(paste("the whole file spans %s bp over %d contigs,",
      "past `span_limit` (%s bp): pass a region, or raise `span_limit`"),
      format(span, big.mark = ",", scientific = FALSE), nrow(regions),
      format(span_limit, big.mark = ",", scientific = FALSE)),
      class = "ggjbrowse_too_large", limit = "span_limit", call = NULL)
  }
  args <- list(adapterConfig = adapter, regions = engine_regions(regions),
    byteLimit = byte_limit, spanLimit = span_limit)
  args$assembly <- assembly
  if (isTRUE(attr(regions, "whole"))) args$rowLimit <- WHOLE_FILE_ROWS
  args$display <- display_config(display)
  t <- call_engine("alignments", args)
  x <- alignments_result(t, assembly)
  x$junctions$region <- junction_region(x$junctions,
    unlist(t$readAs$refName), regions)
  x$junctions$refName <- NULL
  x <- in_one_call(x, regions, t, width)
  x$sequence <- if (has_sequence(assembly)) region_bases(assembly, regions)
  x
}

# Each region's reference bases, NA for a region the assembly cannot read,
# so the sequence band draws from the result alone
region_bases <- function(assembly, regions) {
  vapply(seq_len(nrow(regions)), function(i) tryCatch(
    jbrowse_sequence(assembly, regions[i, c("chrom", "start", "end")]),
    error = function(e) NA_character_), "")
}

alignments_result <- function(t, assembly) {
  sections <- as_frame(t$sections)
  frame <- function(table) {
    d <- as_frame(table)
    if (!is.null(d$section)) {
      d$section <- factor(d$section, levels = sections$label)
    }
    d
  }
  reads <- frame(t$reads)
  reads <- cbind(read = reads$key + 1L, reads[names(reads) != "key"])
  names(reads)[names(reads) == "insertSize"] <- "insert_size"
  names(reads)[names(reads) == "pairOrientation"] <- "pair_orientation"
  on_read <- function(table) {
    d <- as_frame(table)
    d$read <- reads$read[d$read + 1L]
    d$section <- NULL
    d
  }
  mismatches <- on_read(t$mismatches)
  chains <- frame(t$chains)
  overlaps <- frame(t$overlaps)
  gaps <- on_read(t$gaps)
  insertions <- on_read(t$insertions)
  clips <- on_read(t$clips)
  modifications <- on_read(t$modifications)
  qualities <- on_read(t$qualities)
  coverage <- frame(t$coverage)
  alleles <- frame(t$alleleCoverage)
  modification_coverage <- frame(t$modificationCoverage)
  interbase_coverage <- frame(t$interbaseCoverage)
  indicators <- frame(t$indicators)
  junctions <- frame(t$junctions)
  arcs <- frame(t$arcs)
  arc_ticks <- frame(t$arcTicks)
  read_key <- frame(t$legend)
  if (nrow(overlaps)) {
    read_key <- rbind(read_key, data.frame(
      label = "Pair/chain reads overlap here", color = overlaps$color[1]))
  }
  read_key <- structure(read_key, title = t$legendTitle)
  arc_key <- structure(frame(t$arcLegend), title = t$arcLegendTitle)
  engine_key <- function(key, title) {
    key_frame(unlist(key$label) %||% character(0),
      unlist(key$color) %||% character(0), title)
  }
  base_key <- engine_key(t$baseLegend, t$baseLegendTitle)
  interbase_key <- engine_key(t$interbaseLegend, t$interbaseLegendTitle)
  junction_key <- engine_key(t$junctionLegend, t$junctionLegendTitle)
  structure(
    list(reads = reads, mismatches = mismatches, gaps = gaps, chains = chains,
      overlaps = overlaps,
      insertions = insertions, clips = clips, modifications = modifications,
      qualities = qualities,
      coverage = coverage, allele_coverage = alleles,
      modification_coverage = modification_coverage,
      interbase_coverage = interbase_coverage, indicators = indicators,
      junctions = junctions, arcs = arcs, arc_ticks = arc_ticks,
      sections = sections, display = t$config,
      keys = list(reads = read_key, bases = base_key, arcs = arc_key,
        junctions = junction_key, interbase = interbase_key),
      ramps = ramp_frames(t$ramps), y_scale = y_scale(t$coverageScale),
      fields = list(),
      read_pitch = t$readPitch, glyph = t$glyph, arc_span_max = t$arcSpanMax),
    primary = "reads",
    class = c("jbrowse_alignments", "jbrowse_result"))
}

#' Read a region of the reference
#'
#' The bases of each region, read through JBrowse's own FASTA or 2bit
#' reader.
#' @param assembly an indexed FASTA, plain or bgzipped, or a 2bit, path or
#'   URL; or a [jbrowse_assembly()], whose aliases a region may then name a
#'   contig by
#' @param region the regions, as a reader's `region` takes them: loci such as
#'   `"chr12:6,534,501-6,534,520"`, contig names for whole contigs, a data
#'   frame of `chrom`, `start` and `end` (0-based, half-open), or a GRanges
#' @return a string of each region's bases, `end - start` of them
#' @export
#' @examples
#' fa <- system.file("extdata", "gapdh_hg38.fa.gz", package = "ggjbrowse")
#' jbrowse_sequence(fa, "chr12:6,534,501-6,534,520")
jbrowse_sequence <- function(assembly, region) {
  assembly <- assembly_config(assembly)
  regions <- regions_of(region, function()
    assembly_contigs(assembly))
  vapply(seq_len(nrow(regions)), function(i) {
    args <- list(region = list(refName = regions$chrom[i],
      start = regions$start[i], end = regions$end[i]), assembly = assembly)
    in_engine(function(ctx) ctx$call("jbrowse.sequence",
      as_engine_json(args), await = TRUE))
  }, "")
}

# The 0-based region a junction draws in: the one holding its start, else its
# end, else the first on its contig that it crosses
junction_region <- function(j, contigs, regions) {
  region <- ifelse(j$startRegion >= 0, j$startRegion, j$endRegion)
  for (k in which(region < 0)) {
    region[k] <- c(which(contigs == j$refName[k] & j$start[k] < regions$end &
      j$end[k] > regions$start), 1L)[1] - 1L
  }
  as.integer(region)
}

# A band's y axis as the display autoscales it: its domain, its type
# (`"linear"`, `"log"` or `"symlog"`), the symlog's constant, and its ticks
y_scale <- function(s) {
  list(domain = unlist(s$domain), type = s$type %||% s$scaleType,
    symlog_constant = s$symlogConstant, ticks = unlist(s$ticks))
}

paints_bisulfite <- function(display) {
  is.list(display$baseColor) && identical(display$baseColor$field, "bisulfite")
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
config_color <- function(color) in_json(color, c("domain", "range", "labels"))

# A display config as the engine reads it, each slot that lists values a JSON
# array however many it holds; NULL for an empty one
display_config <- function(display) {
  if (!length(display)) return(NULL)
  display$color <- config_color(display$color)
  if (is.character(display$filter)) display$filter <- I(display$filter)
  display$facet <- in_json(display$facet, c("domain", "hidden"))
  display$modifications <- in_json(display$modifications,
    "shownModifications")
  display
}

#' @export
print.jbrowse_alignments <- function(x, ...) {
  cat(sprintf(
    "<jbrowse_alignments> %s, %d reads in %d sections, %d mismatches, %d gaps, %d insertions, %d clips, %d indicators, %d modifications, %d junctions, %d arcs\n",
    print_loci(x), nrow(x$reads),
    nrow(x$sections), nrow(x$mismatches), nrow(x$gaps),
    nrow(x$insertions), nrow(x$clips), nrow(x$indicators),
    nrow(x$modifications), nrow(x$junctions), nrow(x$arcs)))
  invisible(x)
}
