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
#' @param region the region: a locus as JBrowse's location box writes it
#'   (`"ctgA:1,001-9,000"`, 1-based and inclusive), several of them, which
#'   the engine reads as one view and draws end to end on one axis as a
#'   JBrowse view lays them out, a data frame of `chrom`, `start` and `end`, or a contig name
#'   with `start` and `end`. With `NULL`, the default, the reader reads the
#'   whole file, each contig it lists, and refuses with a message past 5
#'   million rows in any one table
#' @param start,end with a contig name as `region`, the region's bounds,
#'   0-based half-open
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
#'     (1540 by default), `readName`, `tagFilters` (a list of
#'     `list(tag = , value = )`), `spliced`, `svSignal` (SV evidence: a split
#'     read, an insertion or deletion of 50 bp or more, or a discordant pair),
#'     `singletons` and `split` (each `"only"` or `"exclude"`)
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
#'   * `readConnections`: the read-connection band, `"arc"` for each mate
#'     pair's and split read's arc or `"cloud"` for the read cloud; unset
#'     reads none
#'   * `showCoverage`, `showPileup`, `showSashimiArcs`, and any other slot
#'     such as `showOrdinaryPairs`, `arcColor` or `minSashimiScore`. A band
#'     set not to show (`showCoverage`, `showPileup`, `showMismatches`,
#'     `showInterbaseIndicators`, `showSashimiArcs`) reads as empty tables
#' @param byte_limit the most compressed bytes the regions may need, refused
#'   off the index before any read is fetched; the whole table is held in
#'   memory
#' @param span_limit the most bp the regions may span in all, refused before
#'   anything is fetched, since the coverage band holds every base of them
#' @return a `jbrowse_alignments` list of data frames:
#'   * `reads`: `read` (the row number the other tables' `read` column names
#'     it by, kept through filtering, within this call), `id` (the record's
#'     own id, the same in every call, where its mates and supplementary
#'     segments share its `name`), `start`, `end`, `row`, `strand`, `mapq`,
#'     `flags`, `name`,
#'     `insert_size`, `pair_orientation`, `color` (the fill the display paints,
#'     as `#rrggbb`) and `category` (the bucket the display files the read
#'     under: `fwdStrand`, `mapq`, `unmappedMate`, ...). Over several
#'     regions a read is listed once in each region it reaches, its copies
#'     sharing `key` and `row`
#'   * `mismatches`: `position`, `base`, `read`, `row`, `frequency` (how many
#'     reads carry that base there, over the depth), `quality`, `fade`
#'   * `gaps`: `start`, `end`, `type` (`"deletion"` or `"skip"`), `read`, `row`,
#'     `fade`
#'   * `chains`: under `unit = "chain"`, each chain of more than one read:
#'     `start`, `end`, `row`, and the `color` and `opacity` of the hairline
#'     across it; empty otherwise
#'   * `overlaps`: under `unit = "chain"`, each span where two reads of one row
#'     both align: `start`, `end`, `row` and the neutral `color` the display
#'     fills it with, which is neither read's color; empty otherwise
#'   * `insertions`: `position`, `length`, `sequence`, `read`, `row`, `fade`
#'   * `clips`: `position`, `length`, `type` (`"soft"` or `"hard"`), `read`,
#'     `row`, `fade`
#'   * `modifications`: each base modification call the layer paints:
#'     `position`, `read`, `row`, `type` (the code, `m` for 5mC), `name`,
#'     `probability` (1 for a bisulfite call), `unmodified` (whether the
#'     call is that the base is not modified), and the `color` and `opacity`
#'     the display paints it at; empty without `display$baseColor`
#'   * `qualities`: under `baseColor` `"baseQuality"`, each aligned base:
#'     `position`, `read`, `row`, `score` (Phred) and the `color` the display
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
#'   * `coverage_scale`: the band's depth axis as the display autoscales it,
#'     shared by every section: `domain`, `type` and `ticks`
#'   * `junctions`: the reads' splice gaps merged, each once over every
#'     region: `start`, `end`,
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
#'   * `keys`: the display's keys, each titled as the browser titles it:
#'     `reads`, `bases` (empty while a modifications or bisulfite layer greys
#'     the bases), `junctions`, `interbase` and `arcs`
#'
#'   Every table but the keys carries `section`, a factor in stacking order,
#'   and each section's rows start at 0, and every mark table carries `color`,
#'   the color the display paints the mark. `read` is a row number of
#'   `reads`. `fade` is the share of the depth a mark's site holds, zeroed
#'   under the display's depth-dependent threshold, which the pileup fades a
#'   sub-pixel mark by.
#'   `display` is the display's config as JBrowse resolves it, every slot
#'   under JBrowse's own name and the schema's defaults filled in. `read_pitch` gives the
#'   `featureHeight` a read draws at and the `featureSpacing` under it.
#' @export
#' @examples
#' bam <- system.file("extdata", "gapdh_direct_rna.bam",
#'   package = "ggjbrowse")
#' a <- jbrowse_alignments(bam, "chr12", 6534200, 6538700,
#'   display = list(color = list(field = "strand")))
#' head(a$reads)
#' a$legend
jbrowse_alignments <- function(file, region = NULL, start = NULL, end = NULL,
                               display = NULL, assembly = NULL, index = NULL,
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
  regions <- reader_regions(adapter, region, start, end)
  args <- list(adapterConfig = adapter, regions = engine_regions(regions),
    byteLimit = byte_limit, spanLimit = span_limit)
  args$assembly <- assembly
  if (is.null(region)) args$rowLimit <- WHOLE_FILE_ROWS
  args$display <- display_config(display)
  t <- call_engine("alignments", args)
  x <- alignments_result(t, assembly)
  x$junctions$region <- junction_region(x$junctions,
    unlist(t$readAs$refName), regions)
  x$junctions$refName <- NULL
  in_one_call(x, regions, t$missingContigs)
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
  reads$read <- seq_len(nrow(reads))
  names(reads)[names(reads) == "insertSize"] <- "insert_size"
  names(reads)[names(reads) == "pairOrientation"] <- "pair_orientation"
  on_read <- function(table) {
    d <- frame(table)
    d$read <- d$read + 1
    d$row <- reads$row[d$read]
    d
  }
  painted <- function(d) {
    d$color <- jbrowse_color(d$color)
    d
  }
  reads <- painted(reads)
  mismatches <- painted(on_read(t$mismatches))
  chains <- painted(frame(t$chains))
  overlaps <- painted(frame(t$overlaps))
  gaps <- painted(on_read(t$gaps))
  insertions <- painted(on_read(t$insertions))
  clips <- painted(on_read(t$clips))
  modifications <- painted(on_read(t$modifications))
  qualities <- painted(on_read(t$qualities))
  coverage <- painted(frame(t$coverage))
  alleles <- painted(frame(t$alleleCoverage))
  modification_coverage <- painted(frame(t$modificationCoverage))
  interbase_coverage <- painted(frame(t$interbaseCoverage))
  indicators <- painted(frame(t$indicators))
  junctions <- painted(frame(t$junctions))
  arcs <- painted(frame(t$arcs))
  arc_ticks <- painted(frame(t$arcTicks))
  legend <- frame(t$legend)
  if (nrow(overlaps)) {
    legend <- rbind(legend, data.frame(label = "Pair/chain reads overlap here",
      color = overlaps$color[1]))
  }
  legend <- structure(legend, title = t$legendTitle, complete = TRUE)
  arc_legend <- structure(frame(t$arcLegend), title = t$arcLegendTitle)
  engine_key <- function(key, title, complete = FALSE) {
    key_frame(unlist(key$label) %||% character(0),
      unlist(key$color) %||% character(0), title, complete)
  }
  base_key <- engine_key(t$baseLegend, t$baseLegendTitle, complete = TRUE)
  interbase_key <- engine_key(t$interbaseLegend, t$interbaseLegendTitle,
    complete = TRUE)
  junction_key <- engine_key(t$junctionLegend, t$junctionLegendTitle)
  structure(
    list(reads = reads, mismatches = mismatches, gaps = gaps, chains = chains,
      overlaps = overlaps,
      insertions = insertions, clips = clips, modifications = modifications,
      qualities = qualities,
      coverage = coverage, allele_coverage = alleles,
      modification_coverage = modification_coverage,
      interbase_coverage = interbase_coverage, indicators = indicators,
      junctions = junctions,
      legend = legend, arcs = arcs, arc_ticks = arc_ticks,
      arc_legend = arc_legend, arc_span_max = t$arcSpanMax,
      read_connections = t$config$readConnections,
      sections = sections,
      keys = list(reads = legend, bases = base_key, arcs = arc_legend,
        junctions = junction_key, interbase = interbase_key),
      ramps = ramp_frames(t$ramps), coverage_scale = t$coverageScale,
      reference = if (has_sequence(assembly))
        list(assembly = assembly, chrom = t$readAs$originalRefName[[1]]),
      settings = t$settings, display = t$config, read_pitch = t$readPitch,
      glyph = t$glyph),
    primary = "reads",
    class = c("jbrowse_alignments", "jbrowse_result"))
}

#' Read a region of the reference
#'
#' The bases of a region, read through JBrowse's own FASTA or 2bit reader.
#' @param assembly an indexed FASTA, plain or bgzipped, or a 2bit, path or
#'   URL; or a [jbrowse_assembly()], whose aliases then name `chrom`
#' @param chrom,start,end the region, 0-based half-open
#' @return one string of `end - start` bases
#' @export
#' @examples
#' fa <- system.file("extdata", "gapdh_hg38.fa.gz", package = "ggjbrowse")
#' jbrowse_sequence(fa, "chr12", 6534500, 6534520)
jbrowse_sequence <- function(assembly, chrom, start, end) {
  args <- list(region = list(refName = chrom, start = start, end = end),
    assembly = assembly_config(assembly))
  in_engine(function(ctx) ctx$call("jbrowse.sequence",
    as_engine_json(args), await = TRUE))
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
