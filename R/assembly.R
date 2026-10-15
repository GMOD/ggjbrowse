#' An assembly: the reference sequence and its contigs' aliases
#'
#' The assembly a JBrowse config's `assemblies` entry writes. Each reader
#' takes it as `assembly`, as each track of a JBrowse config names its
#' assembly:
#' * its sequence is the reference a CRAM decodes against, the bisulfite
#'   layer calls cytosines against, a BAM without MD tags finds its
#'   mismatches against, and a CDS translates against
#' * for files that spell one contig several ways (`chr7` in a BAM,
#'   `NC_000007.13` in a GFF, `7` in the FASTA), a region may name a contig by
#'   any of its aliases, and the engine reads each file under the alias that
#'   file uses, as JBrowse renames a track's regions
#' * the reader names the contig in its tables as the region does, so
#'   [jbrowse_stack()] draws files that spell it differently on one axis
#'
#' A reader also takes the sequence's path alone as `assembly`, an assembly
#' with no aliases, or a JBrowse config's own assembly as a list:
#' `list(name = "hg19", uri = "hg19.fa.gz", refNameAliases = list(uri =
#' "hg19_aliases.txt"))`.
#' @param sequence the reference: an indexed FASTA, plain or bgzipped, or a
#'   2bit, path or URL
#' @param aliases the alias table, path or URL: tab-separated, a row per
#'   contig listing its names, the first column the sequence's, as UCSC's
#'   `chromAlias.txt` writes it; or NCBI's `sequence_report.tsv`
#' @param name the assembly's name
#' @return a `jbrowse_assembly`: the list a config's `assemblies` entry writes
#' @export
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' aliases <- tempfile(fileext = ".txt")
#' writeLines("chr12\t12\tNC_000012.12", aliases)
#' hg38 <- jbrowse_assembly(e("gapdh_hg38.fa.gz"), aliases)
#' a <- jbrowse_alignments(e("gapdh_direct_rna.bam"), "12:6,534,201-6,538,700",
#'   assembly = hg38)
#' a$regions$chrom
jbrowse_assembly <- function(sequence = NULL, aliases = NULL,
                             name = "assembly") {
  config <- list(name = name, uri = sequence)
  if (!is.null(aliases)) config$refNameAliases <- list(uri = aliases)
  structure(config, class = "jbrowse_assembly")
}

#' @export
print.jbrowse_assembly <- function(x, ...) {
  cat(sprintf("<jbrowse_assembly> %s: sequence %s, aliases %s\n", x$name,
    x$uri %||% "none", x$refNameAliases$uri %||% "none"))
  invisible(x)
}

sequence_adapter <- function(reference) {
  a <- guess_adapter(reference, NULL, list(
    list(pattern = "\\.2bit$", type = "TwoBitAdapter",
      slot = "twoBitLocation"),
    list(pattern = "\\.gz$", type = "BgzipFastaAdapter",
      slot = "fastaLocation"),
    list(pattern = "", type = "IndexedFastaAdapter", slot = "fastaLocation")))
  if (a$type == "TwoBitAdapter") return(a)
  a$faiLocation <- location(beside(reference, ".fai"))
  if (a$type == "BgzipFastaAdapter") {
    a$gziLocation <- location(beside(reference, ".gzi"))
  }
  a
}

#' The contigs of an assembly or a file
#'
#' Each contig from 0 to its length, one shape whatever `x` is: a region
#' frame any reader takes, to read a whole contig or genome.
#' * an assembly lists the contigs its sequence names
#' * a linear file (BAM, CRAM, BigWig, VCF, GFF3, GTF, BED, BigBed) lists the
#'   contigs it names, each to the length the file, or `assembly`, states,
#'   and whether its index shows records there
#' * a pairwise alignment file (PAF, PIF, delta, chain, MashMap) lists each
#'   side's contigs in the order the file first names them, with the length
#'   it states: the regions a whole-genome [jbrowse_dotplot()] or
#'   [jbrowse_synteny()] lays end to end where it is given none
#' @param x a [jbrowse_assembly()], an assembly as a JBrowse config's
#'   `assemblies` entry writes it, or a FASTA or 2bit path or URL; or a
#'   linear or pairwise file, a path or URL, or a JBrowse track's `adapter`
#'   config as a list
#' @param assembly for a linear file, the assembly whose contigs' lengths
#'   stand where the file states none
#' @param index the file's index, if not beside it under its usual name
#' @return a data frame of `assembly` (the assembly's name, the side's for a
#'   pairwise file, which a file read from its path names `"query"` or
#'   `"target"`, and `NA` for a linear file read with no `assembly`),
#'   `chrom`, `start` and `end`, 0-based half-open, `end` `NA` where a
#'   pairwise file does not state the length; a linear file adds `has_data`,
#'   `FALSE` where its index shows no records and `NA` where it cannot tell
#' @export
#' @examples
#' e <- function(f) system.file("extdata", f, package = "ggjbrowse")
#' jbrowse_contigs(e("gapdh_hg38.fa.gz"))
#' jbrowse_contigs(e("volvox-sorted.bam"))
#' head(jbrowse_contigs(e("peach_grape.paf.gz")))
jbrowse_contigs <- function(x, assembly = NULL, index = NULL) {
  kind <- contigs_kind(x)
  if (kind == "assembly") {
    config <- assembly_config(x)
    contigs <- region_frame(call_engine("regions", list(assembly = config)))
    return(cbind(assembly = rep(config$name %||% NA_character_,
      nrow(contigs)), contigs))
  }
  if (kind == "pairwise") {
    sides <- contigs_of(comparative_args(x, index, NULL)$adapterConfig)
    return(data.frame(assembly = sides$assembly, chrom = sides$chrom,
      start = 0, end = sides$length))
  }
  config <- assembly_config(assembly)
  contigs <- file_contigs(adapter_of(x, function() linear_adapter(x, index)),
    config)
  cbind(assembly = rep(config$name %||% NA_character_, nrow(contigs)),
    contigs)
}

PAIRWISE_ADAPTERS <- c("PAFAdapter", "PairwiseIndexedPAFAdapter",
  "DeltaAdapter", "ChainAdapter", "MashMapAdapter", "MCScanAnchorsAdapter",
  "MCScanSimpleAnchorsAdapter", "BlastTabularAdapter", "AllVsAllPAFAdapter")

# Whether `x` is an assembly, a pairwise alignment file or a linear file:
# a list by its adapter `type`, which an assembly config lacks, and a path by
# its extension
contigs_kind <- function(x) {
  if (inherits(x, "jbrowse_assembly")) return("assembly")
  if (is.list(x)) {
    if (is.null(x$type)) return("assembly")
    return(if (x$type %in% PAIRWISE_ADAPTERS) "pairwise" else "linear")
  }
  name <- tolower(bare_name(x))
  if (grepl("\\.(fa|fasta|fna|fas)(\\.gz)?$|\\.2bit$", name)) "assembly"
  else if (grepl("\\.(paf(\\.gz)?|pif\\.gz|delta|chain|out)$", name)) "pairwise"
  else "linear"
}

# A linear file's adapter, by the reader whose extension it has
linear_adapter <- function(file, index) {
  name <- tolower(bare_name(file))
  if (grepl("\\.(bam|cram)$", name)) alignments_adapter(file, index)
  else if (grepl("\\.(bw|bigwig)$", name)) {
    list(type = "BigWigAdapter", bigWigLocation = location(file))
  } else features_adapter(file, index)
}

# An assembly config as the engine reads it, built afresh at each read: a
# sequence's path or URL alone as jbrowse_assembly() writes it, the sequence
# and alias files a config names by `uri` as the adapters sequence_adapter()
# and alias_adapter() make of them, each location carrying its file's
# modification time so a rewritten file reads anew, and any other `uri` a
# path or URL made readable, as a shorthand the engine expands as JBrowse's
# config does
assembly_config <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  if (is.character(assembly) && length(assembly) == 1) {
    assembly <- jbrowse_assembly(assembly)
  }
  if (!is.list(assembly)) {
    stop("an assembly is a jbrowse_assembly(), a JBrowse assembly config as ",
      "a list, or a FASTA or 2bit path", call. = FALSE)
  }
  assembly <- unclass(assembly)
  if (is.character(assembly$uri) && is.null(assembly$sequence)) {
    assembly$sequence <- list(type = "ReferenceSequenceTrack",
      trackId = paste0(assembly$name, "-ReferenceSequenceTrack"),
      adapter = sequence_adapter(assembly$uri))
    assembly$uri <- NULL
  }
  aliases <- assembly$refNameAliases
  if (is.character(aliases)) aliases <- list(uri = aliases)
  if (is.character(aliases$uri) && is.null(aliases$adapter)) {
    assembly$refNameAliases <- list(adapter = alias_adapter(aliases$uri))
  }
  readable <- function(x) {
    if (!is.list(x)) return(x)
    x <- lapply(unclass(x), readable)
    if (is.character(x$uri) && length(x$uri) == 1) {
      x$uri <- shorthand_uri(x$uri)
    }
    if (is.character(x[["aliases"]])) x$aliases <- I(x[["aliases"]])
    x
  }
  readable(assembly)
}

# The adapter an alias file is read with, by its name as the engine's
# refNameAliasAdapterType() names it
alias_adapter <- function(aliases) {
  ncbi <- grepl("sequence_report\\.tsv(\\.gz)?$", bare_name(aliases),
    ignore.case = TRUE)
  list(type = if (ncbi) "NcbiSequenceReportAliasAdapter"
    else "RefNameAliasAdapter", location = location(aliases))
}

# Whether an assembly config names a sequence
has_sequence <- function(config) {
  !is.null(config$uri) || !is.null(config$sequence)
}
