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
#' a$region$chrom
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

#' The contigs of an assembly
#'
#' Each contig the assembly's sequence lists, named as its file names it,
#' from 0 to its length: a `region` for any reader, to read a file whole.
#' @param assembly a [jbrowse_assembly()], an assembly as a JBrowse config's
#'   `assemblies` entry writes it, or the sequence's path or URL alone
#' @return a data frame of `chrom`, `start` and `end`
#' @export
#' @examples
#' fa <- system.file("extdata", "gapdh_hg38.fa.gz", package = "ggjbrowse")
#' jbrowse_regions(jbrowse_assembly(fa))
jbrowse_regions <- function(assembly) {
  region_frame(call_engine("regions",
    list(assembly = assembly_config(assembly))))
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
  ncbi <- grepl("sequence_report\\.tsv(\\.gz)?$", sub("\\?.*$", "", aliases),
    ignore.case = TRUE)
  list(type = if (ncbi) "NcbiSequenceReportAliasAdapter"
    else "RefNameAliasAdapter", location = location(aliases))
}

# Whether an assembly config names a sequence
has_sequence <- function(config) {
  !is.null(config$uri) || !is.null(config$sequence)
}
