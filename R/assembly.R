#' An assembly: the reference sequence and its contigs' aliases
#'
#' The assembly a JBrowse config's `assemblies` entry writes, for files that
#' spell one contig several ways (`chr7` in a BAM, `NC_000007.13` in a GFF,
#' `7` in the FASTA). Each reader takes it as `assembly`, as each track of a
#' JBrowse config names its assembly:
#' * a region may name a contig by any of its aliases, and the engine reads
#'   each file under the alias that file uses, as JBrowse renames a track's
#'   regions
#' * the reader names the contig in its tables as the region does, so
#'   [jbrowse_stack()] draws files that spell it differently on one axis
#' * its sequence is the `reference` where a reader takes one and none is
#'   given, which a CRAM and the bisulfite layer read against
#'
#' A JBrowse config's own assembly, as a list, serves the same:
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
    x$uri %||% x$sequence$adapter$uri %||% "none",
    x$refNameAliases$uri %||% "none"))
  invisible(x)
}

#' The contigs of an assembly
#'
#' Each contig the assembly's sequence lists, named as its file names it,
#' from 0 to its length: a `region` for any reader, to read a file whole.
#' @param assembly a [jbrowse_assembly()], or an assembly as a JBrowse
#'   config's `assemblies` entry writes it
#' @return a data frame of `chrom`, `start` and `end`
#' @export
#' @examples
#' fa <- system.file("extdata", "gapdh_hg38.fa.gz", package = "ggjbrowse")
#' jbrowse_regions(jbrowse_assembly(fa))
jbrowse_regions <- function(assembly) {
  r <- call_engine("regions", list(assembly = assembly_config(assembly)))
  data.frame(chrom = as.character(unlist(r$refName)),
    start = as.numeric(unlist(r$start)), end = as.numeric(unlist(r$end)))
}

# An assembly config as the engine reads it: each `uri` in it a path or URL
# made readable, as a shorthand the engine expands as JBrowse's config does
assembly_config <- function(assembly) {
  if (is.null(assembly)) return(NULL)
  if (!is.list(assembly)) {
    stop("an assembly is a jbrowse_assembly() or a JBrowse assembly config ",
      "as a list", call. = FALSE)
  }
  readable <- function(x) {
    if (!is.list(x)) return(x)
    x <- lapply(unclass(x), readable)
    if (is.character(x$uri) && length(x$uri) == 1) {
      x$uri <- shorthand_uri(x$uri)
    }
    x
  }
  readable(assembly)
}

# Whether an assembly config names a sequence
has_sequence <- function(assembly) {
  !is.null(assembly$uri) || !is.null(assembly$sequence)
}

# The engine args a reader adds for its reference: `assembly`, and the
# sequence adapter where `reference` names one
reference_args <- function(reference, assembly) {
  args <- list()
  if (!is.null(assembly)) args$assembly <- assembly_config(assembly)
  if (!is.null(reference)) args$sequenceAdapter <- sequence_adapter(reference)
  args
}
