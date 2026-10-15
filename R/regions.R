LOCUS <- "^(.+):([0-9,]+)(-|\\.\\.)([0-9,]+)$"

# A read with no region reads the whole file: each contig the file lists, to
# its length where the file or `assembly` states one, refused past
# WHOLE_FILE_ROWS rows in any one table, which a region read has no need of.
# A contig the index shows empty is left out; one it cannot tell about stays.
WHOLE_FILE_ROWS <- 5e6

# Each contig a file lists, from 0 to the length the file or `assembly`
# states, and whether the file's index shows records there: `NA` where it
# cannot tell
file_contigs <- function(adapter, assembly = NULL) {
  args <- list(adapterConfig = adapter)
  args$assembly <- assembly
  r <- call_engine("fileRegions", args)
  d <- region_frame(r)
  d$has_data <- as.logical(unlist(r$hasData) %||% logical(0))
  d
}

# The engine's regions as `chrom`, `start` and `end`
region_frame <- function(r) {
  data.frame(chrom = as.character(r$refName), start = as.numeric(r$start),
    end = as.numeric(r$end))
}

# The regions a reader reads, 0-based half-open: those `region` names, or
# each contig the file lists but those its index shows empty where it names
# none. A whole contig, or the whole file, is marked `whole`, which the
# engine caps at WHOLE_FILE_ROWS rows.
reader_regions <- function(adapter, region, assembly = NULL) {
  contigs <- memo(function() file_contigs(adapter, assembly))
  if (is.null(region)) {
    all <- contigs()
    all <- all[!all$has_data %in% FALSE, c("chrom", "start", "end")]
    rownames(all) <- NULL
    return(structure(all, whole = TRUE))
  }
  regions_of(region, function() {
    file <- contigs()
    if (is.null(assembly)) return(file)
    rbind(file[c("chrom", "start", "end")], region_frame(call_engine("regions",
      list(assembly = assembly))))
  })
}

# Regions as the engine's `regions` argument names them
engine_regions <- function(regions) {
  data.frame(refName = regions$chrom, start = regions$start, end = regions$end)
}

# Regions as 0-based half-open `chrom`, `start` and `end`, from loci as
# JBrowse's location box writes them ("ctgA:1,001-9,000", 1-based and
# inclusive), contig names alone for the whole contig at the length
# `contigs()` lists, a data frame of `chrom`, `start` and `end`, or a GRanges
regions_of <- function(region, contigs = NULL) {
  if (inherits(region, "GRanges")) {
    d <- as.data.frame(region)
    return(data.frame(chrom = as.character(d$seqnames), start = d$start - 1,
      end = as.numeric(d$end)))
  }
  if (is.data.frame(region)) return(region_table(region))
  if (is.factor(region)) region <- as.character(region)
  if (!is.character(region) || !length(region)) {
    rlang::abort(paste("a region is a locus such as \"ctgA:1,001-9,000\", a",
      "contig name, a data frame of chrom, start and end (0-based,",
      "half-open), or a GRanges"), class = "ggjbrowse_bad_region")
  }
  m <- regmatches(region, regexec(LOCUS, region))
  locus <- lengths(m) > 0
  number <- function(i) as.numeric(gsub(",", "", vapply(m[locus], `[`, "", i)))
  out <- data.frame(chrom = region, start = NA_real_, end = NA_real_)
  out$chrom[locus] <- vapply(m[locus], `[`, "", 2)
  out$start[locus] <- number(3) - 1
  out$end[locus] <- number(5)
  if (all(locus)) return(out)
  known <- if (is.null(contigs)) NULL else contigs()
  at <- match(region[!locus], known$chrom)
  if (anyNA(at)) {
    rlang::abort(sprintf(paste("\"%s\" is neither a locus such as",
      "\"ctgA:1,001-9,000\" nor a contig of the file%s"),
      region[!locus][is.na(at)][1], if (length(known$chrom)) paste0(
        ", whose contigs are ", paste(utils::head(unique(known$chrom), 20),
          collapse = ", ")) else ""), class = "ggjbrowse_bad_region")
  }
  out$start[!locus] <- known$start[at]
  out$end[!locus] <- known$end[at]
  structure(out, whole = TRUE)
}

# A data frame of regions as `chrom`, `start` and `end`, refused where it
# reads as 1-based: a GRanges' `seqnames`, or a `width` of end - start + 1
region_table <- function(region) {
  width <- if ("width" %in% names(region)) region[["width"]]
  closed <- !is.null(width) && nrow(region) &&
    isTRUE(all(width == region[["end"]] - region[["start"]] + 1))
  if ("seqnames" %in% names(region) || closed) {
    rlang::abort(paste0("a region frame is 0-based and half-open, chrom, ",
      "start and end, as BED writes it, but this one reads as 1-based (",
      if (closed) "its width is end - start + 1" else "it has seqnames",
      "): pass the GRanges itself, or subtract 1 from start"),
      class = "ggjbrowse_bad_region")
  }
  missing <- setdiff(c("chrom", "start", "end"), names(region))
  if (length(missing)) {
    rlang::abort(paste0("a region frame names chrom, start and end, ",
      "0-based and half-open; this one lacks ", paste(missing,
        collapse = ", ")), class = "ggjbrowse_bad_region")
  }
  data.frame(chrom = as.character(region[["chrom"]]),
    start = as.numeric(region[["start"]]), end = as.numeric(region[["end"]]))
}

# Several regions end to end on one axis: each spans its width in bp, a gap
# of `gap_frac` of the axis between them; the first keeps its own
# coordinates, so one region is drawn at its genomic position. JBrowse's
# views lay regions out with no gap and paint a 3 px bar over each seam; a
# stack keeps a gap, so labels either side of it stay apart, and paints the
# bar in it (genomic_x()), while a comparative axis lays its regions as the
# synteny and dotplot workers do, with none.
region_layout <- function(regions, gap_frac = 0.012) {
  w <- regions$end - regions$start
  gap <- if (gap_frac > 0) max(1, round(sum(w) * gap_frac)) else 0
  regions$offset <- regions$start[1] + utils::head(cumsum(c(0, w + gap)),
    length(w))
  regions$cum_start <- regions$offset
  regions$cum_end <- regions$offset + w
  regions
}

# Each region's missing-contig message, "" or NULL where the file has the
# region's contig: an error where no region's contig is the file's, and else
# one warning naming the contigs that read empty
check_missing <- function(messages, chrom) {
  missing <- vapply(messages, function(m) m %||% "", "")
  if (length(missing) && all(missing != "")) {
    rlang::abort(missing[1], class = "ggjbrowse_missing_contig", call = NULL)
  }
  lacking <- unique(chrom[missing != ""])
  said <- if (length(lacking) == 1) missing[missing != ""][1]
    else if (length(lacking) > 1) {
      sprintf(paste("%d contigs are not the file's and read empty: %s;",
        "regions$missing says why"), length(lacking), paste(c(
          utils::head(lacking, 5), if (length(lacking) > 5)
            sprintf("and %d more", length(lacking) - 5)), collapse = ", "))
    }
  if (length(said)) warning(said, call. = FALSE)
  structure(missing, warning = said)
}

# A result the engine `t` read over every region in one call, laid out for
# `width` px: each row of a table the engine numbers by its 0-based `region`
# gets its 1-based `region`, and its contig where the row has coordinates of
# its own, `regions` lays the regions on one axis, each with its
# missing-contig message, and `warnings` holds what the engine warned of
in_one_call <- function(x, regions, t, width) {
  missing <- check_missing(lapply(seq_len(nrow(regions)), function(i)
    t$missingContigs[[i]]), regions$chrom)
  cls <- class(x)
  x <- unclass(x)
  for (name in names(x)) {
    d <- x[[name]]
    if (!is.data.frame(d) || is.null(d$region)) next
    i <- as.integer(d$region) + 1L
    d$region <- NULL
    if (!is.null(d$start) || !is.null(d$position)) {
      d$chrom <- regions$chrom[i]
    }
    d$region <- i
    x[[name]] <- d
  }
  layout <- region_layout(regions[c("chrom", "start", "end")])
  layout$missing <- as.vector(missing)
  x$regions <- layout
  x$width <- width
  x$warnings <- unique(c(character(0), unlist(t$warnings),
    attr(missing, "warning")))
  structure(x, class = cls)
}

# The x range a result's regions span on its axis
xlim_of <- function(x) c(min(x$regions$cum_start), max(x$regions$cum_end))

# A result's tables as its geoms draw them, on the axis that lays several
# regions end to end, each table marked so it moves once. One region's
# tables draw where they are.
on_axis <- function(x) {
  if (inherits(x, "jbrowse_synteny")) return(synteny_on_axis(x))
  if (inherits(x, "jbrowse_dotplot")) return(dotplot_on_axis(x))
  if (inherits(x, "jbrowse_multiway")) return(multiway_on_axis(x))
  if (inherits(x, "jbrowse_variants")) x <- in_row_order(x)
  if (inherits(x, "jbrowse_wiggle")) x <- in_source_order(x)
  cls <- class(x)
  families <- families_of(x)
  x <- unclass(x)
  r <- x$regions
  for (name in names(x)) {
    if (!is.data.frame(x[[name]]) || name == "regions") next
    x[[name]] <- if (name %in% c("arcs", "junctions"))
      arcs_on_axis(x[[name]], r) else frame_on_axis(x[[name]], r)
  }
  structure(joined(x, families), class = cls, on_axis = TRUE)
}

is_on_axis <- function(d) isTRUE(attr(d, "on_axis"))

# Whether `d` names its rows' regions of `regions` in a `region` column
has_regions <- function(d, regions) {
  is.numeric(d$region) && all(d$region %in% seq_len(nrow(regions)))
}

# A read connection or splice junction with each foot on the axis in the
# region the engine names it by (`startRegion`, `endRegion`), so one across
# two regions arcs over the gap between them, and a foot outside every region
# at the edge of the connection's own region. A connection between two
# contigs has no genomic radius, so its dome rises half as high as its feet
# lie apart on the axis, as the display draws it.
arcs_on_axis <- function(d, regions) {
  if (is_on_axis(d) || !has_regions(d, regions)) return(d)
  if (nrow(regions) == 1) return(structure(d, on_axis = TRUE))
  r <- regions
  place <- function(pos, foot) {
    i <- ifelse(foot >= 0, foot + 1L, d$region)
    pmin(pmax(pos, r$start[i]), r$end[i]) - r$start[i] + r$offset[i]
  }
  d$start <- place(d$start, d$startRegion)
  d$end <- place(d$end, d$endRegion)
  if (!is.null(d$y)) {
    feet <- d$startRegion >= 0 & d$endRegion >= 0
    across <- feet & r$chrom[pmax(d$startRegion, 0) + 1L] !=
      r$chrom[pmax(d$endRegion, 0) + 1L]
    d$y[across] <- abs(d$end[across] - d$start[across]) / 2
  }
  structure(d, on_axis = TRUE)
}

# A table's rows on the axis `regions` lays end to end: each moved by its
# region's offset, a span cut to its region and a mark outside it dropped, so
# nothing draws over the gap into the next region
frame_on_axis <- function(d, regions) {
  if (is_on_axis(d) || !has_regions(d, regions)) return(d)
  if (nrow(regions) == 1) return(structure(d, on_axis = TRUE))
  r <- regions
  if (!is.null(d$start)) {
    d <- d[d$end > r$start[d$region] & d$start < r$end[d$region], ,
      drop = FALSE]
    d$start <- pmax(d$start, r$start[d$region])
    d$end <- pmin(d$end, r$end[d$region])
  }
  if (!is.null(d$position)) {
    d <- d[d$position >= r$start[d$region] &
      d$position < r$end[d$region], , drop = FALSE]
  }
  shift <- r$offset[d$region] - r$start[d$region]
  if (!is.null(d$start)) {
    d$start <- d$start + shift
    d$end <- d$end + shift
  }
  if (!is.null(d$position)) d$position <- d$position + shift
  structure(d, on_axis = TRUE)
}
