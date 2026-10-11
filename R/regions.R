# The regions a reader was asked for, 0-based half-open: a contig with
# `start` and `end`, locus strings as JBrowse's location box writes them
# ("ctgA:1,001-9,000", 1-based and inclusive), or a data frame of chrom,
# start and end
regions_of <- function(region, start, end) {
  if (is.data.frame(region)) {
    return(data.frame(chrom = as.character(region$chrom),
      start = as.numeric(region$start), end = as.numeric(region$end)))
  }
  if (!is.null(start) || !is.null(end)) {
    if (is.null(start) || is.null(end)) {
      stop("give both start and end, or a locus such as \"ctgA:1,001-9,000\"",
        call. = FALSE)
    }
    return(data.frame(chrom = as.character(region), start = as.numeric(start),
      end = as.numeric(end)))
  }
  m <- regmatches(region, regexec(
    "^(.+):([0-9,]+)(-|\\.\\.)([0-9,]+)$", region))
  bad <- lengths(m) == 0
  if (any(bad)) {
    stop("cannot read ", paste0("\"", region[bad], "\"", collapse = ", "),
      " as a locus: write \"ctgA:1,001-9,000\", or pass start and end",
      call. = FALSE)
  }
  number <- function(i) as.numeric(gsub(",", "", vapply(m, `[`, "", i)))
  data.frame(chrom = vapply(m, `[`, "", 2), start = number(3) - 1,
    end = number(5))
}

# Several regions end to end on one axis, as a JBrowse view lays them out:
# each spans its width in bp, a gap between them; the first keeps its own
# coordinates, so one region is drawn at its genomic position
region_layout <- function(regions, gap_frac = 0.012) {
  w <- regions$end - regions$start
  gap <- max(1, round(sum(w) * gap_frac))
  regions$offset <- regions$start[1] + utils::head(cumsum(c(0, w + gap)),
    length(w))
  regions$cum_start <- regions$offset
  regions$cum_end <- regions$offset + w
  regions
}

# Read each region with read(chrom, start, end, share, i), `share` being the
# region's fraction of the regions' bp, and answer one result: the region's
# own where there is one, else every region's bound onto one axis. Warns where
# no region's contig is one the file has, which is a misspelt contig rather
# than a region that holds nothing.
read_regions <- function(region, start, end, read) {
  regions <- regions_of(region, start, end)
  w <- regions$end - regions$start
  parts <- lapply(seq_len(nrow(regions)), function(i)
    read(regions$chrom[i], regions$start[i], regions$end[i], w[i] / sum(w), i))
  missing <- vapply(parts, function(p) p$missing_contig %||% "", "")
  if (all(missing != "")) warning(missing[1], call. = FALSE)
  if (length(parts) == 1) parts[[1]] else bind_regions(region_layout(regions),
    parts)
}

# Bind each region's result into one of the same class on the cumulative
# axis region_layout() defines, which its geoms draw end to end. A table's
# start, end and position move by their region's offset and each row gains
# `.region`; a span is cut to its region and a mark outside it dropped, so
# nothing draws over the gap into the next region. A read, feature,
# subfeature or variant number indexes the bound table, and a line's run
# stays its own. A key lists each entry once, the score domain spans every
# region's, and the result's region is the whole axis.
bind_regions <- function(regions, parts) {
  indexes <- c(read = "reads", feature = "features",
    subfeature = "subfeatures", variant = "variants")
  listing <- c("legend", "arc_legend", "sections", "sources", "rows", "rules")
  counts <- stats::setNames(numeric(length(indexes) + 1), c(indexes, "run"))
  for (i in seq_along(parts)) {
    lo <- regions$start[i]
    hi <- regions$end[i]
    shift <- regions$offset[i] - lo
    x <- parts[[i]]
    for (name in setdiff(names(x), listing)) {
      d <- x[[name]]
      if (!is.data.frame(d)) next
      if (!is.null(d$start)) {
        d <- d[d$end > lo & d$start < hi, , drop = FALSE]
        d$start <- pmax(d$start, lo) + shift
        d$end <- pmin(d$end, hi) + shift
      }
      if (!is.null(d$position)) {
        d <- d[d$position >= lo & d$position < hi, , drop = FALSE]
        d$position <- d$position + shift
      }
      for (column in intersect(names(indexes), names(d))) {
        d[[column]] <- d[[column]] + counts[[indexes[[column]]]]
      }
      if (!is.null(d$run)) d$run <- d$run + counts[["run"]]
      d$.region <- rep(i, nrow(d))
      x[[name]] <- d
    }
    for (table in intersect(indexes, names(x))) {
      counts[[table]] <- counts[[table]] + nrow(parts[[i]][[table]])
    }
    counts[["run"]] <- counts[["run"]] + max(c(0, parts[[i]]$marks$run))
    parts[[i]] <- x
  }
  bound <- function(pick) {
    first <- pick(parts[[1]])
    d <- do.call(rbind, lapply(parts, pick))
    rownames(d) <- NULL
    for (a in setdiff(names(attributes(first)),
      c("names", "row.names", "class"))) {
      attr(d, a) <- attr(first, a)
    }
    d
  }
  out <- parts[[1]]
  for (name in names(out)) {
    if (is.data.frame(out[[name]]) && !name %in% listing[-(1:2)]) {
      out[[name]] <- bound(function(x) x[[name]])
      if (name %in% listing) out[[name]] <- unique(out[[name]])
    }
  }
  for (key in names(out$keys)) {
    out$keys[[key]] <- unique(bound(function(x) x$keys[[key]]))
  }
  if (!is.null(out$domain)) out$domain <- range(lapply(parts, `[[`, "domain"))
  if (!is.null(out$height)) out$height <- max(sapply(parts, `[[`, "height"))
  out$missing_contig <- ""
  out$region <- list(chrom = "", start = min(regions$cum_start),
    end = max(regions$cum_end))
  out$regions <- regions
  out
}
