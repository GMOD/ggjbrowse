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
# own where there is one, else every region's bound. Every row of a table
# with a position names its contig and region. Warns where no region's
# contig is one the file has, which is a misspelt contig rather than a region
# that holds nothing.
read_regions <- function(region, start, end, read) {
  regions <- regions_of(region, start, end)
  w <- regions$end - regions$start
  parts <- lapply(seq_len(nrow(regions)), function(i)
    placed(read(regions$chrom[i], regions$start[i], regions$end[i],
      w[i] / sum(w), i), regions$chrom[i], i))
  missing <- vapply(parts, function(p) p$missing_contig %||% "", "")
  if (all(missing != "")) warning(missing[1], call. = FALSE)
  if (length(parts) == 1) parts[[1]] else bind_regions(region_layout(regions),
    parts)
}

# The tables that list a result's keys, sections and rows rather than marks
# on the reference
LISTING <- c("legend", "arc_legend", "sections", "sources", "rows", "rules")

placed <- function(x, chrom, i) {
  for (name in setdiff(names(x), LISTING)) {
    d <- x[[name]]
    if (is.data.frame(d) && (!is.null(d$start) || !is.null(d$position))) {
      d$chrom <- rep(chrom, nrow(d))
      d$.region <- rep(i, nrow(d))
      x[[name]] <- d
    }
  }
  x
}

# Bind each region's result into one of the same class, its rows at their
# genomic positions. A read, feature, subfeature or variant number indexes
# the bound table, and a line's run stays its own. A key lists each entry
# once, the score domain spans every region's, and the result's region is
# the axis region_layout() lays them on, which on_axis() moves them onto.
bind_regions <- function(regions, parts) {
  indexes <- c(read = "reads", feature = "features",
    subfeature = "subfeatures", variant = "variants")
  counts <- stats::setNames(numeric(length(indexes) + 1), c(indexes, "run"))
  for (i in seq_along(parts)) {
    x <- parts[[i]]
    for (name in setdiff(names(x), LISTING)) {
      d <- x[[name]]
      if (!is.data.frame(d)) next
      for (column in intersect(names(indexes), names(d))) {
        d[[column]] <- d[[column]] + counts[[indexes[[column]]]]
      }
      if (!is.null(d$run)) d$run <- d$run + counts[["run"]]
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
    if (is.data.frame(out[[name]]) && !name %in% LISTING[-(1:2)]) {
      out[[name]] <- bound(function(x) x[[name]])
      if (name %in% LISTING) out[[name]] <- unique(out[[name]])
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

# A result's tables as its geoms draw them, on the axis that lays several
# regions end to end: each row moved by its region's offset, a span cut to
# its region and a mark outside it dropped, so nothing draws over the gap
# into the next region. One region's result draws as it is, and a result
# already on its axis stays there.
on_axis <- function(x) {
  r <- x$regions
  if (is.null(r) || isTRUE(x$on_axis)) return(x)
  x$on_axis <- TRUE
  for (name in names(x)) {
    d <- x[[name]]
    if (!is.data.frame(d) || is.null(d$.region)) next
    if (!is.null(d$start)) {
      d <- d[d$end > r$start[d$.region] & d$start < r$end[d$.region], ,
        drop = FALSE]
      d$start <- pmax(d$start, r$start[d$.region])
      d$end <- pmin(d$end, r$end[d$.region])
    }
    if (!is.null(d$position)) {
      d <- d[d$position >= r$start[d$.region] &
        d$position < r$end[d$.region], , drop = FALSE]
    }
    shift <- r$offset[d$.region] - r$start[d$.region]
    if (!is.null(d$start)) {
      d$start <- d$start + shift
      d$end <- d$end + shift
    }
    if (!is.null(d$position)) d$position <- d$position + shift
    x[[name]] <- d
  }
  x
}
