LOCUS <- "^(.+):([0-9,]+)(-|\\.\\.)([0-9,]+)$"

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
  m <- regmatches(region, regexec(LOCUS, region))
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

# Each region's contig as `reference` spells it: `reference_chrom` for every
# region, or one per region, else the region's own
reference_chroms <- function(reference_chrom, regions) {
  if (is.null(reference_chrom)) return(regions$chrom)
  if (!length(reference_chrom) %in% c(1, nrow(regions))) {
    stop("give reference_chrom once, or once per region (", nrow(regions),
      ")", call. = FALSE)
  }
  rep_len(as.character(reference_chrom), nrow(regions))
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

# Read each region with read(chrom, start, end, share, i), `share` being the
# region's fraction of the regions' bp, and answer one result (join_regions())
read_regions <- function(region, start, end, read) {
  regions <- regions_of(region, start, end)
  w <- regions$end - regions$start
  join_regions(regions, lapply(seq_len(nrow(regions)), function(i)
    read(regions$chrom[i], regions$start[i], regions$end[i], w[i] / sum(w), i)))
}

# One result from each region's: the region's own where there is one, else
# every region's bound. Every row of a table with a position names its contig
# and region. Warns where no region's contig is one the file has, which is a
# misspelt contig rather than a region that holds nothing.
join_regions <- function(regions, parts) {
  parts <- lapply(seq_along(parts), function(i)
    placed(parts[[i]], regions$chrom[i], i))
  missing <- vapply(parts, function(p) p$missing_contig %||% "", "")
  if (all(missing != "")) warning(missing[1], call. = FALSE)
  if (length(parts) == 1) parts[[1]] else bind_regions(region_layout(regions),
    parts)
}

# The tables that list a result's keys, sections and rows rather than marks
# on the reference. bind_regions() binds the first three, each row once.
LISTING <- c("legend", "arc_legend", "info_fields", "sections", "sources",
  "rows", "rules", "bands")

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
    counts[["run"]] <- counts[["run"]] + max(c(-1, parts[[i]]$marks$run)) + 1
    parts[[i]] <- x
  }
  bound <- function(pick) {
    first <- pick(parts[[1]])
    frames <- lapply(parts, pick)
    columns <- unique(unlist(lapply(frames, names)))
    d <- do.call(rbind, lapply(frames, function(f) {
      for (column in setdiff(columns, names(f))) f[[column]] <- rep(NA, nrow(f))
      f[columns]
    }))
    rownames(d) <- NULL
    for (a in setdiff(names(attributes(first)),
      c("names", "row.names", "class"))) {
      attr(d, a) <- attr(first, a)
    }
    d
  }
  out <- parts[[1]]
  for (name in names(out)) {
    if (is.data.frame(out[[name]]) && !name %in% LISTING[-(1:3)]) {
      out[[name]] <- bound(function(x) x[[name]])
      if (name %in% LISTING) out[[name]] <- unique(out[[name]])
    }
  }
  for (key in names(out$keys)) {
    out$keys[[key]] <- unique(bound(function(x) x$keys[[key]]))
  }
  if (!is.null(out$domain)) out$domain <- range(lapply(parts, `[[`, "domain"))
  if (!is.null(out$coverage_scale)) {
    domain <- range(lapply(parts, function(p) p$coverage_scale$domain))
    out$coverage_scale$domain <- domain
    out$coverage_scale$ticks <- unique(c(domain[1], domain[2]))
  }
  if (!is.null(out$height)) out$height <- max(sapply(parts, `[[`, "height"))
  out$missing_contig <- ""
  out$region <- list(chrom = "", start = min(regions$cum_start),
    end = max(regions$cum_end))
  out$regions <- regions
  out
}

# A result's tables as its geoms draw them, on the axis that lays several
# regions end to end. One region's result draws as it is.
on_axis <- function(x) {
  if (inherits(x, "jbrowse_synteny")) return(synteny_on_axis(x))
  if (inherits(x, "jbrowse_dotplot")) return(dotplot_on_axis(x))
  if (inherits(x, "jbrowse_multiway")) return(multiway_on_axis(x))
  if (inherits(x, "jbrowse_variants")) x <- in_row_order(x)
  if (inherits(x, "jbrowse_wiggle")) x <- in_source_order(x)
  r <- x$regions
  if (is.null(r)) return(x)
  for (name in names(x)) {
    if (!is.data.frame(x[[name]])) next
    x[[name]] <- if (name %in% c("arcs", "junctions"))
      arcs_on_axis(x[[name]], r) else frame_on_axis(x[[name]], r)
  }
  x
}

# A read connection or splice junction with each foot on the axis by the
# region holding it, so a split read or splice across two regions arcs over
# the gap between them, and a foot outside every region at its own region's
# edge. Each region's fetch reports
# a connection reaching the other, so it is kept once.
arcs_on_axis <- function(d, regions) {
  if (is.null(d$.region) || !is.null(d$.on_axis)) return(d)
  r <- regions
  place <- function(pos) {
    i <- vapply(seq_along(pos), function(k) {
      hit <- which(r$chrom == d$chrom[k] & pos[k] >= r$start &
        pos[k] <= r$end)
      if (length(hit)) hit[1] else d$.region[k]
    }, 0)
    pmin(pmax(pos, r$start[i]), r$end[i]) - r$start[i] + r$offset[i]
  }
  d$start <- place(d$start)
  d$end <- place(d$end)
  d <- d[!duplicated(d[setdiff(names(d), ".region")]), , drop = FALSE]
  d$.on_axis <- rep(TRUE, nrow(d))
  d
}

# A table's rows on the axis `regions` lays end to end: each moved by its
# region's offset, a span cut to its region and a mark outside it dropped, so
# nothing draws over the gap into the next region. A row keeps `.region` and
# gains `.on_axis`, so a table already moved stays where it is.
frame_on_axis <- function(d, regions) {
  if (is.null(regions) || is.null(d$.region) || !is.null(d$.on_axis)) {
    return(d)
  }
  r <- regions
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
  d$.on_axis <- rep(TRUE, nrow(d))
  d
}
