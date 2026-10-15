LOCUS <- "^(.+):([0-9,]+)(-|\\.\\.)([0-9,]+)$"

# A read with no region reads the whole file: each contig the file lists, from
# 0 to the end its index addresses, refused past WHOLE_FILE_ROWS rows in any
# one table, which a region read has no need of
WHOLE_FILE_ROWS <- 5e6

file_regions <- function(adapter) {
  region_frame(call_engine("fileRegions", list(adapterConfig = adapter)))
}

# The engine's regions as `chrom`, `start` and `end`
region_frame <- function(r) {
  data.frame(chrom = as.character(r$refName), start = as.numeric(r$start),
    end = as.numeric(r$end))
}

# The regions a reader reads: those `region` names, or each contig the file
# lists where it names none
reader_regions <- function(adapter, region, start, end) {
  if (is.null(region)) file_regions(adapter) else regions_of(region, start, end)
}

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

# Read each region with read(chrom, start, end, i), and answer one result
# (join_regions())
read_regions <- function(regions, read) {
  join_regions(regions, lapply(seq_len(nrow(regions)), function(i)
    read(regions$chrom[i], regions$start[i], regions$end[i], i)))
}

# One result from each region's: the region's own where there is one, else
# every region's bound. Every row of a table with a position names its contig
# and region. Warns where no region's contig is one the file has, which is a
# misspelt contig rather than a region that holds nothing.
join_regions <- function(regions, parts) {
  parts <- lapply(seq_along(parts), function(i)
    placed(parts[[i]], regions$chrom[i], i))
  warn_missing(lapply(parts, `[[`, "missing_contig"))
  if (length(parts) == 1) parts[[1]] else bind_regions(region_layout(regions),
    parts)
}

# Each region's missing-contig message, "" or NULL where the file has the
# region's contig; warns with the first where no region's contig is the file's
warn_missing <- function(messages) {
  missing <- vapply(messages, function(m) m %||% "", "")
  if (all(missing != "")) warning(missing[1], call. = FALSE)
  missing
}

# A result the engine read over every region in one call, as join_regions()
# answers one: each row of a table the engine numbers by its 0-based `region`
# names its contig and region, and several regions lie on the axis
# region_layout() lays them on
in_one_call <- function(x, regions, missing_contigs) {
  missing <- warn_missing(lapply(seq_len(nrow(regions)), function(i)
    missing_contigs[[i]]))
  for (name in setdiff(names(x), LISTING)) {
    d <- x[[name]]
    if (!is.data.frame(d) || is.null(d$region)) next
    i <- as.integer(d$region) + 1L
    d$region <- NULL
    d$chrom <- regions$chrom[i]
    d$.region <- i
    x[[name]] <- d
  }
  if (nrow(regions) == 1) {
    x$region <- as.list(regions)
    x$missing_contig <- missing
    return(x)
  }
  on_layout(x, region_layout(regions))
}

# A result over several regions, its region the axis `regions` lays them on
on_layout <- function(x, regions) {
  x$region <- list(chrom = "", start = min(regions$cum_start),
    end = max(regions$cum_end))
  x$missing_contig <- ""
  x$regions <- regions
  x
}

# The tables that list a result's keys and fields rather than marks on the
# reference, which bind_regions() binds, each row once, and those every
# region shares, which it takes from the first
KEY_TABLES <- c("legend", "arc_legend", "info_fields", "format_fields",
  "attribute_fields", "annotation_fields")
SHARED_TABLES <- c("sections", "sources", "rows", "rules", "bands")
LISTING <- c(KEY_TABLES, SHARED_TABLES)

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
# genomic positions. A read or variant number indexes the bound table. A key
# lists each entry once, and the result's region is the axis region_layout()
# lays them on, which on_axis() moves them onto.
bind_regions <- function(regions, parts) {
  indexes <- c(read = "reads", variant = "variants")
  counts <- stats::setNames(numeric(length(indexes)), indexes)
  for (i in seq_along(parts)) {
    x <- parts[[i]]
    for (name in setdiff(names(x), LISTING)) {
      d <- x[[name]]
      if (!is.data.frame(d)) next
      for (column in intersect(names(indexes), names(d))) {
        d[[column]] <- d[[column]] + counts[[indexes[[column]]]]
      }
      x[[name]] <- d
    }
    for (table in intersect(indexes, names(x))) {
      counts[[table]] <- counts[[table]] + nrow(parts[[i]][[table]])
    }
    parts[[i]] <- x
  }
  bound <- function(pick) {
    first <- pick(parts[[1]])
    frames <- lapply(parts, pick)
    columns <- unique(unlist(lapply(frames, names)))
    nested <- Filter(function(column)
      any(vapply(frames, function(f) is.list(f[[column]]), NA)), columns)
    d <- do.call(rbind, lapply(frames, function(f) {
      for (column in setdiff(columns, names(f))) f[[column]] <- rep(NA, nrow(f))
      for (column in nested) f[[column]] <- as.list(f[[column]])
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
    if (is.data.frame(out[[name]]) && !name %in% SHARED_TABLES) {
      out[[name]] <- bound(function(x) x[[name]])
      if (name %in% KEY_TABLES) out[[name]] <- unique(out[[name]])
    }
  }
  for (key in names(out$keys)) {
    out$keys[[key]] <- unique(bound(function(x) x$keys[[key]]))
  }
  if (!is.null(out$coverage_scale)) {
    domain <- range(lapply(parts, function(p) p$coverage_scale$domain))
    out$coverage_scale$domain <- domain
    out$coverage_scale$ticks <- unique(c(domain[1], domain[2]))
  }
  on_layout(out, regions)
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
  same <- setdiff(names(d), c(".region", "region", "startRegion", "endRegion"))
  d <- d[!duplicated(d[same]), , drop = FALSE]
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
