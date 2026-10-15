# The engine is @jbrowse/tables' bare-host bundle, run in one V8 context per
# session. Every byte it reads comes through read_range(), which JavaScript
# calls synchronously; jsonlite sends the raw body back base64-encoded.

read_range <- function(url, start, end) {
  if (startsWith(url, "file://")) {
    path <- utils::URLdecode(sub("\\?.*$", "", sub("^file://", "", url)))
    path <- sub("^/([A-Za-z]:)", "\\1", path)
    size <- file.size(path)
    last <- if (end < 0) size - 1 else min(end, size - 1)
    con <- file(path, "rb")
    on.exit(close(con))
    seek(con, start)
    bytes <- readBin(con, "raw", max(0, last - start + 1))
    return(list(status = if (end < 0 && start == 0) 200L else 206L,
      size = size, bytes = bytes))
  }
  h <- curl::new_handle()
  if (end >= 0 || start > 0) {
    curl::handle_setheaders(h, Range = sprintf("bytes=%.0f-%s", start,
      if (end < 0) "" else sprintf("%.0f", end)))
  }
  r <- curl::curl_fetch_memory(url, handle = h)
  range <- curl::parse_headers_list(r$headers)[["content-range"]]
  size <- if (is.null(range)) length(r$content) else
    as.numeric(sub(".*/", "", range))
  list(status = r$status_code, size = size, bytes = r$content)
}

engine <- local({
  ctx <- NULL
  function(reset = FALSE) {
    if (reset) {
      ctx <<- NULL
      deferred$warnings <- character(0)
      return(invisible())
    }
    if (is.null(ctx)) {
      fresh <- V8::v8()
      fresh$eval(paste(
        "var jbrowseHost = { readRange: function(url, start, end) {",
        "return console.r.call('ggjbrowse:::read_range_deferring',",
        "{ url: url, start: start, end: end }) } };",
        "var ggjbrowseQueued = { warnings: [], errors: [] };",
        "function ggjbrowseQueue(kind) { return function() {",
        "ggjbrowseQueued[kind].push(Array.prototype.join.call(arguments, ' '))",
        "} }",
        "console.warn = ggjbrowseQueue('warnings');",
        "console.error = ggjbrowseQueue('errors')"))
      fresh$source(system.file("js", "jbrowse-tables.js", package = "ggjbrowse"))
      ctx <<- fresh
    }
    ctx
  }
})

#' Start the engine afresh
#'
#' The engine caches what it reads, keyed by each file's URL. A local file's
#' URL carries its modification time and size, so an overwritten file is read
#' again; a remote file changed on its server is not, until this is called.
#' @export
jbrowse_reset <- function() engine(reset = TRUE)

# V8 raises console.warn as an R warning and console.error as an R error from
# inside its own frames, and a handler that unwinds through them (tryCatch,
# options(warn = 2)) leaves the context broken, so the next call aborts R.
# The engine's warnings and errors, and read_range()'s warnings, wait in a
# queue that in_engine() raises once its call has returned: an error stops
# the call, as it did when V8 raised it.
deferred <- new.env()
deferred$warnings <- character(0)

read_range_deferring <- function(url, start, end) {
  withCallingHandlers(read_range(url, start, end), warning = function(w) {
    deferred$warnings <- c(deferred$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}

# Runs use(<the engine's context>) and raises what the engine queued meanwhile
in_engine <- function(use) {
  ctx <- engine()
  on.exit({
    queued <- ctx$get("ggjbrowseQueued")
    ctx$eval("ggjbrowseQueued = { warnings: [], errors: [] }")
    warnings <- unique(c(deferred$warnings, unlist(queued$warnings)))
    deferred$warnings <- character(0)
    for (w in warnings) warning(w, call. = FALSE)
    if (length(queued$errors)) {
      stop(paste(unique(unlist(queued$errors)), collapse = "\n"),
        call. = FALSE)
    }
  })
  use(ctx)
}

as_engine_json <- function(args) {
  objects <- function(x) {
    if (!is.list(x) || inherits(x, "AsIs")) return(x)
    if (!length(x)) return(structure(list(), names = character()))
    x[] <- lapply(x, objects)
    x
  }
  as.character(jsonlite::toJSON(objects(args), auto_unbox = TRUE))
}

call_engine <- function(table, args) {
  wire <- in_engine(function(ctx) {
    json <- ctx$call(sprintf("jbrowse.%s", table), as_engine_json(args),
      await = TRUE)
    bytes <- ctx$get("jbrowseBytes")
    ctx$eval("jbrowseBytes = null")
    list(json = json, bytes = bytes)
  })
  bytes <- rawConnection(wire$bytes)
  on.exit(close(bytes))
  from_wire(jsonlite::fromJSON(wire$json, simplifyVector = TRUE,
    simplifyDataFrame = FALSE), bytes)
}

# The engine's tables as R values. The engine writes the numbers of a column
# as bytes, which `bytes`, a connection, holds: `$i32`, `$f64` and `$lgl`
# name where in it a column's values start and how many there are, and NA
# crosses as R's own NA. Repeated text comes as `$dict`, each row's code in
# its `levels`, and a nested column as `$list`, its rows' `lengths` (-1 for a
# missing row) and their `values` end to end, which reads as a list column, a
# missing row NA.
from_wire <- function(x, bytes) {
  if (!is.list(x)) return(x)
  kind <- names(x)
  if (length(kind) == 1 && startsWith(kind, "$")) {
    spec <- x[[1]]
    return(switch(kind,
      "$i32" = read_column(bytes, spec, "integer", 4L),
      "$f64" = read_column(bytes, spec, "double", 8L),
      "$lgl" = as.logical(read_column(bytes, spec, "integer", 4L)),
      "$dict" = dict_column(spec$levels, from_wire(spec$codes, bytes)),
      "$list" = list_column(from_wire(spec$lengths, bytes),
        from_wire(spec$values, bytes)),
      stop("the engine sent a column as ", kind, ", which this ggjbrowse ",
        "cannot read", call. = FALSE)))
  }
  x[] <- lapply(x, from_wire, bytes)
  x
}

read_column <- function(bytes, at, what, size) {
  seek(bytes, at[[1]])
  readBin(bytes, what, n = at[[2]], size = size, endian = "little")
}

dict_column <- function(levels, codes) {
  codes[codes < 0L] <- NA_integer_
  as.character(unlist(levels))[codes + 1L]
}

list_column <- function(lengths, values) {
  lengths <- as.integer(unlist(lengths))
  if (is.list(values)) values <- unlist(values) %||% logical(0)
  rows <- split(values, factor(rep(seq_along(lengths), pmax(lengths, 0L)),
    levels = seq_along(lengths)))
  names(rows) <- NULL
  rows[lengths < 0] <- list(values[NA_integer_])
  rows
}

# A file's location for the engine: a URL as given, and a local path as a
# file:// URL, each segment percent-encoded so a %, ? or # in a name stays in
# the name, with the file's modification time and size as a query, which
# read_range() drops, so a file overwritten in place is a new cache key.
location <- function(path) {
  uri <- if (grepl("^[a-z]+://", path)) path else {
    full <- normalizePath(path, winslash = "/", mustWork = TRUE)
    info <- file.info(full)
    encoded <- paste(vapply(strsplit(full, "/", fixed = TRUE)[[1]],
      utils::URLencode, "", reserved = TRUE, repeated = TRUE), collapse = "/")
    sprintf("file://%s%s?v=%s-%.0f", if (startsWith(encoded, "/")) "" else "/",
      encoded, format(as.numeric(info$mtime), digits = 17), info$size)
  }
  list(uri = uri, locationType = "UriLocation")
}

# The adapter config of the first of `formats` whose `pattern` matches
# `file`'s name, any URL query aside: its `type`, the file's location in its
# `slot`, and where it names an `index` type, the index (`index`, else the
# file's name and that type's extension), CSI where its name says so
guess_adapter <- function(file, index, formats, expected = NULL) {
  name <- tolower(sub("\\?.*$", "", file))
  for (format in formats) {
    if (!grepl(format$pattern, name)) next
    adapter <- list(type = format$type)
    adapter[[format$slot]] <- location(file)
    if (!is.null(format$index)) {
      index <- index %||% beside(file, paste0(".", tolower(format$index)))
      adapter$index <- list(location = location(index),
        indexType = if (grepl("\\.csi(\\?|$)", index)) "CSI"
          else format$index)
    }
    return(adapter)
  }
  stop("cannot tell the format of ", file, ": expected ", expected,
    call. = FALSE)
}

# The file beside `file` its name and `suffix` name, a URL's query kept after
# it, as a signed URL's index carries the same query
beside <- function(file, suffix) {
  sub("^([^?]*)", paste0("\\1", suffix), file)
}

# The adapter config the engine reads `file` through: a track's own adapter
# config as given, each of its locations a path or URL made readable, or else
# the one `guess` makes of a path or URL
adapter_of <- function(file, guess) {
  if (!is.list(file)) return(guess())
  readable <- function(x) {
    if (!is.list(x)) return(x)
    path <- x[["uri"]] %||% x[["localPath"]]
    single <- is.character(path) && length(path) == 1
    if (single && is.null(x[["type"]])) return(location(path))
    x <- lapply(x, readable)
    if (single && !is.null(x[["uri"]])) x$uri <- shorthand_uri(path)
    x
  }
  readable(file)
}

# An adapter's own `uri`, the shorthand the engine expands into its file and
# index locations as JBrowse's config does: a local path as a file URI with no
# version query, since the engine names the index by appending to it
shorthand_uri <- function(path) {
  if (grepl("^[a-z]+://", path)) path
  else sub("\\?v=.*$", "", location(path)$uri)
}

# A table's columns as a data frame, an empty column read as an empty vector
# and a nested column kept a list column
as_frame <- function(columns) {
  columns <- lapply(columns, function(x)
    if (is.null(x) || (is.list(x) && !length(x))) logical(0) else x)
  n <- unique(lengths(columns))
  if (length(n) > 1) stop("the engine's columns differ in length", call. = FALSE)
  structure(columns, names = names(columns) %||% character(0),
    row.names = .set_row_names(if (length(n)) n else 0L), class = "data.frame")
}

`%||%` <- function(a, b) if (is.null(a)) b else a
