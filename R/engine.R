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

call_engine <- function(table, args) {
  json <- in_engine(function(ctx) ctx$call(sprintf("jbrowse.%s", table),
    as.character(jsonlite::toJSON(args, auto_unbox = TRUE)), await = TRUE))
  jsonlite::fromJSON(json, simplifyVector = TRUE, simplifyDataFrame = FALSE)
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

# An adapter config for `file`, from the first of `formats` whose `pattern`
# matches its name without any URL query; the error `expected` when none does.
# A format names its adapter `type`, the `slot` holding the file's location,
# an optional `index` (its `suffix` beside the file, its default `kind`, and a
# `slot` to hold it in place of an `index` block) and `extra` slots, a function
# of the file.
guess_adapter <- function(file, index, formats, expected = NULL) {
  name <- tolower(sub("\\?.*$", "", file))
  for (format in formats) {
    if (!grepl(format$pattern, name)) next
    adapter <- list(type = format$type)
    adapter[[format$slot]] <- location(file)
    if (!is.null(format$extra)) adapter <- c(adapter, format$extra(file))
    if (!is.null(format$index)) {
      index <- index %||% paste0(file, format$index$suffix)
      if (!is.null(format$index$slot)) {
        adapter[[format$index$slot]] <- location(index)
      } else {
        adapter$index <- list(location = location(index),
          indexType = if (grepl("\\.csi$", index)) "CSI" else format$index$kind)
      }
    }
    return(adapter)
  }
  stop("cannot tell the format of ", file, ": expected ", expected, call. = FALSE)
}

# The adapter config the engine reads `file` through: a track's own adapter
# config as given, each of its locations a path or URL made readable, or else
# the one `guess` makes of a path or URL
adapter_of <- function(file, guess) {
  if (!is.list(file)) return(guess())
  readable <- function(x) {
    if (!is.list(x)) return(x)
    path <- x[["uri"]] %||% x[["localPath"]]
    if (is.character(path) && length(path) == 1) location(path)
    else lapply(x, readable)
  }
  readable(file)
}

# A table's columns as a data frame, an empty column read as an empty vector
as_frame <- function(columns) {
  as.data.frame(lapply(columns, function(x)
    if (is.list(x)) unlist(x) %||% logical(0) else x),
    stringsAsFactors = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
