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
      return(invisible())
    }
    if (is.null(ctx)) {
      fresh <- V8::v8()
      fresh$eval(paste(
        "var jbrowseHost = { readRange: function(url, start, end) {",
        "return console.r.call('ggjbrowse:::read_range',",
        "{ url: url, start: start, end: end }) } }"))
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

call_engine <- function(table, args) {
  json <- engine()$call(sprintf("jbrowse.%s", table),
    as.character(jsonlite::toJSON(args, auto_unbox = TRUE)), await = TRUE)
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

# A table's columns as a data frame, an empty column read as an empty vector
as_frame <- function(columns) {
  as.data.frame(lapply(columns, function(x)
    if (is.list(x)) unlist(x) %||% logical(0) else x),
    stringsAsFactors = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
