# The engine is @jbrowse/tables' bare-host bundle, run in one V8 context per
# session. Every byte it reads comes through read_range(), which JavaScript
# calls synchronously; jsonlite sends the raw body back base64-encoded.

read_range <- function(url, start, end) {
  path <- sub("^file://", "", url)
  if (!grepl("^[a-z]+://", path)) {
    path <- utils::URLdecode(path)
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
  function() {
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

call_engine <- function(table, args) {
  json <- engine()$eval(
    sprintf("jbrowse.%s(%s)", table,
      jsonlite::toJSON(jsonlite::toJSON(args, auto_unbox = TRUE))),
    await = TRUE)
  jsonlite::fromJSON(json, simplifyVector = TRUE, simplifyDataFrame = FALSE)
}

location <- function(path) {
  uri <- if (grepl("^[a-z]+://", path)) path else
    paste0("file://", normalizePath(path, mustWork = TRUE))
  list(uri = uri, locationType = "UriLocation")
}

# A table's columns as a data frame, an empty column read as an empty vector
as_frame <- function(columns) {
  as.data.frame(lapply(columns, function(x)
    if (is.list(x)) unlist(x) %||% logical(0) else x),
    stringsAsFactors = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
