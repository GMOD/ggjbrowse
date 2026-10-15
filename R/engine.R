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
      if (!identical(Sys.getpid(), deferred$pid) || forked()) {
        rlang::abort(paste("the engine cannot start in a forked worker",
          "(parallel::mclapply(), future's multicore, BiocParallel's",
          "MulticoreParam), where V8 hangs: read in a socket cluster, such",
          "as parallel::makeCluster() or future's multisession, whose",
          "workers load ggjbrowse afresh"), class = "ggjbrowse_forked")
      }
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

# The process that loaded ggjbrowse, the only one engine() starts V8 in: a
# forked child inherits a started engine, but one starting there hangs
.onLoad <- function(libname, pkgname) {
  deferred$pid <- Sys.getpid()
}

# Whether this process is a worker parallel's fork made, where a V8 the parent
# loaded hangs starting, even when ggjbrowse loaded only after the fork
forked <- function() get("isChild", asNamespace("parallel"))()

read_range_deferring <- function(url, start, end) {
  withCallingHandlers(read_range(url, start, end), warning = function(w) {
    deferred$warnings <- c(deferred$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}

# Runs use(<the engine's context>) and raises what the engine queued meanwhile.
# V8 raises a rejected call's error with JavaScript's own "Error: " prefix,
# which R would print after its own. `display` is the config the call passed,
# which an error about it names.
in_engine <- function(use, display = NULL) {
  ctx <- engine()
  on.exit({
    queued <- ctx$get("ggjbrowseQueued")
    ctx$eval("ggjbrowseQueued = { warnings: [], errors: [] }")
    warnings <- unique(c(deferred$warnings, unlist(queued$warnings)))
    deferred$warnings <- character(0)
    for (w in warnings) warning(w, call. = FALSE)
    if (length(queued$errors)) {
      engine_error(paste(unique(unlist(queued$errors)), collapse = "\n"),
        display)
    }
  })
  tryCatch(use(ctx), "C++Error" = function(e) {
    engine_error(sub("^Error: ", "", conditionMessage(e)), display)
  })
}

TOO_LARGE <- paste0("^(.*) (spans?|needs?) ([0-9]+) (bp|bytes) of (\\w+)",
  "( in total)?, over the ([0-9]+)-(bp|byte) limit; read a smaller region ",
  "or raise (spanLimit|byteLimit)$")

# The engine's error as R raises it: a region over a limit as
# `ggjbrowse_too_large`, naming R's argument and the loci as the location box
# writes them, and a display slot the config schema rejects as
# `ggjbrowse_bad_display`, naming the slot
engine_error <- function(message, display = NULL) {
  big <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE)
  m <- regmatches(message, regexec(TOO_LARGE, message))[[1]]
  if (length(m)) {
    argument <- c(spanLimit = "span_limit", byteLimit = "byte_limit")[[m[10]]]
    loci <- m[2]
    at <- regmatches(loci, gregexpr("[^ ,]+:[0-9]+-[0-9]+", loci))[[1]]
    if (length(at)) {
      parts <- regmatches(at, regexec("^(.+):([0-9]+)-([0-9]+)$", at))
      loci <- paste(vapply(parts, function(p) sprintf("%s:%s-%s", p[2],
        big(as.numeric(p[3]) + 1), big(p[4])), ""), collapse = ", ")
    }
    rlang::abort(sprintf("%s %s %s %s of %s%s, past `%s` (%s): read a smaller region or raise `%s`",
      loci, m[3], big(m[4]), m[5], m[6], m[7], argument,
      paste(big(m[8]), if (m[9] == "bp") "bp" else "bytes"), argument),
      class = "ggjbrowse_too_large", limit = argument, call = NULL)
  }
  if (grepl("-row limit", message, fixed = TRUE)) {
    rlang::abort(message, class = "ggjbrowse_too_large", limit = "rows",
      call = NULL)
  }
  union <- regmatches(message, regexec("No matching type for union (\\S+)",
    message))[[1]]
  if (length(union)) {
    slot <- Filter(function(s) grepl(tolower(s), tolower(union[2]),
      fixed = TRUE), names(display))
    value <- if (length(slot)) unlist(display[[slot[1]]])
    rlang::abort(if (length(slot)) {
      sprintf(paste("`display$%s` is %s, which the display's config schema",
        "does not accept (no type of its union %s matches)"), slot[1],
        paste(if (is.character(value)) encodeString(value, quote = "\"")
          else format(value), collapse = ", "), union[2])
    } else {
      sprintf(paste("the display config holds a value its schema does not",
        "accept (no type of its union %s matches)"), union[2])
    }, class = "ggjbrowse_bad_display", call = NULL)
  }
  rlang::abort(message, call = NULL)
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
  }, args$display)
  bytes <- rawConnection(wire$bytes)
  on.exit(close(bytes))
  x <- from_wire(jsonlite::fromJSON(wire$json, simplifyVector = TRUE,
    simplifyDataFrame = FALSE), bytes)
  for (w in unique(unlist(if (is.list(x)) x$warnings))) {
    warning(w, call. = FALSE)
  }
  x
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
  name <- tolower(bare_name(file))
  for (format in formats) {
    if (!grepl(format$pattern, name)) next
    adapter <- list(type = format$type)
    adapter[[format$slot]] <- location(file)
    if (!is.null(format$index)) {
      index <- index %||% beside(file, paste0(".", tolower(format$index)))
      adapter$index <- list(location = location(index),
        indexType = if (grepl("\\.csi$", tolower(bare_name(index)))) "CSI"
          else format$index)
    }
    return(adapter)
  }
  stop("cannot tell the format of ", file, ": expected ", expected,
    call. = FALSE)
}

# The file beside `file` its name and `suffix` name, a URL's query and
# fragment kept after it, as a signed URL's index carries the same query; a
# local path's `?` and `#` are part of its name
beside <- function(file, suffix) {
  if (!grepl("^[a-z]+://", file)) return(paste0(file, suffix))
  sub("^([^?#]*)", paste0("\\1", suffix), file)
}

# A path or URL as its extension reads, a URL's query and fragment aside
bare_name <- function(file) {
  if (grepl("^[a-z]+://", file)) sub("[?#].*$", "", file) else file
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
