#!/usr/bin/env Rscript
# Builds docs/tutorials/<name>.md from tutorials/<name>.Rmd.
#
# A source is R Markdown whose `<!-- jbrowse ... -->` comments pull sections
# out of the matching JBrowse website tutorial, at the jbrowse-components
# commit tutorials/JBROWSE_COMPONENTS pins. The script expands them, knits the
# result and writes the markdown and its figures.
#
#   Rscript scripts/tutorials.R [name ...]   expand and knit
#   Rscript scripts/tutorials.R --expand     write the expanded .Rmd to a
#                                            temporary file and print its path
#   Rscript scripts/tutorials.R --check      knit to a temporary directory and
#                                            fail if the markdown differs
#   Rscript scripts/tutorials.R --status     list website commits since the pin
#
# The website pages come from the checkout $JBROWSE_COMPONENTS names (default
# ~/src/jbrowse-components, as scripts/gallery.R has it), or from GitHub when
# that checkout lacks the pin. --status needs the checkout.
#
# A directive is YAML inside the comment:
#
#   <!-- jbrowse
#   section: Producing the BAM     # a heading's text, or `intro`
#   title: Aligning the reads      # renames the heading; `false` drops it
#   drop: [bwameth.py index]       # drops each block holding the text
#   keep: [MethylDackel extract]   # keeps a block the defaults would drop
#   replace:                       # each `from` must occur exactly once
#     - {from: JBrowse makes, to: The engine makes}
#   -->
#
# A directive with `assign` writes an R chunk instead, each variable taking a
# value out of the track config whose trackId starts its path, so the R code
# draws with the settings the website's track loads. A number in the path
# steps into an array, from 1 (`multi.adapter.subadapters.2.uri`):
#
#   <!-- jbrowse
#   chunk: track                   # the chunk's label
#   assign:
#     gff: dtu_muscle_vs_liver.uri
#     display: dtu_muscle_vs_liver.displayDefaults
#   -->
#
# Every `drop`, `keep` and `replace` must match, so a website edit that removes
# the text fails the build instead of passing silently. By default a pulled
# section loses what only the browser uses: JSON configs, `jbrowse` CLI
# commands, `<Figure>` and `<Video>` tags, `<details>` holding a config, and
# the paragraph ending in a colon that introduced one of them.

SITE <- "https://jbrowse.org/jb2"
REPO <- "GMOD/jbrowse-components"
SRC_DIR <- "tutorials"
OUT_DIR <- file.path("docs", "tutorials")

pin <- function() trimws(readLines(file.path(SRC_DIR, "JBROWSE_COMPONENTS"), 1))

checkout <- function() {
  Sys.getenv("JBROWSE_COMPONENTS", path.expand("~/src/jbrowse-components"))
}

website_cache <- new.env()

website_file <- function(path) {
  key <- path
  if (!is.null(website_cache[[key]])) return(website_cache[[key]])
  sha <- pin()
  repo <- checkout()
  text <- NULL
  if (dir.exists(repo)) {
    text <- suppressWarnings(tryCatch(
      system2("git", c("-C", repo, "show", paste0(sha, ":", path)),
        stdout = TRUE, stderr = FALSE),
      error = function(e) NULL))
    if (!is.null(attr(text, "status"))) text <- NULL
  }
  if (is.null(text)) {
    url <- sprintf("https://raw.githubusercontent.com/%s/%s/%s", REPO, sha, path)
    text <- tryCatch(readLines(url, warn = FALSE, encoding = "UTF-8"),
      error = function(e) NULL)
  }
  website_cache[[key]] <- text
  text
}

tutorial_lines <- function(name) {
  lines <- website_file(sprintf("website/docs/tutorials/%s.md", name))
  if (is.null(lines)) stop("no website tutorial named ", name, call. = FALSE)
  lines
}

split_frontmatter <- function(lines) {
  if (length(lines) && lines[1] == "---") {
    end <- which(lines == "---")[2]
    list(meta = yaml::yaml.load(paste(lines[2:(end - 1)], collapse = "\n")),
      body = lines[-seq_len(end)])
  } else {
    list(meta = list(), body = lines)
  }
}

# A block is a fence, a component tag, a <details>, an aside, a comment line or
# a run of non-blank lines.
parse_blocks <- function(lines) {
  blocks <- list()
  i <- 1
  n <- length(lines)
  add <- function(kind, from, to, info = "") {
    blocks[[length(blocks) + 1]] <<- list(kind = kind, info = info,
      lines = lines[from:to])
  }
  while (i <= n) {
    line <- lines[i]
    if (!nzchar(trimws(line))) {
      i <- i + 1
      next
    }
    j <- i
    if (grepl("^\\s*(```|~~~)", line)) {
      marker <- sub("^\\s*((```|~~~)`*).*", "\\1", line)
      repeat {
        j <- j + 1
        if (j > n || startsWith(trimws(lines[j]), marker)) break
      }
      add("fence", i, min(j, n), trimws(sub("^\\s*(```|~~~)`*", "", line)))
    } else if (grepl("^#{2,6} ", line)) {
      add("heading", i, i)
    } else if (grepl("^<(Figure|Video)\\b", line)) {
      while (j < n && !grepl("/>\\s*$", lines[j])) j <- j + 1
      add("media", i, j)
    } else if (grepl("^<details", line)) {
      while (j < n && !grepl("</details>", lines[j])) j <- j + 1
      add("details", i, j)
    } else if (grepl("^:::[a-z]", line)) {
      while (j < n && !grepl("^:::\\s*$", lines[j + 1])) j <- j + 1
      add("aside", i, min(j + 1, n))
      j <- j + 1
    } else if (grepl("^<!--.*-->\\s*$", line)) {
      add("comment", i, i)
    } else {
      while (j < n && nzchar(trimws(lines[j + 1])) &&
        !grepl("^\\s*(```|~~~)|^#{2,6} |^<(Figure|Video|details)\\b|^:::",
          lines[j + 1])) {
        j <- j + 1
      }
      add("text", i, j)
    }
    i <- j + 1
  }
  blocks
}

heading_level <- function(b) nchar(sub(" .*", "", b$lines[1]))

heading_text <- function(b) {
  trimws(sub("\\s*\\{#[^}]*\\}\\s*$", "", sub("^#+ ", "", b$lines[1])))
}

heading_id <- function(b) {
  explicit <- regmatches(b$lines[1], regexpr("\\{#[^}]*\\}", b$lines[1]))
  if (length(explicit)) gsub("^\\{#|\\}$", "", explicit) else slug(heading_text(b))
}

# GitHub's and Starlight's slugger: lowercase, punctuation out, spaces to dashes
slug <- function(text) {
  text <- tolower(text)
  text <- gsub("[^[:alnum:] _-]", "", text, perl = TRUE)
  gsub(" ", "-", text)
}

# The blocks of a heading's section, down to the next heading at its level or
# above. `intro` is everything before the first heading.
section_blocks <- function(blocks, section) {
  heads <- which(vapply(blocks, function(b) b$kind == "heading", TRUE))
  if (identical(section, "intro")) {
    end <- if (length(heads)) heads[1] - 1 else length(blocks)
    return(list(head = NULL, body = blocks[seq_len(end)]))
  }
  hit <- heads[vapply(blocks[heads], heading_text, "") == section]
  if (length(hit) != 1) {
    stop(sprintf("the website tutorial has %d sections headed '%s'",
      length(hit), section), call. = FALSE)
  }
  level <- heading_level(blocks[[hit]])
  later <- heads[heads > hit]
  later <- later[vapply(blocks[later], heading_level, 0) <= level]
  end <- if (length(later)) later[1] - 1 else length(blocks)
  list(head = blocks[[hit]], body = blocks[seq_len(end - hit) + hit])
}

browser_only <- function(b) {
  switch(b$kind,
    fence = grepl("^json", b$info) || (b$info %in% c("bash", "sh") &&
      any(grepl("^(jbrowse|jb|npx)\\b",
        trimws(b$lines[-c(1, length(b$lines))])))),
    media = TRUE,
    comment = TRUE,
    details = any(grepl("^\\s*```json", b$lines)),
    FALSE)
}

block_text <- function(b) paste(b$lines, collapse = "\n")

matches <- function(blocks, text) {
  which(vapply(blocks, function(b) grepl(text, block_text(b), fixed = TRUE),
    TRUE))
}

one_match <- function(blocks, text, what, section) {
  hit <- matches(blocks, text)
  if (length(hit) != 1) {
    stop(sprintf("%s '%s' matches %d blocks of section '%s'", what, text,
      length(hit), section), call. = FALSE)
  }
  hit
}

pull <- function(spec, tutorial, level_offset) {
  blocks <- parse_blocks(split_frontmatter(tutorial_lines(tutorial))$body)
  section <- spec[["section"]]
  s <- section_blocks(blocks, section)
  body <- s$body
  dropped <- vapply(body, browser_only, TRUE)
  for (k in seq_along(body)) {
    nxt <- k + 1
    while (nxt <= length(body) && body[[nxt]]$kind == "comment") nxt <- nxt + 1
    if (body[[k]]$kind == "text" && nxt <= length(body) && dropped[nxt] &&
      grepl(":\\s*$", block_text(body[[k]]))) {
      dropped[k] <- TRUE
    }
  }
  if (length(spec[["only"]])) {
    dropped[] <- TRUE
    for (text in spec[["only"]]) dropped[one_match(body, text, "only", section)] <- FALSE
  }
  for (text in spec[["keep"]]) dropped[one_match(body, text, "keep", section)] <- FALSE
  for (text in spec[["drop"]]) dropped[one_match(body, text, "drop", section)] <- TRUE
  body <- body[!dropped]
  for (text in spec[["drop_lines"]]) {
    k <- one_match(body, text, "drop_lines", section)
    hit <- grepl(text, body[[k]]$lines, fixed = TRUE)
    if (sum(hit) != 1) {
      stop(sprintf("drop_lines '%s' matches %d lines", text, sum(hit)),
        call. = FALSE)
    }
    body[[k]]$lines <- body[[k]]$lines[!hit]
  }
  for (r in spec[["replace"]]) {
    count <- vapply(body, function(b) lengths(regmatches(block_text(b),
      gregexpr(r$from, block_text(b), fixed = TRUE))), 0)
    if (sum(count) != 1) {
      stop(sprintf("replace '%s' matches %d times in section '%s'", r$from,
        sum(count), section), call. = FALSE)
    }
    k <- which(count == 1)
    body[[k]]$lines <- strsplit(sub(r$from, r$to, block_text(body[[k]]),
      fixed = TRUE), "\n")[[1]]
  }
  out <- vapply(body, function(b) {
    if (b$kind == "heading") {
      level <- heading_level(b) + level_offset
      return(paste(strrep("#", level), heading_text(b)))
    }
    if (b$kind == "aside") return(aside(b))
    if (b$kind %in% c("text", "details")) return(links(block_text(b), tutorial))
    block_text(b)
  }, "")
  text <- paste(out, collapse = "\n\n")
  head <- NULL
  if (!is.null(s$head) && !isFALSE(spec[["title"]])) {
    level <- heading_level(s$head) + level_offset
    title <- if (is.character(spec[["title"]])) spec[["title"]] else heading_text(s$head)
    head <- paste(strrep("#", level), title)
  }
  paste(c(head, text), collapse = "\n\n")
}

# A Starlight aside becomes a GitHub alert.
aside <- function(b) {
  type <- sub("^:::([a-z]+).*", "\\1", b$lines[1])
  alert <- c(note = "NOTE", tip = "TIP", caution = "CAUTION",
    danger = "WARNING")[[type]]
  inner <- b$lines[-c(1, length(b$lines))]
  paste(c(sprintf("> [!%s]", alert), paste0(">", ifelse(nzchar(inner),
    paste0(" ", inner), ""))), collapse = "\n")
}

page_title <- function(path) {
  path <- sub("/$", "", sub("#.*", "", path))
  for (file in sprintf("website%s%s", path, c(".md", "/index.md", ".mdx"))) {
    lines <- website_file(file)
    if (!is.null(lines)) return(split_frontmatter(lines)$meta$title)
  }
  basename(path)
}

# Site-relative links go to jbrowse.org; an anchor stays local when the
# expanded page has that heading, which resolve_anchors() decides at the end.
links <- function(text, tutorial) {
  empty <- gregexpr("\\[\\]\\((/docs/[^)\\s]*)\\)", text, perl = TRUE)
  regmatches(text, empty) <- lapply(regmatches(text, empty), function(m)
    vapply(m, function(link) {
      path <- sub("^\\[\\]\\(([^)]*)\\)$", "\\1", link)
      sprintf("[%s](%s)", page_title(path), path)
    }, ""))
  text <- gsub("\\]\\((/[^)\\s]*)\\)", sprintf("](%s\\1)", SITE), text,
    perl = TRUE)
  gsub("\\]\\(#([^)\\s]*)\\)", sprintf("](#@@%s@@\\1)", tutorial), text,
    perl = TRUE)
}

resolve_anchors <- function(text) {
  lines <- strsplit(text, "\n")[[1]]
  fence <- cumsum(grepl("^\\s*```", lines)) %% 2 == 1
  heads <- lines[!fence & grepl("^#{1,6} ", lines)]
  ids <- slug(trimws(sub("^#+ ", "", heads)))
  pattern <- "\\(#@@([a-z0-9_]+)@@([^)]*)\\)"
  m <- gregexpr(pattern, text, perl = TRUE)
  regmatches(text, m) <- lapply(regmatches(text, m), function(found)
    vapply(found, function(f) {
      tutorial <- sub(pattern, "\\1", f, perl = TRUE)
      anchor <- sub(pattern, "\\2", f, perl = TRUE)
      if (anchor %in% ids) sprintf("(#%s)", anchor)
      else sprintf("(%s/docs/tutorials/%s/#%s)", SITE, tutorial, anchor)
    }, ""))
  text
}

# The website tutorial's JSON fences that configure a track, a list of them
# per trackId, since a page may show one track's config more than once
track_configs <- function(tutorial) {
  blocks <- parse_blocks(split_frontmatter(tutorial_lines(tutorial))$body)
  configs <- list()
  for (b in blocks) {
    if (b$kind != "fence" || !grepl("^json", b$info)) next
    config <- tryCatch(jsonlite::fromJSON(paste(b$lines[-c(1,
      length(b$lines))], collapse = "\n"), simplifyDataFrame = FALSE,
      simplifyMatrix = FALSE), error = function(e) NULL)
    id <- config$trackId
    if (is.null(id)) next
    configs[[id]] <- c(configs[[id]], list(config))
  }
  configs
}

# The value at `path` in a track config: a name steps into an object, a number
# into an array, counting from 1 as R does
config_value <- function(config, path) {
  for (step in path) {
    if (is.null(config)) return(NULL)
    config <- if (grepl("^[0-9]+$", step)) {
      i <- as.integer(step)
      if (i <= length(config)) config[[i]]
    } else config[[step]]
  }
  config
}

# An R expression that evaluates to `x`, indented to sit `indent` spaces in
# and wrapped one entry a line when it overruns the line
r_literal <- function(x, indent = 0) {
  key <- function(k) if (make.names(k) == k) k else sprintf("`%s`", k)
  scalar <- function(v) {
    if (is.character(v)) encodeString(v, quote = "\"")
    else if (is.logical(v)) as.character(v)
    else format(v, digits = 15)
  }
  if (is.null(x)) return("NULL")
  if (is.atomic(x) && is.null(names(x))) {
    items <- vapply(x, scalar, "")
    if (length(items) == 1) return(items)
  } else {
    if (is.atomic(x)) x <- as.list(x)
    items <- vapply(seq_along(x), function(i) {
      value <- r_literal(x[[i]], indent + 2)
      if (is.null(names(x)) || !nzchar(names(x)[i])) value
      else paste(key(names(x)[i]), "=", value)
    }, "")
  }
  open <- if (is.list(x)) "list(" else "c("
  flat <- paste0(open, paste(items, collapse = ", "), ")")
  if (indent + nchar(flat) <= 78 && !any(grepl("\n", items))) return(flat)
  pad <- strrep(" ", indent + 2)
  paste0(open, "\n", pad, paste(items, collapse = paste0(",\n", pad)), "\n",
    strrep(" ", indent), ")")
}

# A chunk assigning R variables from the website's track configs, each value
# named by its trackId and a path into the config: `gff: dtu.uri`
assign_chunk <- function(spec, tutorial) {
  configs <- track_configs(tutorial)
  lines <- vapply(names(spec[["assign"]]), function(name) {
    path <- strsplit(spec[["assign"]][[name]], ".", fixed = TRUE)[[1]]
    fences <- configs[[path[1]]]
    if (is.null(fences)) {
      stop(sprintf("no track config with trackId %s in the website tutorial",
        path[1]), call. = FALSE)
    }
    values <- Filter(Negate(is.null), lapply(fences, config_value, path[-1]))
    if (!length(values)) {
      stop(sprintf("track %s has no %s", path[1], spec[["assign"]][[name]]),
        call. = FALSE)
    }
    if (length(unique(values)) > 1) {
      stop(sprintf("the website's configs of track %s disagree on %s",
        path[1], spec[["assign"]][[name]]), call. = FALSE)
    }
    paste(name, "<-", r_literal(values[[1]], 0))
  }, "")
  header <- paste0("```{r", if (!is.null(spec[["chunk"]]))
    paste0(" ", spec[["chunk"]]), "}")
  paste(c(header, lines, "```"), collapse = "\n")
}

expand <- function(name) {
  source_file <- file.path(SRC_DIR, paste0(name, ".Rmd"))
  doc <- split_frontmatter(readLines(source_file, encoding = "UTF-8"))
  tutorial <- doc$meta$jbrowse_tutorial %||% name
  lines <- doc$body
  out <- character()
  i <- 1
  while (i <= length(lines)) {
    if (grepl("^<!-- jbrowse\\b", lines[i])) {
      j <- i
      while (!grepl("-->\\s*$", lines[j])) j <- j + 1
      yaml_text <- sub("^<!-- jbrowse", "", paste(lines[i:j], collapse = "\n"))
      spec <- yaml::yaml.load(sub("-->\\s*$", "", yaml_text))
      if (is.character(spec)) spec <- list(section = spec)
      from <- spec[["tutorial"]] %||% tutorial
      out <- c(out, if (!is.null(spec[["assign"]])) assign_chunk(spec, from)
        else pull(spec, from, spec[["level"]] %||% 0))
      i <- j + 1
    } else {
      out <- c(out, lines[i])
      i <- i + 1
    }
  }
  body <- resolve_anchors(paste(out, collapse = "\n"))
  if (grepl("`r ", gsub("```[^`]*```", "", body))) {
    stop(name, ": pulled text holds `r , which knitr would run", call. = FALSE)
  }
  header <- c(
    sprintf("<!-- Generated by scripts/tutorials.R from %s and %s's website/docs/tutorials/%s.md at %s. Edit those, not this file. -->",
      source_file, REPO, tutorial, substr(pin(), 1, 10)),
    "",
    paste("#", doc$meta$title))
  paste(c(header, "", sub("^\\s+", "", body)), collapse = "\n")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

knit <- function(name, rmd, out_dir) {
  dir.create(file.path(out_dir, "figures"), recursive = TRUE,
    showWarnings = FALSE)
  input <- file.path(tempdir(), paste0(name, ".Rmd"))
  writeLines(rmd, input, useBytes = TRUE)
  output <- normalizePath(file.path(out_dir, paste0(name, ".md")),
    mustWork = FALSE)
  knitr::opts_knit$set(base.dir = normalizePath(out_dir))
  knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
    fig.path = sprintf("figures/%s-", name), dev = "png", dpi = 110,
    fig.width = 10, fig.height = 5, message = FALSE, error = FALSE)
  env <- new.env(parent = globalenv())
  knitr::knit(input, output, envir = env, quiet = TRUE, encoding = "UTF-8")
  output
}

names_of <- function(args) {
  args <- args[!startsWith(args, "--")]
  if (length(args)) return(args)
  sub("\\.Rmd$", "", list.files(SRC_DIR, "\\.Rmd$"))
}

main <- function(args) {
  tutorials <- names_of(args)
  if ("--status" %in% args) {
    for (name in tutorials) {
      meta <- split_frontmatter(readLines(file.path(SRC_DIR,
        paste0(name, ".Rmd"))))$meta
      page <- sprintf("website/docs/tutorials/%s.md",
        meta$jbrowse_tutorial %||% name)
      log <- system2("git", c("-C", checkout(), "log", "--oneline",
        paste0(pin(), "..origin/main"), "--", page), stdout = TRUE)
      cat(sprintf("%s: %d website commits since the pin\n", name, length(log)))
      if (length(log)) cat(paste0("  ", log, "\n"), sep = "")
    }
    return(invisible())
  }
  stale <- character()
  for (name in tutorials) {
    message(name)
    rmd <- expand(name)
    if ("--expand" %in% args) {
      expanded <- file.path(tempdir(), paste0(name, ".Rmd"))
      writeLines(rmd, expanded, useBytes = TRUE)
      message("  ", expanded)
      next
    }
    if ("--check" %in% args) {
      scratch <- file.path(tempdir(), "tutorials")
      fresh <- readLines(knit(name, rmd, scratch), encoding = "UTF-8")
      committed <- file.path(OUT_DIR, paste0(name, ".md"))
      if (!file.exists(committed) ||
        !identical(fresh, readLines(committed, encoding = "UTF-8"))) {
        stale <- c(stale, name)
      }
    } else {
      knit(name, rmd, OUT_DIR)
    }
  }
  if (length(stale)) {
    stop("out of date, rerun scripts/tutorials.R: ",
      paste(stale, collapse = ", "), call. = FALSE)
  }
}

main(commandArgs(trailingOnly = TRUE))
