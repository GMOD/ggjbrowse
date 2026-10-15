# How does this work?

- **In:** a BAM, CRAM, feature, wiggle, VCF or PAF file, and a region.
- **Out:** the data frame JBrowse 2's display draws from, then a ggplot2 layer.

## The engine

- `inst/js/jbrowse-tables.js` is `@jbrowse/tables` from
  [jbrowse-components](https://github.com/GMOD/jbrowse-components), bundled
  to run without a DOM or a web worker. Its first line names the commit.
- `engine()` in `R/engine.R` starts one [V8](https://CRAN.R-project.org/package=V8)
  context per R session and sources the bundle into it. The context stays
  alive between calls, so the engine's caches do too.
- `jbrowse_reset()` drops the context, which clears the caches.
- The engine runs on V8 through the `V8` package, which handles the 3 MB bundle
  at full speed with full ES support. QuickJS could host it the same way.
- The bundle is committed, so installing needs no JavaScript toolchain.
- JBrowse's own tests hold the bundle, run with no web globals, to the
  browser's engine.
- [`inst/COPYRIGHTS`](../inst/COPYRIGHTS) names the bundle's source and the
  license and copyright of each package bundled in it.
- [CONTRIBUTING.md](../CONTRIBUTING.md#updating-the-engine) covers rebuilding
  it and the CI that takes each new release.

## A call

- `jbrowse_alignments()` and the other `jbrowse_*()` functions turn a path or
  URL into an adapter config (`guess_adapter()`) and a region into a locus.
- `call_engine()` sends them to `jbrowse.<table>` as JSON and waits on the
  promise (`await = TRUE`).
- The engine returns columns as JSON. `as_frame()` makes each table a data
  frame.
- `autoplot()` and the geoms in `R/geoms-alignments.R` draw the frames.

## Reading files

- The JavaScript cannot touch a disk or a network. It calls
  `jbrowseHost.readRange(url, start, end)`, which V8 routes back to R's
  `read_range()`.
- `read_range()` seeks and reads a local file, or sends a `curl` range request
  for a remote one. The JavaScript waits on the answer, so the call is
  synchronous.
- jsonlite sends the bytes to the JavaScript base64-encoded.
- `location()` puts a local file's modification time and size in its URL, so
  the engine rereads a file overwritten in place. A remote file changed on its server
  needs `jbrowse_reset()`.

## Warnings and errors

- The engine's `console.warn` and `console.error` go to a queue, and
  nothing raises yet.
- `in_engine()` raises the queue once the V8 call has returned. An R handler
  that unwinds through V8's frames breaks the context, so nothing may raise
  inside the call.

## Who decides what

- The engine decides whatever JBrowse decides from the track's config, such as
  colors and the display's geometry.
- R decides whatever depends on the figure, such as pixels per base, inches,
  label sizes and which ggplot layer draws what.
- [`agent-docs/engine-owns-display-decisions.md`](../agent-docs/engine-owns-display-decisions.md)
  lists the rules R still copies.

## Updating the bundle

- Run `pnpm install` in a jbrowse-components checkout, then
  `scripts/build-js.sh` here.
- The script looks in `../jbrowse-components`, or in `$JBROWSE_COMPONENTS`.
