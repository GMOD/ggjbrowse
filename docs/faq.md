# FAQ

## How does it work?

- **In V8:** JBrowse's `@jbrowse/tables` package. It takes a region through
  the same fetching and layout as a display and returns the result as
  tables. JBrowse's canvas drawing and UI do not run.
- **In R:** the readers pass file bytes to that package, and ggplot2 draws the
  tables it returns.
- `engine()` in `R/engine.R` starts one [V8](https://CRAN.R-project.org/package=V8)
  context per R session and sources the bundle into it. The context stays
  alive between calls, so the engine's caches do too.
- `jbrowse_reset()` drops the context, which clears the caches.

## What happens in one call?

One `jbrowse_alignments(bam, "chr12:6,534,201-6,538,700")` call, from path to plot:

```mermaid
sequenceDiagram
    participant R as R reader
    participant V as V8 context
    participant H as read_range()
    participant F as File or URL

    R->>R: regions_of() parses the locus
    R->>R: adapter_of() and location() build the adapter config
    R->>V: call_engine("alignments", JSON)
    loop each byte range the engine needs
        V->>H: jbrowseHost.readRange(url, start, end)
        H->>F: seek and readBin, or curl range request
        F-->>H: bytes
        H-->>V: bytes, base64 through jsonlite
    end
    V-->>R: columns as JSON
    R->>R: as_frame() builds the data frames
    R->>R: in_engine() raises queued warnings and errors
    R->>R: autoplot() hands each frame to its geom
```

## Who decides what?

- The engine decides whatever JBrowse decides from the track's config, such as
  colors and the display's geometry.
- R decides whatever depends on the figure, such as pixels per base, inches,
  label sizes and which ggplot layer draws what.
- [agent-docs/engine-owns-display-decisions.md](../agent-docs/engine-owns-display-decisions.md)
  lists the rules R still copies.

## Why V8 and not QuickJS?

- JBrowse's engine is a 3 MB bundle built for a browser-like host. V8 runs it
  at full speed with full ES support, and the `V8` package is on CRAN.
- QuickJS could host the same bundle through the same `readRange` callback.
  Nothing in the engine depends on V8.

## Is the engine a reimplementation of JBrowse?

- No. `inst/js/jbrowse-tables.js` is `@jbrowse/tables`, built from
  jbrowse-components' own source and run with no web globals. JBrowse's tests
  hold it to the browser's output.

## How does it compare with Rsamtools?

- Rsamtools cannot open a CRAM. `jbrowse_alignments()` reads one when given
  its reference as `assembly`.
- Rsamtools parses faster. On a local 50 kb region with 9,596 reads,
  `scanBam()` took 0.02 s and `jbrowse_alignments()` 0.37 s; on a 73-read
  region, 0.02 s and 0.07 s (median of 5, ggjbrowse 0.7.0, Rsamtools 2.26.0).
- The two do different work. `scanBam()` returns raw fields, and
  `jbrowse_alignments()` also packs pileup rows and computes mismatches,
  coverage, junctions and colors.
- The first call in a session adds about 0.2 s to start the engine.

## Why is the bundle committed?

- Installing the package needs no Node or pnpm.
- The first line of the file names the release, or the jbrowse-components
  commit, it came from.
- [`inst/COPYRIGHTS`](../inst/COPYRIGHTS) names the license and copyright of
  each package bundled in it.

## How does it read a file in R?

- The JavaScript cannot touch a disk or a network. It calls
  `jbrowseHost.readRange(url, start, end)`, which V8 routes back to R's
  `read_range()`.
- `read_range()` seeks and reads a local file, or sends a `curl` range request
  for a remote one. The JavaScript waits on the answer, so the call is
  synchronous.
- jsonlite sends the bytes to the JavaScript base64-encoded.

## Why does a changed remote file still show old data?

- The engine caches every read by URL. A local file's modification time and
  size are part of its URL, so an overwritten local file rereads.
- A remote file has no such key. `jbrowse_reset()` clears the caches.

## Can it read a file without an index?

- A BAM or CRAM needs its index, and a bgzipped GFF3, GTF, BED or VCF needs a
  tabix or CSI index. Put it beside the file under its usual name, or pass
  `index =`.
- A plain GFF3, GTF, BED or VCF has no index. The engine reads the whole file.

## A CRAM gives an error about the reference

- A CRAM stores reads against the FASTA it was encoded with. Pass that FASTA
  or 2bit as `assembly =`.

## The files spell a contig differently

- Pass each reader `assembly = jbrowse_assembly(fasta, aliases)`, the FASTA
  and an alias table listing each contig's names on a row, such as UCSC's
  `chromAlias.txt`. A region may name the contig as any file does. The
  engine reads each file under that file's spelling, and the reader names the
  contig in its tables as the region does.
- A JBrowse config's assembly works as is: `assembly = list(name = "hg19",
  uri = "hg19.fa.gz", refNameAliases = list(uri = "hg19_aliases.txt"))`.

## Do R warnings from inside the engine work with `tryCatch`?

- Yes, after the call. The engine's `console.warn` and `console.error` go to
  a queue, and `in_engine()` raises it once V8 returns. A handler that unwinds
  through V8's frames would break the context.

## Which coordinates do I give?

- A locus string is 1-based and inclusive, as in JBrowse's location box.
- `region = "chr12", start, end` is 0-based half-open.

## Can I restyle a figure?

- Yes. `autoplot()` returns an ordinary ggplot, so `+` adds layers, scales and
  themes. Each table is also a data frame to plot yourself with
  `geom_jbrowse_*()`.
- `autoplot(a, "pileup")` draws one band, and `jbrowse_stack()` stacks bands
  and tracks on one axis.

## How do I update the engine?

- Run `pnpm install` in a jbrowse-components checkout, then
  `scripts/build-js.sh` here. The script looks in `../jbrowse-components`, or
  in `$JBROWSE_COMPONENTS`.
- [CONTRIBUTING.md](../CONTRIBUTING.md#updating-the-engine) covers the CI that
  takes each new release.
