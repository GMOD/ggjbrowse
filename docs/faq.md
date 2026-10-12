# FAQ

## Why V8 and not QuickJS?

- JBrowse's engine is a 3 MB bundle built for a browser-like host. V8 runs it
  at full speed with full ES support, and the `V8` package is on CRAN.
- QuickJS could host the same bundle through the same `readRange` callback.
  Nothing in the engine depends on V8.

## Is the engine a reimplementation of JBrowse?

- No. `inst/js/jbrowse-tables.js` is `@jbrowse/tables`, built from
  jbrowse-components' own source. JBrowse's tests hold it to the browser's
  output.
- R draws the tables. Where the browser decides a color or a geometry, the
  engine returns it.

## How does it compare with Rsamtools?

- Rsamtools cannot open a CRAM. `jbrowse_alignments()` reads one when given
  its `reference`.
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

## How does it read a file in R?

- The JavaScript calls `jbrowseHost.readRange(url, start, end)` and R's
  `read_range()` answers, from disk or over HTTP.
- See [how-does-this-work.md](how-does-this-work.md#reading-files).

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
  or 2bit as `reference =`.

## The BAM and the FASTA spell a contig differently

- If the BAM names a contig `1` and the FASTA `chr1`, pass `reference_chrom` with
  the FASTA's spelling.

## Do R warnings from inside the engine work with `tryCatch`?

- Yes, after the call. The engine queues its warnings and errors, and
  `in_engine()` raises them once V8 returns. A handler that unwinds through
  V8's frames would break the context.

## Which coordinates do I give?

- A locus string is 1-based and inclusive, as in JBrowse's location box.
- `region = "ctgA", start, end` is 0-based half-open.

## Can I restyle a figure?

- Yes. `autoplot()` returns an ordinary ggplot, so `+` adds layers, scales and
  themes. Each table is also a data frame to plot yourself with
  `geom_jbrowse_*()`.
- `autoplot(a, "pileup")` draws one band, and `jbrowse_stack()` stacks bands
  and tracks on one axis.

## How do I update the engine?

- See [CONTRIBUTING.md](../CONTRIBUTING.md#updating-the-engine).
