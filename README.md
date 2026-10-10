# ggjbrowse

Genome browser tracks as data frames for ggplot2. `jbrowse_alignments()` reads
a region of a BAM or CRAM through the [JBrowse 2](https://jbrowse.org) genome
browser's own data engine, run in V8, and returns what its alignments display
draws: each read on its pileup row, with its mismatches, deletions and splice
gaps, and the coverage. A plot is an ordinary ggplot: themes, scales,
`patchwork` and `ggsave()` all work.

```r
remotes::install_github("GMOD/ggjbrowse")
```

```r
library(ggjbrowse)
library(patchwork)

a <- jbrowse_alignments("rnaseq.bam", "ctgA", 400, 1300)
autoplot(a, "coverage") / autoplot(a, color_by = "strand")
```

![RNA-seq reads over two exons of volvox ctgA, split at their introns](man/figures/spliced.png)

The tables are plain data frames, so a figure can be any ggplot over them:

```r
library(ggplot2)

cram <- jbrowse_alignments("reads.cram", "chr1", 1e6, 1.002e6,
  reference = "hg38.fa")
ggplot(cram$mismatches, aes(position, frequency, colour = base)) +
  geom_point() +
  scale_colour_manual(values = jbrowse_colors()$bases)
```

A file can be a path or a URL; a remote file is read in byte ranges through
`curl`. CRAM needs the reference, as a FASTA with its `.fai` or a 2bit. There
is no Bioconductor dependency: the BAM, CRAM and FASTA readers are JBrowse's.

## The engine

`inst/js/jbrowse-tables.js` is
[`@jbrowse/tables`](https://github.com/GMOD/jbrowse-components/tree/main/packages/tables)
bundled for a bare JavaScript host, committed so installing needs no
JavaScript toolchain. Its first line names the jbrowse-components release and
commit it was built from. JBrowse's own tests hold the bundle, run with no web
globals, to the browser's engine.

To rebuild it from a jbrowse-components checkout beside this repo:

```console
cd ../jbrowse-components && pnpm install && cd -
scripts/build-js.sh
```

`scripts/check-against-jbrowse.sh` rebuilds, installs and tests in one step,
against `$JBROWSE_COMPONENTS` or `../jbrowse-components`. The `Bundle` workflow
runs it nightly against jbrowse-components `main`, and takes each new
jbrowse-components release into the committed bundle once R CMD check passes.

## Not yet

Feature tracks (GFF3, BED, BigBed gene models), BigWig and VCF; the
alignments display's color schemes beyond strand; insertions, clips and
modifications on reads.

## See also

- [JBrowse 2](https://jbrowse.org)
- [JBrowseR](https://github.com/GMOD/JBrowseR) - JBrowse 2 as an R htmlwidget
- [ggbandage](https://github.com/GMOD/ggbandage) and
  [ggtubemap](https://github.com/GMOD/ggtubemap) - pangenome graph figures as
  ggplot2 layers
