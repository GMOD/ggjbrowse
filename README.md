# ggjbrowse


> [!CAUTION]
> Very early results from work in progress!

Genome browser tracks as data frames for ggplot2. `jbrowse_alignments()` reads
a region of a BAM or CRAM through the [JBrowse 2](https://jbrowse.org) genome
browser's own data engine, run in V8, and returns what its alignments display
draws: each read on its pileup row, filled as the display's `color` setting
paints it, with its mismatches, deletions and splice gaps, the coverage, and
the display's key. A plot is an ordinary ggplot: themes, scales, `patchwork`
and `ggsave()` all work.

```r
remotes::install_github("GMOD/ggjbrowse")
```

```r
library(ggjbrowse)
library(patchwork)

a <- jbrowse_alignments("rnaseq.bam", "ctgA", 400, 1300, color = "strand")
autoplot(a, "coverage") / autoplot(a)
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

`color` and `filter` are the display's own settings objects, so what JBrowse
can color or filter by, this can: `color = "mapq"`,
`color = list(field = "insertSize", domain = c(200, 600))`,
`filter = list(flagExclude = 1796, tagFilters = list(list(tag = "HP", value = "1")))`.
Each read comes back with the fill the display paints it and the bucket it
files it under, and `a$legend` is the display's key.

A file can be a path or a URL; a remote file is read in byte ranges through
`curl`. CRAM needs the reference, as a FASTA with its `.fai` or a 2bit. There
is no Bioconductor dependency: the BAM, CRAM and FASTA readers are JBrowse's.

## The engine

`inst/js/jbrowse-tables.js` is
[`@jbrowse/tables`](https://github.com/GMOD/jbrowse-components/tree/main/packages/tables)
bundled for a bare JavaScript host, committed so installing needs no
JavaScript toolchain. Its first line names the release, or the
jbrowse-components commit it was built from. JBrowse's own tests hold the
bundle, run with no web globals, to the browser's engine.

To rebuild it from a jbrowse-components checkout beside this repo:

```console
cd ../jbrowse-components && pnpm install && cd -
scripts/build-js.sh
```

`scripts/check-against-jbrowse.sh` rebuilds, installs and tests in one step,
against `$JBROWSE_COMPONENTS` or `../jbrowse-components`; jbrowse-components'
downstream canary runs it against its main nightly, so a change there that
breaks this package shows up there. The `Bundle` workflow takes each new
`@jbrowse/tables` release from npm (`scripts/take-release.sh`) into the
committed bundle once R CMD check passes.

## Not yet

Feature tracks (GFF3, BED, BigBed gene models), BigWig and VCF; insertions,
clips, modifications and per-base colors on reads; a legend for a linear
color ramp. A CRAN release, which needs the bundle's own dependencies credited
(`cph` and a LICENSE.note) and its source named.

## See also

- [JBrowse 2](https://jbrowse.org)
- [JBrowseR](https://github.com/GMOD/JBrowseR) - JBrowse 2 as an R htmlwidget
- [ggbandage](https://github.com/GMOD/ggbandage) - pangenome graphs as ggplot2 layers
- [ggtubemap](https://github.com/GMOD/ggtubemap) - pangenome tubemaps as ggplot2 layers
