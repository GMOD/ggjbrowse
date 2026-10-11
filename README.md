# ggjbrowse


> [!CAUTION]
> Very early results from work in progress!

Genome browser tracks as data frames for ggplot2. `jbrowse_alignments()` reads
a region of a BAM or CRAM through the [JBrowse 2](https://jbrowse.org) genome
browser's own data engine, run in V8, and returns what its alignments display
draws: each read on its pileup row, filled as the display's `color` setting
paints it, with its mismatches, deletions, insertions, clips and splice gaps,
the coverage, the splice junctions and the display's key. `autoplot()` draws
the browser's track from them, and a plot is an ordinary ggplot: themes,
scales, `patchwork` and `ggsave()` all work.

```r
remotes::install_github("GMOD/ggjbrowse")
```

```r
library(ggjbrowse)

a <- jbrowse_alignments("rnaseq.bam", "ctgA", 400, 1300, color = "strand")
autoplot(a)
```

![RNA-seq reads over two exons of volvox ctgA: junction arcs, the coverage band and the pileup, split at the introns](man/figures/spliced.png)

The track is the display's stack: splice-junction arcs where the reads
splice, each counted and colored by strand; the coverage band with every
mismatch position's alleles stacked in their base colors; and the pileup,
each read filled as the browser fills it, with its mismatches, deletions,
insertions and soft or hard clips, each faded as the browser fades it where
it is narrower than a pixel and rarer than the depth there warrants, so
sequencing errors recede and shared variants stand. Zoomed in far enough, the
bases and insertion lengths are written out.

```r
l <- jbrowse_alignments("long-reads.bam", "ctgA", 1000, 2700, color = "strand")
autoplot(l)
```

![Long reads over a volvox structural variant: soft clips at the breakpoints, insertions along the reads](man/figures/long-reads.png)

```r
autoplot(jbrowse_alignments("reads.bam", "ctgA", 1000, 1120, color = "mapq"))
```

![120 bp of volvox reads, colored by mapping quality, with the mismatched bases written out](man/figures/zoom.png)

```r
m <- jbrowse_alignments("ont.bam", "20", 10000, 22000, reference = "hg38.fa",
  base_color = "modifications", modifications = list(unmodified = "calls"))
autoplot(m)
```

![Nanopore reads over 12 kb of chr20, each CpG call painted red where methylated and blue where not](man/figures/methylation.png)

Each band is its own ggplot too (`autoplot(a, "pileup")`,
`"coverage"`, `"junctions"`), to restyle or stack with patchwork.

The tables are plain data frames, so a figure can be any ggplot over them:

```r
library(ggplot2)

cram <- jbrowse_alignments("reads.cram", "chr1", 1e6, 1.002e6,
  reference = "hg38.fa")
ggplot(cram$mismatches, aes(position, frequency, colour = base)) +
  geom_point() +
  scale_colour_manual(values = jbrowse_colors()$bases)
```

`color`, `filter`, `base_color`, `modifications` and `sort` are the
display's own settings objects, so what JBrowse can color, filter, paint or
sort by, this can: `color = "mapq"`,
`color = list(field = "insertSize", domain = c(200, 600))`,
`filter = list(flagExclude = 1796, tagFilters = list(list(tag = "HP", value = "1")))`,
`base_color = "modifications"`, `sort = list(type = "basePair", pos = 1055)`.
Each read comes back on the row the display lays it on, with the fill the
display paints it and the bucket it files it under, each modification call
with its color and opacity, the coverage band's stacks as the display stacks
them, and `a$legend` is the display's key.

## Feature tracks

`jbrowse_features()` reads a GFF3, GTF, BED or BigBed region through
JBrowse's feature display and returns what it draws: each feature on its
row, its gene glyph's boxes, intron lines, strand chevrons and arrows in
their colors, its transcripts, and the labels the display shows at that
zoom. `display` takes the display's own settings (`color`, `displayMode`,
`showLabels`, `geneGlyphMode`, `showOnlyGenes`, ...), and `width` is the
pixel width the display lays the region out across, which decides which
names fit. Draw it over reads with `features`:

```r
genes <- jbrowse_features("genes.gff3.gz", "ctgA", 1000, 10000,
  display = list(showOnlyGenes = TRUE))
reads <- jbrowse_alignments("rnaseq.bam", "ctgA", 1000, 10000,
  color = "strand")
autoplot(reads, features = genes)
```

![EDEN's three isoforms and two smaller genes over spliced RNA-seq reads on volvox ctgA](man/figures/genes-over-reads.png)

The y axis of a feature track is the display's pixels. Labels scale with the
panel they land in, so they keep to the room the layout left them at any
output size.

A file can be a path or a URL; a remote file is read in byte ranges through
`curl`. CRAM needs the reference, as a FASTA with its `.fai` or a 2bit. There
is no Bioconductor dependency: the BAM, CRAM, FASTA, GFF3, BED and BigBed
readers are JBrowse's.

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

BigWig and VCF; per-base quality and bisulfite layers; read-pair arcs; a
legend for a linear color ramp; facets on a feature track; the amino-acid
overlay on a zoomed-in CDS. A CRAN release, which needs the bundle's own dependencies credited
(`cph` and a LICENSE.note) and its source named.

## See also

- [JBrowse 2](https://jbrowse.org)
- [JBrowseR](https://github.com/GMOD/JBrowseR) - JBrowse 2 as an R htmlwidget
- [ggbandage](https://github.com/GMOD/ggbandage) - pangenome graphs as ggplot2 layers
- [ggtubemap](https://github.com/GMOD/ggtubemap) - pangenome tubemaps as ggplot2 layers
