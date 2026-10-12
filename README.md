# ggjbrowse


> [!CAUTION]
> Very early results from work in progress!

Genome browser tracks as ggplot2 layers.

- **In:** a region of a BAM, CRAM, GFF3, GTF, BED, BigBed, BigWig or VCF, or
  two regions a PAF aligns.
- **Out:** the tables the [JBrowse 2](https://jbrowse.org) genome browser's
  own display draws from, run in V8: each read on its pileup row, each gene's
  glyph, each bin's score, each sample's genotype, each synteny ribbon, with
  the color the display paints it.
- **Drawn with:** a geom per table, so a figure is an ordinary ggplot that
  `+` changes piece by piece. `autoplot()` composes the geoms into the
  browser's track.

```r
remotes::install_github("GMOD/ggjbrowse")
```

See [Docs](#docs) for how the engine works.

Every example here runs as written, on test data installed with the
package: clips of public human data, and JBrowse's simulated volvox genome
where an example needs a BigWig or a synteny alignment.
`scripts/extdata-*.sh` make the human clips.

```r
library(ggjbrowse)
e <- function(f) system.file("extdata", f, package = "ggjbrowse")

a <- jbrowse_alignments(e("gapdh_direct_rna.bam"), "chr12:6,534,201-6,538,700")
autoplot(a)
```

A region is a locus as JBrowse's location box writes it, 1-based:
`"chr12:6,534,201-6,538,700"`, or `"chr12", 6534200, 6538700` in 0-based
half-open coordinates.

![Nanopore direct RNA reads from NA12878 across GAPDH: junction arcs counting each splice, the coverage band and the reads split at the introns](man/figures/gapdh.png)

The track is the display's stack: splice-junction arcs where the reads
splice, each counted and colored by strand; the coverage band with every
mismatch position's alleles stacked in their base colors; and the pileup,
each read filled as the browser fills it, with its mismatches, deletions,
insertions and soft or hard clips, each faded as the browser fades it where
it is narrower than a pixel and rarer than the depth there warrants, so
sequencing errors recede and shared variants stand. Zoomed in far enough, the
bases and insertion lengths are written out.

```r
l <- jbrowse_alignments(e("colo829_tumor_chr3.bam"),
  "chr3:25,357,600-25,361,000",
  display = list(color = list(field = "strand")))
autoplot(l)
```

![COLO829 nanopore reads over two chr3 breakpoints: reads clipped where the tumor chromosome continues on chr12 and chr10](man/figures/long-reads.png)

```r
d <- jbrowse_alignments(e("hg008t_cuzd1_deletion.bam"),
  "chr10:122,831,700-122,840,800",
  display = list(sort = list(type = "basePair", pos = 122836300)))
autoplot(d)
```

![HG008-T HiFi reads over CUZD1, sorted by the base at the center so the reads with the 1.8 kb somatic deletion gather in one block](man/figures/deletion.png)

```r
autoplot(jbrowse_alignments(e("hg008t_cuzd1_deletion.bam"),
  "chr10:122,835,300-122,835,420",
  display = list(sort = list(type = "basePair", pos = 122835344))))
```

![120 bp at the deletion's left edge: most reads draw the deletion as a line from the left edge, and a heterozygous G is written out as a base](man/figures/zoom.png)

```r
m <- jbrowse_alignments(e("methylation_clip.bam"), "20:10,001-22,000",
  reference = e("hg38_clip.fa"),
  display = list(baseColor = list(field = "modifications"),
    modifications = list(unmodified = "all")))
autoplot(m)
```

![Nanopore reads over 12 kb of chr20, each CpG call painted red where methylated and blue where not](man/figures/methylation.png)

`readConnections = "arc"` adds the display's read-pair band: each mate
pair's and split read's arc, colored by insert size and pair orientation.
`"cloud"` plots each abnormal pair at its insert size instead. A `facet` splits
the reads into the display's stacked sections, each packed separately:

```r
inv <- jbrowse_alignments(e("hg00151_chr1_inversion.bam"),
  "chr1:197,786,900-197,789,700",
  display = list(readConnections = "arc", color = list(field = "strand"),
    facet = list(field = "split")))
autoplot(inv)
```

![HG00151 nanopore reads over a 1.2 kb chr1 inversion, split by whether a read has an SA tag: the split reads' arcs join the breakpoints, a reverse-strand core running between forward flanks](man/figures/inversion.png)

Each band is its own ggplot too (`autoplot(a, "pileup")`,
`"coverage"`, `"junctions"`, `"arcs"`), to restyle or stack with patchwork.

## The grammar

`autoplot()` only composes pieces you can use yourself:

- **Data:** `ggplot(a)` holds a result, and each table is a data frame
  (`a$reads`, `a$mismatches`, `a$coverage`, ...).
- **Geoms:** one per table: `geom_jbrowse_reads()`,
  `geom_jbrowse_mismatches()`, `geom_jbrowse_coverage()`,
  `geom_jbrowse_arcs()`, `geom_jbrowse_glyphs()`, `geom_jbrowse_wiggle()`,
  `geom_jbrowse_genotypes()` and the rest. Each draws as the display does
  and takes any column of its table in `aes()`.
- **Scales:** each table's `color` column is the display's color, which
  ggplot paints through an identity scale. Map another column and ggplot's
  own scales take over. `scale_fill_jbrowse(a)` adds the display's key,
  `scale_y_jbrowse(w)` a quantitative display's score axis, and
  `scale_x_jbrowse(a)` the axis of several regions laid end to end, which
  `ggplot(a)` adds unasked, moving the result's tables in your own layers
  onto it.
- **Facets:** the `section` column is the display's facet, so
  `facet_grid(section ~ .)` splits a plot the way `facet =` splits the track.

```r
library(ggplot2)

ggplot(a) +
  geom_jbrowse_reads(aes(fill = mapq)) +
  geom_jbrowse_mismatches() +
  scale_fill_viridis_c() +
  scale_y_reverse()
```

The tables work without the geoms as well:

```r
cram <- jbrowse_alignments(e("volvox-sorted.cram"), "ctgA:1,001-3,000",
  reference = e("volvox.fa"))
ggplot(cram$mismatches, aes(position, frequency, colour = color)) +
  geom_point() +
  scale_colour_jbrowse(cram, "bases")
```

`display` is the display's config, written as a JBrowse track config writes
it, so what JBrowse can color, filter, paint, sort or split by, this can:
`color = list(field = "mapq")`,
`color = list(field = "insertSize", domain = c(200, 600))`,
`filter = list(flagExclude = 1796, tagFilters = list(list(tag = "HP", value = "1")))`,
`baseColor = list(field = "modifications")`,
`sort = list(type = "basePair", pos = 1055)`,
`facet = list(field = "tags.HP")`, or any other slot such as
`showOrdinaryPairs = FALSE`. As in a track config, a bare string `color` is
one constant color: `color = "purple"`.
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
names fit. Draw it above reads with `above`:

```r
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"),
  "chr12:6,534,201-6,538,700")
autoplot(a, above = genes)
```

![GENCODE's GAPDH transcripts over the direct RNA reads, each read's exons lining up under the transcripts' boxes](man/figures/genes-over-reads.png)

The y axis of a feature track is the display's pixels. Labels scale with the
panel they land in, so they keep to the room the layout left them at any
output size.

## Quantitative tracks

`jbrowse_wiggle()` reads a BigWig through JBrowse's quantitative display: the
bins the file serves at that zoom, the marks the display draws from them in
its colors, and the score axis it autoscales to. `display` takes its settings
(`mark`, `aggregate`, `extent`, `color`, `scales`, ...), and `width` picks the
zoom level that answers. Zoomed out, a min-max extent draws each bin's max and
min around its mean:

```r
bw <- e("volvox_microarray_posneg.bw")
autoplot(jbrowse_wiggle(bw, "ctgA:1-50,000", width = 100))
autoplot(jbrowse_wiggle(bw, "ctgA:1-50,000", width = 100,
  display = list(mark = "line", interpolate = "linear")))
autoplot(jbrowse_wiggle(bw, "ctgA:1-50,000", display = list(mark = "span")))
```

![A signed volvox microarray track three ways: min-max bars, an interpolated line over its min-max band, and density](man/figures/wiggle.png)

A track of several BigWigs, read through a `MultiWiggleAdapter` config,
draws a row per source in `w$sources`' order. Reorder or filter that table
to rearrange them, as for a variant matrix's rows:

- `as.matrix(w)` is the matrix the display's "Cluster rows" clusters: a row
  per source, a column per pixel.
- `w$tree <- hclust(dist(as.matrix(w)), "average")` puts the rows in the
  tree's order, and `autoplot()` draws the tree beside them.

## Variant matrices

`jbrowse_variants()` reads a multi-sample VCF through JBrowse's variant
display: the records its filters keep, each sample's genotype at each, and
the fill the display paints it, its hue shaded by the alt dosage. `display`
takes the display's settings (`unit = "haplotype"` for a row per haplotype,
`color`, `filter`, `minorAlleleFrequencyFilter`, ...), and `samples` a table
of the samples' attributes, which `rowColor` marks the rows by:

```r
v <- jbrowse_variants(e("1000g_dcst2.vcf.gz"), "1:155,000,001-155,050,000",
  samples = e("1000g_dcst2_samples.tsv"),
  display = list(rowColor = list(field = "superpopulation")))
autoplot(v, layout = "columns")
```

![100 people of the 1000 Genomes, 20 from each superpopulation, genotyped at 249 variants across DCST2 (GRCh37), a column per record, each row marked by its superpopulation](man/figures/variants.png)

The matrix draws its rows in `v$rows`' order: reorder or filter that table
to rearrange them. A tree in `v$tree` puts the rows in its leaf order, and
`autoplot()` draws it beside them (`geom_jbrowse_tree()`). Any clustering
that gives an `hclust`, or anything `as.hclust()` takes, labelled by row
name, will do.

`as.matrix(v)` is the genotype matrix the display's "Cluster rows by
genotype" clusters, a column per ALT allele with `NA` for a no-call, and
this is the R script its dialog writes:

```r
v$tree <- hclust(dist(as.matrix(v)), "average")
autoplot(v, layout = "columns")
```

![The same people clustered: most East Asians gather near the top, the Africans and Europeans further down, the tree beside the rows](man/figures/variants-clustered.png)

## A view of several tracks

`jbrowse_stack()` stacks tracks as a JBrowse view does: each track's plot,
or an alignments track's bands, top to bottom in the order given, on one x
axis, with the keys collected. Several loci draw end to end, as a view lays
out several regions, each tick labelled in its own region's positions:

```r
loci <- c("ctgA:1,001-4,000", "ctgA:6,001-9,000")
genes <- jbrowse_features(e("volvox.sort.gff3.gz"), loci)
signal <- jbrowse_wiggle(e("volvox_microarray_posneg.bw"), loci)
reads <- jbrowse_alignments(e("volvox-sorted.bam"), loci,
  display = list(color = list(field = "strand")))
jbrowse_stack(genes, signal, reads)
```

![Volvox genes, a BigWig and strand-colored reads over two regions of ctgA laid end to end](man/figures/stack.png)

Each band is an ordinary ggplot. Restyle every band with patchwork's `&`,
`jbrowse_stack(genes, signal, reads) & theme_classic()`, or draw one on its
own and add to it: `autoplot(reads, "pileup") + theme_classic()`. A ggplot
of your own stacks as a band on the same axis.
A track's settings decide its bands: an alignments track read with
`display = list(showCoverage = FALSE)` stacks no coverage band.

Each band takes the pixels the display gives it. Reads draw at the display's
`featureHeight` (7, or 3 and 1 for its compact and super-compact settings),
and `jbrowse_height()` gives the figure height that keeps those proportions,
so a compact pileup saves as a shorter figure:

```r
compact <- jbrowse_alignments(e("volvox-sorted.bam"), loci,
  display = list(featureHeight = 3, color = list(field = "strand")))
fig <- jbrowse_stack(genes, signal, compact)
ggsave("stack.png", fig, width = 10, height = jbrowse_height(fig))
```

A file can be a path or a URL; a remote file is read in byte ranges through
`curl`. CRAM needs the reference: a FASTA with its `.fai`, a bgzipped one
with its `.fai` and `.gzi`, or a 2bit. There is no Bioconductor dependency:
the BAM, CRAM, FASTA, GFF3, BED, BigBed, BigWig and VCF readers are
JBrowse's.

## Synteny

`jbrowse_synteny()` reads a pairwise alignment (PAF, indexed PIF, MUMmer
delta, UCSC chain, MashMap) through JBrowse's linear synteny display: the
blocks two rows of regions share, and the ribbons the band draws between
them in the view's colors and opacity, with each CIGAR insertion, deletion
and skip wider than a pixel drawn in its color. `view` takes the
`LinearSyntenyView`'s settings (`color`, `opacity`, `cigarMode`,
`minAlignmentLength`, `drawLocationMarkers`, ...). Tracks over either row
stack above and below the band:

```r
paf <- e("volvox_inv_indels.paf")
s <- jbrowse_synteny(paf, "ctgA:1-49,186", "ctgA:1-50,001",
  view = list(color = list(field = "strand")))
genes <- jbrowse_features(e("volvox.sort.gff3.gz"), "ctgA:1-50,001",
  display = list(showOnlyGenes = TRUE))
autoplot(s, below = genes)
```

![A simulated volvox with two inversions over volvox ctgA: forward blocks in tomato, the inverted ones crossing in blue, the genes under the target row](man/figures/synteny.png)

Zoomed in, the CIGAR shows: insertions in yellow, deletions in blue over the
red match block, as the band paints them:

```r
autoplot(jbrowse_synteny(paf, "ctgA:26,801-30,000", "ctgA:27,201-30,500",
  width = 1500))
```

![3 kb of the same alignment: the match block sheared by its insertions and deletions](man/figures/synteny-zoom.png)

`jbrowse_dotplot()` reads the same files through JBrowse's dotplot display,
a line per block, or per CIGAR step where a block is wide enough on the
plot, and `view` takes the `DotplotView`'s settings:

```r
autoplot(jbrowse_dotplot(paf, "ctgA:1-49,186", "ctgA:1-50,001",
  view = list(color = list(field = "strand"))))
```

![The same alignment as a dotplot: a diagonal of forward blocks, two short reverse ones](man/figures/dotplot.png)

A continuous field such as `list(field = "identity")` keys as the display's
bar, `0%` to `100%`, through `guides_jbrowse()`; a threshold field such as
`mapq` keys its bins as chips.

The tables are `s$alignments`, `s$ribbons` (each foot's span on its row, and
`thin` where the band draws a one-pixel line) and `s$offscreen` (blocks whose
mate the other row does not show), and `geom_jbrowse_ribbons()` and
`geom_jbrowse_dotplot()` draw them. `geom_jbrowse_offscreen_mates()` draws
`s$offscreen_strip`, the marks and contig names the band hangs off a row's
edge for those blocks, and `autoplot()` keys each contig with the line the
strip's hover prints, from `s$offscreen_contigs`. A file guessed from its extension names
its sides `query` and `target`; a JBrowse adapter config names assemblies,
and `assemblies = c(query, target)` picks which goes on which row.
`swap = TRUE` puts the file's target on the top row instead. With no regions
given, a row shows every contig of its assembly, in natural order, and a
contig named alone shows whole, so a whole-genome dotplot is one call:

```r
autoplot(jbrowse_dotplot(e("peach_grape.paf.gz")))
```

![Peach against grape, every chromosome of each along its axis, named](man/figures/whole-genome.png)

### Several genomes

`jbrowse_multiway()` reads an all-vs-all PAF (sequences named
`sample#haplotype#contig`), an indexed PIF or a multi-genome adapter config
through JBrowse's multi-way synteny display: an anchor lane over the locus
given, a lane per other genome in the alignments' own coordinates, each at
the zoom and orientation the display gives it, and the ribbons between
adjacent lanes. `lanes` picks and orders the genomes, or pins a lane to a
contig, and `display` takes the display's settings (`ribbonColor`,
`drawCurves`, `bridgeSkippedLanes`, ...). Tracks over the anchor lane stack
above it, as `autoplot(m, above = genes)`:

```r
m <- jbrowse_multiway(e("volvox_all_vs_all.paf"), "ctgA:1-50,001", "volvox",
  display = list(ribbonColor = list(field = "strand")))
autoplot(m)
```

![Volvox over its insertion and deletion simulations: three lanes, each with its own ruler, the ribbons between them opening at the indels](man/figures/multiway.png)

## From a JBrowse view

- JBrowse's linear genome view is gaining an **Export R script** menu item,
  not yet in a release. It writes one ggjbrowse reader per open track, then
  `jbrowse_stack()` and `ggsave(height = jbrowse_height(fig))`.
- Each figure on the right is that script, downloaded from the view on the
  left and run unchanged with `Rscript`.

<table>
<tr><th width="50%">JBrowse</th><th width="50%">R</th></tr>
<tr><td colspan="2"><b>Genes over a BAM</b> <code>ctgA:1,000-9,000</code></td></tr>
<tr><td width="50%"><img src="gallery/genes-and-reads-jbrowse.webp" alt="Genes over a BAM in JBrowse" width="100%"></td><td width="50%"><img src="gallery/genes-and-reads-r.webp" alt="Genes over a BAM drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Base-level zoom, reads colored by MAPQ</b> <code>ctgA:1,000-1,150</code></td></tr>
<tr><td width="50%"><img src="gallery/mapq-zoom-jbrowse.webp" alt="Base-level zoom, reads colored by MAPQ in JBrowse" width="100%"></td><td width="50%"><img src="gallery/mapq-zoom-r.webp" alt="Base-level zoom, reads colored by MAPQ drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Compact reads over a deep region</b> <code>ctgA:1-20,000</code></td></tr>
<tr><td width="50%"><img src="gallery/compact-reads-jbrowse.webp" alt="Compact reads over a deep region in JBrowse" width="100%"></td><td width="50%"><img src="gallery/compact-reads-r.webp" alt="Compact reads over a deep region drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Super-compact long reads</b> <code>ctgA:1,000-2,700</code></td></tr>
<tr><td width="50%"><img src="gallery/super-compact-reads-jbrowse.webp" alt="Super-compact long reads in JBrowse" width="100%"></td><td width="50%"><img src="gallery/super-compact-reads-r.webp" alt="Super-compact long reads drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>RNA-seq with splice junctions</b> <code>ctgA:1,000-10,000</code></td></tr>
<tr><td width="50%"><img src="gallery/rnaseq-sashimi-jbrowse.webp" alt="RNA-seq with splice junctions in JBrowse" width="100%"></td><td width="50%"><img src="gallery/rnaseq-sashimi-r.webp" alt="RNA-seq with splice junctions drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Reads split into strand sections</b> <code>ctgA:2,000-5,000</code></td></tr>
<tr><td width="50%"><img src="gallery/facet-strand-jbrowse.webp" alt="Reads split into strand sections in JBrowse" width="100%"></td><td width="50%"><img src="gallery/facet-strand-r.webp" alt="Reads split into strand sections drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Paired reads across an inversion</b> <code>ctgA:18,000-32,000</code></td></tr>
<tr><td width="50%"><img src="gallery/inversion-arcs-jbrowse.webp" alt="Paired reads across an inversion in JBrowse" width="100%"></td><td width="50%"><img src="gallery/inversion-arcs-r.webp" alt="Paired reads across an inversion drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Long reads split across a structural variant</b> <code>ctgA:1,800-3,000 ctgA:33,900-35,100</code></td></tr>
<tr><td width="50%"><img src="gallery/long-reads-sv-jbrowse.webp" alt="Long reads split across a structural variant in JBrowse" width="100%"></td><td width="50%"><img src="gallery/long-reads-sv-r.webp" alt="Long reads split across a structural variant drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>CpG methylation from MM/ML tags</b> <code>20:10,000-22,000</code></td></tr>
<tr><td width="50%"><img src="gallery/methylation-jbrowse.webp" alt="CpG methylation from MM/ML tags in JBrowse" width="100%"></td><td width="50%"><img src="gallery/methylation-r.webp" alt="CpG methylation from MM/ML tags drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Multi-BigWig rows</b> <code>ctgA:1-50,000</code></td></tr>
<tr><td width="50%"><img src="gallery/multi-bigwig-jbrowse.webp" alt="Multi-BigWig rows in JBrowse" width="100%"></td><td width="50%"><img src="gallery/multi-bigwig-r.webp" alt="Multi-BigWig rows drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>200 BigWig rows</b> <code>ctgA:1-50,000</code></td></tr>
<tr><td width="50%"><img src="gallery/many-rows-200-jbrowse.webp" alt="200 BigWig rows in JBrowse" width="100%"></td><td width="50%"><img src="gallery/many-rows-200-r.webp" alt="200 BigWig rows drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>1,000 BigWig rows as density</b> <code>ctgA:1-50,000</code></td></tr>
<tr><td width="50%"><img src="gallery/many-rows-1000-jbrowse.webp" alt="1,000 BigWig rows as density in JBrowse" width="100%"></td><td width="50%"><img src="gallery/many-rows-1000-r.webp" alt="1,000 BigWig rows as density drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Quantitative renderings</b> <code>ctgA:1-50,000</code></td></tr>
<tr><td width="50%"><img src="gallery/wiggle-renderings-jbrowse.webp" alt="Quantitative renderings in JBrowse" width="100%"></td><td width="50%"><img src="gallery/wiggle-renderings-r.webp" alt="Quantitative renderings drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Two regions end to end</b> <code>ctgA:1,001-4,000 ctgA:6,001-9,000</code></td></tr>
<tr><td width="50%"><img src="gallery/two-regions-jbrowse.webp" alt="Two regions end to end in JBrowse" width="100%"></td><td width="50%"><img src="gallery/two-regions-r.webp" alt="Two regions end to end drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Variants beside genes</b> <code>ctgA:1-25,000</code></td></tr>
<tr><td width="50%"><img src="gallery/variants-and-genes-jbrowse.webp" alt="Variants beside genes in JBrowse" width="100%"></td><td width="50%"><img src="gallery/variants-and-genes-r.webp" alt="Variants beside genes drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Genes colored by a config expression</b> <code>ctgA:1-25,000</code></td></tr>
<tr><td width="50%"><img src="gallery/jexl-colors-jbrowse.webp" alt="Genes colored by a config expression in JBrowse" width="100%"></td><td width="50%"><img src="gallery/jexl-colors-r.webp" alt="Genes colored by a config expression drawn by ggjbrowse" width="100%"></td></tr>
<tr><td colspan="2"><b>Structural variants across 1,800 dogs</b> <code>chr6:47,340,000-47,420,000</code></td></tr>
<tr><td width="50%"><img src="gallery/dog-sv-matrix-jbrowse.webp" alt="Structural variants across 1,800 dogs in JBrowse" width="100%"></td><td width="50%"><img src="gallery/dog-sv-matrix-r.webp" alt="Structural variants across 1,800 dogs drawn by ggjbrowse" width="100%"></td></tr>
</table>

## The engine

- `inst/js/jbrowse-tables.js` is
  [`@jbrowse/tables`](https://github.com/GMOD/jbrowse-components/tree/main/packages/tables)
  bundled for a bare JavaScript host, committed so installing needs no
  JavaScript toolchain.
- JBrowse's own tests hold the bundle, run with no web globals, to the
  browser's engine.
- [CONTRIBUTING.md](CONTRIBUTING.md#updating-the-engine) covers rebuilding it
  and the CI that takes each new release.

## Not yet

The variant display's row arrangements other than clustering (bands, a
saved order) and its zoomed-out density fade; the bisulfite layer; facets on a feature
track; the amino-acid
overlay on a zoomed-in CDS. A CRAN release, which needs the bundle's own dependencies credited
(`cph` and a LICENSE.note) and its source named.

## Docs

- [docs/api.md](docs/api.md): every reader, helper and geom, and what each
  returns
- [docs/how-does-this-work.md](docs/how-does-this-work.md): the V8 engine, how
  it reads files through R, and which side decides colors and geometry
- [docs/dataflow.md](docs/dataflow.md): a diagram of one call from path to plot
- [docs/faq.md](docs/faq.md): V8 versus QuickJS, caching, indexes, coordinates
- [agent-docs/engine-owns-display-decisions.md](agent-docs/engine-owns-display-decisions.md):
  the plan to move the rules R still copies into the engine
- [CONTRIBUTING.md](CONTRIBUTING.md): development, rebuilding the engine and CI

## See also

- [JBrowse 2](https://jbrowse.org)
- [JBrowseR](https://github.com/GMOD/JBrowseR) - JBrowse 2 as an R htmlwidget
- [ggbandage](https://github.com/GMOD/ggbandage) - pangenome graphs as ggplot2 layers
- [ggtubemap](https://github.com/GMOD/ggtubemap) - pangenome tubemaps as ggplot2 layers

## Academic use

Written with [NHGRI](http://genome.gov) funding as part of
[JBrowse](http://jbrowse.org). If you use this in a publication, please cite the
most recent JBrowse paper at [jbrowse.org](http://jbrowse.org).

## License

Apache License 2.0 © [Colin Diesh](https://github.com/cmdcolin)
