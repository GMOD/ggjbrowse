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
- **Also:** reads CRAM, which Rsamtools cannot, and needs no Bioconductor
  packages.
- **Drawn with:** a geom per table, so a figure is an ordinary ggplot that
  `+` changes piece by piece. `autoplot()` composes the geoms into the
  browser's track.

```r
remotes::install_github("GMOD/ggjbrowse")
```

## Gallery

- Each row draws one figure twice from data installed with the package, with
  JBrowse's own renderer (`jb2export` from `@jbrowse/img`) on the left and
  ggjbrowse on the right.
- Each title links to the figure's R and `jb2export` recipes in
  [docs/gallery.md](docs/gallery.md).

<!-- gallery -->
<table>
<tr><th width="50%">JBrowse</th><th width="50%">R</th></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#genes_and_reads">GENCODE genes over direct RNA reads</a></b> <code>chr12:6,534,201-6,538,700</code></td></tr>
<tr><td width="50%"><img src="gallery/genes_and_reads-jbrowse.webp" alt="GENCODE genes over direct RNA reads in JBrowse. Junction arcs count each splice, from 31 to 68 reads, over the coverage band and the reads split at every intron." width="100%"></td><td width="50%"><img src="gallery/genes_and_reads-r.webp" alt="GENCODE genes over direct RNA reads drawn by ggjbrowse. Junction arcs count each splice, from 31 to 68 reads, over the coverage band and the reads split at every intron." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#long_reads_sv">COLO829 der(3) split reads chained across its two chr3 breakpoints</a></b> <code>chr3:25,358,500-25,359,800</code></td></tr>
<tr><td width="50%"><img src="gallery/long_reads_sv-jbrowse.webp" alt="COLO829 der(3) split reads chained across its two chr3 breakpoints in JBrowse. Coverage halves at the first breakpoint, and each read's segments share a row, overlapping in dark grey and colored by strand past it, joined by duplication-type arcs." width="100%"></td><td width="50%"><img src="gallery/long_reads_sv-r.webp" alt="COLO829 der(3) split reads chained across its two chr3 breakpoints drawn by ggjbrowse. Coverage halves at the first breakpoint, and each read's segments share a row, overlapping in dark grey and colored by strand past it, joined by duplication-type arcs." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#deletion">Reads sorted by the base at the center, over a somatic deletion</a></b> <code>chr10:122,831,700-122,840,800</code></td></tr>
<tr><td width="50%"><img src="gallery/deletion-jbrowse.webp" alt="Reads sorted by the base at the center, over a somatic deletion in JBrowse. Coverage halves across the deletion, and the reads sorted to the top span it with dark deletion bars, each labelled with its length." width="100%"></td><td width="50%"><img src="gallery/deletion-r.webp" alt="Reads sorted by the base at the center, over a somatic deletion drawn by ggjbrowse. Coverage halves across the deletion, and the reads sorted to the top span it with dark deletion bars, each labelled with its length." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#compact_reads">Compact reads over the same deletion</a></b> <code>chr10:122,831,700-122,840,800</code></td></tr>
<tr><td width="50%"><img src="gallery/compact_reads-jbrowse.webp" alt="Compact reads over the same deletion in JBrowse. At 3 pixels a row the pileup takes far less height, the deletion's dark bars banded through it." width="100%"></td><td width="50%"><img src="gallery/compact_reads-r.webp" alt="Compact reads over the same deletion drawn by ggjbrowse. At 3 pixels a row the pileup takes far less height, the deletion's dark bars banded through it." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#base_zoom">Base-level zoom at the deletion's left edge</a></b> <code>chr10:122,835,300-122,835,420</code></td></tr>
<tr><td width="50%"><img src="gallery/base_zoom-jbrowse.webp" alt="Base-level zoom at the deletion's left edge in JBrowse. At 120 bp the pileup writes out each mismatch's base, and every read sorted to the top carries a G a few bases before its deletion." width="100%"></td><td width="50%"><img src="gallery/base_zoom-r.webp" alt="Base-level zoom at the deletion's left edge drawn by ggjbrowse. At 120 bp the pileup writes out each mismatch's base, and every read sorted to the top carries a G a few bases before its deletion." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#inversion">Split reads chained across an inversion, with their arcs</a></b> <code>chr1:197,786,900-197,789,700</code></td></tr>
<tr><td width="50%"><img src="gallery/inversion-jbrowse.webp" alt="Split reads chained across an inversion, with their arcs in JBrowse. The split reads, in their own section above the unsplit ones, are colored inverted between the two breakpoints, where inversion-type arcs join their segments." width="100%"></td><td width="50%"><img src="gallery/inversion-r.webp" alt="Split reads chained across an inversion, with their arcs drawn by ggjbrowse. The split reads, in their own section above the unsplit ones, are colored inverted between the two breakpoints, where inversion-type arcs join their segments." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#methylation">CpG methylation from MM/ML tags</a></b> <code>20:10,001-22,000</code></td></tr>
<tr><td width="50%"><img src="gallery/methylation-jbrowse.webp" alt="CpG methylation from MM/ML tags in JBrowse. Each read's 5mC calls in red where methylated and blue where not." width="100%"></td><td width="50%"><img src="gallery/methylation-r.webp" alt="CpG methylation from MM/ML tags drawn by ggjbrowse. Each read's 5mC calls in red where methylated and blue where not." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#dnase_rows">DNase in twelve cell lines over the beta-globin locus</a></b> <code>chr11:5,270,001-5,320,000</code></td></tr>
<tr><td width="50%"><img src="gallery/dnase_rows-jbrowse.webp" alt="DNase in twelve cell lines over the beta-globin locus in JBrowse. Most lines share one peak near 5,281,000, and K562 alone shows several more to either side." width="100%"></td><td width="50%"><img src="gallery/dnase_rows-r.webp" alt="DNase in twelve cell lines over the beta-globin locus drawn by ggjbrowse. Most lines share one peak near 5,281,000, and K562 alone shows several more to either side." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#dnase_clustered">DNase in twelve cell lines as density, clustered by signal</a></b> <code>chr11:5,270,001-5,320,000</code></td></tr>
<tr><td width="50%"><img src="gallery/dnase_clustered-jbrowse.webp" alt="DNase in twelve cell lines as density, clustered by signal in JBrowse. K562 splits off alone from the other eleven lines." width="100%"></td><td width="50%"><img src="gallery/dnase_clustered-r.webp" alt="DNase in twelve cell lines as density, clustered by signal drawn by ggjbrowse. K562 splits off alone from the other eleven lines." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#conservation">phyloP conservation across GAPDH as bars, a line and density</a></b> <code>chr12:6,534,201-6,538,700</code></td></tr>
<tr><td width="50%"><img src="gallery/conservation-jbrowse.webp" alt="phyloP conservation across GAPDH as bars, a line and density in JBrowse. The exons stand out as runs of high scores in bars, in a mean line over its min-max band, and in a density strip." width="100%"></td><td width="50%"><img src="gallery/conservation-r.webp" alt="phyloP conservation across GAPDH as bars, a line and density drawn by ggjbrowse. The exons stand out as runs of high scores in bars, in a mean line over its min-max band, and in a density strip." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#two_regions">Genes, conservation and strand-colored reads over two regions end to end</a></b> <code>chr12:6,534,201-6,535,700 chr12:6,537,201-6,538,700</code></td></tr>
<tr><td width="50%"><img src="gallery/two_regions-jbrowse.webp" alt="Genes, conservation and strand-colored reads over two regions end to end in JBrowse. A grey divider separates the two loci, and each tick gives a position in its own region." width="100%"></td><td width="50%"><img src="gallery/two_regions-r.webp" alt="Genes, conservation and strand-colored reads over two regions end to end drawn by ggjbrowse. A grey divider separates the two loci, and each tick gives a position in its own region." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#variants_by_population">100 people's genotypes across DCST2 (GRCh37), colored by superpopulation</a></b> <code>1:155,000,001-155,050,000</code></td></tr>
<tr><td width="50%"><img src="gallery/variants_by_population-jbrowse.webp" alt="100 people's genotypes across DCST2 (GRCh37), colored by superpopulation in JBrowse. 20 people from each superpopulation genotyped at 249 variants, a column per record, the East Asians homozygous for the alt allele at many of the first 40." width="100%"></td><td width="50%"><img src="gallery/variants_by_population-r.webp" alt="100 people's genotypes across DCST2 (GRCh37), colored by superpopulation drawn by ggjbrowse. 20 people from each superpopulation genotyped at 249 variants, a column per record, the East Asians homozygous for the alt allele at many of the first 40." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#variants_clustered">100 people's genotypes across DCST2 clustered, the tree beside the rows</a></b> <code>1:155,000,001-155,050,000</code></td></tr>
<tr><td width="50%"><img src="gallery/variants_clustered-jbrowse.webp" alt="100 people's genotypes across DCST2 clustered, the tree beside the rows in JBrowse. Most East Asians gather near the top, the Africans and Europeans further down." width="100%"></td><td width="50%"><img src="gallery/variants_clustered-r.webp" alt="100 people's genotypes across DCST2 clustered, the tree beside the rows drawn by ggjbrowse. Most East Asians gather near the top, the Africans and Europeans further down." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#collapsed_introns">GAPDH's nine exons side by side, with the junctions between them as arcs</a></b> <code>ENST00000229239.10's exons, each padded by 20 bp</code></td></tr>
<tr><td width="50%"><img src="gallery/collapsed_introns-jbrowse.webp" alt="GAPDH's nine exons side by side, with the junctions between them as arcs in JBrowse. Each junction's read count sits over the arc joining its two exons, and the direct RNA reads run across." width="100%"></td><td width="50%"><img src="gallery/collapsed_introns-r.webp" alt="GAPDH's nine exons side by side, with the junctions between them as arcs drawn by ggjbrowse. Each junction's read count sits over the arc joining its two exons, and the direct RNA reads run across." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#inversion_synteny">HG002's maternal chr8 over its paternal copy at 8p23.1</a></b> <code>chr8_MATERNAL:5,250,001-14,250,000 chr8_PATERNAL:5,250,001-14,250,000</code></td></tr>
<tr><td width="50%"><img src="gallery/inversion_synteny-jbrowse.webp" alt="HG002's maternal chr8 over its paternal copy at 8p23.1 in JBrowse. Forward blocks in tomato flank the 4 Mb inversion, whose reverse ribbon crosses in blue." width="100%"></td><td width="50%"><img src="gallery/inversion_synteny-r.webp" alt="HG002's maternal chr8 over its paternal copy at 8p23.1 drawn by ggjbrowse. Forward blocks in tomato flank the 4 Mb inversion, whose reverse ribbon crosses in blue." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#inversion_dotplot">HG002's maternal chr8 against its paternal copy as a dotplot</a></b> <code>chr8_MATERNAL chr8_PATERNAL</code></td></tr>
<tr><td width="50%"><img src="gallery/inversion_dotplot-jbrowse.webp" alt="HG002's maternal chr8 against its paternal copy as a dotplot in JBrowse. One forward diagonal runs the length of chr8, broken only by the inversion's short reverse segment near 10 Mb." width="100%"></td><td width="50%"><img src="gallery/inversion_dotplot-r.webp" alt="HG002's maternal chr8 against its paternal copy as a dotplot drawn by ggjbrowse. One forward diagonal runs the length of chr8, broken only by the inversion's short reverse segment near 10 Mb." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#whole_genome">Peach against grape, every chromosome</a></b> <code>whole genome</code></td></tr>
<tr><td width="50%"><img src="gallery/whole_genome-jbrowse.webp" alt="Peach against grape, every chromosome in JBrowse. Short diagonal runs scatter over most pairs of peach and grape chromosomes." width="100%"></td><td width="50%"><img src="gallery/whole_genome-r.webp" alt="Peach against grape, every chromosome drawn by ggjbrowse. Short diagonal runs scatter over most pairs of peach and grape chromosomes." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#amylase_synteny">GRCh38 over HG00133 at the amylase genes</a></b> <code>GRCh38#0#chr1:103,520,894-103,832,637</code></td></tr>
<tr><td width="50%"><img src="gallery/amylase_synteny-jbrowse.webp" alt="GRCh38 over HG00133 at the amylase genes in JBrowse. Forward blocks flank a reversed copy in the middle, with alignments of the other copies crossing over it." width="100%"></td><td width="50%"><img src="gallery/amylase_synteny-r.webp" alt="GRCh38 over HG00133 at the amylase genes drawn by ggjbrowse. Forward blocks flank a reversed copy in the middle, with alignments of the other copies crossing over it." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#amylase_zoom">3 kb at the amylase genes, the match block cut by an 80 bp deletion</a></b> <code>GRCh38#0#chr1:103,619,501-103,622,500</code></td></tr>
<tr><td width="50%"><img src="gallery/amylase_zoom-jbrowse.webp" alt="3 kb at the amylase genes, the match block cut by an 80 bp deletion in JBrowse. The deletion shows as a purple wedge narrowing to a point on HG00133's row." width="100%"></td><td width="50%"><img src="gallery/amylase_zoom-r.webp" alt="3 kb at the amylase genes, the match block cut by an 80 bp deletion drawn by ggjbrowse. The deletion shows as a purple wedge narrowing to a point on HG00133's row." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#synteny_genes">HG002's maternal chr8 over its paternal copy, each with its genes</a></b> <code>chr8_MATERNAL:5,250,001-14,250,000 chr8_PATERNAL:5,250,001-14,250,000</code></td></tr>
<tr><td width="50%"><img src="gallery/synteny_genes-jbrowse.webp" alt="HG002's maternal chr8 over its paternal copy, each with its genes in JBrowse. Forward blocks in tomato, the 4 Mb inversion crossing in blue, and the genes above and below." width="100%"></td><td width="50%"><img src="gallery/synteny_genes-r.webp" alt="HG002's maternal chr8 over its paternal copy, each with its genes drawn by ggjbrowse. Forward blocks in tomato, the 4 Mb inversion crossing in blue, and the genes above and below." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#multiway">The TP53 neighbourhood on hg38 over chimp, gorilla and orangutan</a></b> <code>chr17:7,400,001-7,700,000</code></td></tr>
<tr><td width="50%"><img src="gallery/multiway-jbrowse.webp" alt="The TP53 neighbourhood on hg38 over chimp, gorilla and orangutan in JBrowse. Each lane runs in its own genome's coordinates, every liftOver ribbon in hg38's orientation, with gaps where a genome lacks sequence the one above has." width="100%"></td><td width="50%"><img src="gallery/multiway-r.webp" alt="The TP53 neighbourhood on hg38 over chimp, gorilla and orangutan drawn by ggjbrowse. Each lane runs in its own genome's coordinates, every liftOver ribbon in hg38's orientation, with gaps where a genome lacks sequence the one above has." width="100%"></td></tr>
</table>
<!-- /gallery -->

## Quick start

```r
library(ggjbrowse)
e <- function(f) system.file("extdata", f, package = "ggjbrowse")

a <- jbrowse_alignments(e("gapdh_direct_rna.bam"), "chr12:6,534,201-6,538,700")
autoplot(a)
```

![Nanopore direct RNA reads from NA12878 across GAPDH: junction arcs counting each splice, the coverage band and the reads split at the introns](man/figures/gapdh.png)

- **Data:** every example runs as written on clips of public human data
  installed with the package. `scripts/extdata-*.sh` make the clips.
- **Regions:** a locus as JBrowse's location box writes it, 1-based
  (`"chr12:6,534,201-6,538,700"`), or `"chr12", 6534200, 6538700` in 0-based
  half-open coordinates. A vector of several loci draws them end to end.
- **Files:** a path or a URL, which `curl` reads in byte ranges. CRAM needs
  the reference: a FASTA with its `.fai`, a bgzipped one with its `.fai` and
  `.gzi`, or a 2bit.
- **Bands:** the track stacks the display's splice-junction arcs, its
  coverage band with each mismatch position's alleles, and its pileup, each
  read's mismatches, indels and clips faded where the browser fades them.
  `autoplot(a, "pileup")` draws one band alone, as do `"coverage"`,
  `"junctions"` and `"arcs"`.

## Alignments

`display` is the display's config, written as a JBrowse track config writes
it, so `jbrowse_alignments()` colors, filters, sorts and splits reads as
JBrowse does:

- `color = list(field = "mapq")`, or
  `list(field = "insertSize", domain = c(200, 600))`. A bare string such as
  `color = "purple"` is one constant color.
- `filter = list(flagExclude = 1796, tagFilters = list(list(tag = "HP", value = "1")))`
- `sort = list(type = "basePair", pos = 122836300)` sorts by the base at a
  position ([a somatic deletion](docs/gallery.md#deletion)).
- `featureHeight = 3` packs reads compact
  ([compact reads](docs/gallery.md#compact_reads)).
- `baseColor = list(field = "modifications")` paints MM/ML modification calls
  ([CpG methylation](docs/gallery.md#methylation)).
- `readConnections = "arc"` adds the read-pair band, an arc per mate pair and
  split read colored by insert size and pair orientation. `"cloud"`
  plots each abnormal pair at its insert size instead.
- `unit = "chain"` puts a read's mate or split segments on one row, joined by
  a hairline ([COLO829 split reads](docs/gallery.md#long_reads_sv)).
- `facet = list(field = "split")` splits the reads into the display's stacked
  sections, each packed separately
  ([an inversion](docs/gallery.md#inversion)).
- Any other slot, such as `showOrdinaryPairs = FALSE`.

Each read comes back on the row the display lays it on, with the fill the
display paints it, and `a$legend` is the display's key.

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
cram <- jbrowse_alignments(e("gapdh_direct_rna.cram"),
  "chr12:6,534,201-6,535,200", reference = e("gapdh_hg38.fa.gz"))
ggplot(cram$mismatches, aes(position, frequency, colour = color)) +
  geom_point() +
  scale_colour_jbrowse(cram, "bases")
```

## Feature tracks

`jbrowse_features()` reads a GFF3, GTF, BED or BigBed region through
JBrowse's feature display and returns what it draws:

- each feature on its row, its gene glyph's boxes, intron lines, strand
  chevrons and arrows in their colors, and its transcripts
- the labels the display shows at that zoom, scaled to keep to the room the
  layout left them at any output size

`display` takes the display's settings (`color`, `displayMode`,
`showLabels`, `geneGlyphMode`, `showOnlyGenes`, ...), and `width` is the
pixel width the display lays the region out across, which decides which
names fit. `autoplot(reads, above = genes)` draws genes over reads
([GENCODE genes](docs/gallery.md#genes_and_reads)).

## Quantitative tracks

`jbrowse_wiggle()` reads a BigWig through JBrowse's quantitative display: the
bins the file serves at that zoom, the marks the display draws from them, and
the score axis it autoscales to. `display` takes its settings (`mark`,
`aggregate`, `extent`, `color`, `scales`, ...), and `width` picks the zoom
level that answers ([bars, a line and density](docs/gallery.md#conservation)).

- A `MultiWiggleAdapter` config draws a row per source in `w$sources`' order.
  Reorder or filter that table to rearrange them
  ([twelve cell lines](docs/gallery.md#dnase_rows)).
- `as.matrix(w)` is the matrix the display's "Cluster rows" clusters: a row
  per source, a column per pixel.
- `w$tree <- hclust(dist(as.matrix(w)), "average")` puts the rows in the
  tree's order, and `autoplot()` draws the tree beside them
  ([clustered](docs/gallery.md#dnase_clustered)).

## Variant matrices

`jbrowse_variants()` reads a multi-sample VCF through JBrowse's variant
display: the records its filters keep, each sample's genotype at each, and
the fill the display paints it, its hue shaded by the alt dosage
([by superpopulation](docs/gallery.md#variants_by_population)).

- `display` takes the display's settings: `unit = "haplotype"` for a row per
  haplotype, `color`, `filter`, `minorAlleleFrequencyFilter`, ...
- `samples` is a table of the samples' attributes, which `rowColor` marks the
  rows by.
- Reordering or filtering `v$rows` rearranges the matrix's rows.
- `as.matrix(v)` is the genotype matrix the display's "Cluster rows by
  genotype" clusters, a column per ALT allele with `NA` for a no-call.
- A tree in `v$tree` puts the rows in its leaf order, and `autoplot()` draws
  it beside them. The display's dialog writes
  `v$tree <- hclust(dist(as.matrix(v)), "average")`, and any `hclust`, or
  anything `as.hclust()` takes, labelled by row name, will do
  ([clustered](docs/gallery.md#variants_clustered)).

## A view of several tracks

`jbrowse_stack()` stacks tracks as a JBrowse view does: each track's plot, or
an alignments track's bands, top to bottom in the order given, on one x axis,
with the keys collected.

- `jbrowse_stack()` lays several loci end to end, as a view lays out several
  regions, and labels each tick in its own region's positions
  ([two regions](docs/gallery.md#two_regions)).
- A transcript's exons, each with a few bases of flank, lay out side by side
  with the introns collapsed, as JBrowse's **Collapse introns** does, and each
  splice junction arcs from one exon straight to the next
  ([GAPDH's exons](docs/gallery.md#collapsed_introns)).
- Each band is an ordinary ggplot, so patchwork's `&` restyles them all, as
  in `jbrowse_stack(genes, signal, reads) & theme_classic()`. A ggplot of
  your own stacks as a band on the same axis.
- An alignments track read with `display = list(showCoverage = FALSE)`
  stacks no coverage band.

Each band takes the pixels the display gives it. Reads draw at the display's
`featureHeight` (7, or 3 and 1 for its compact and super-compact settings),
and `jbrowse_height()` gives the figure height that keeps those proportions,
so a compact pileup saves as a shorter figure:

```r
loci <- c("chr12:6,534,201-6,535,700", "chr12:6,537,201-6,538,700")
genes <- jbrowse_features(e("gapdh_gencode.gff3.gz"), loci)
reads <- jbrowse_alignments(e("gapdh_direct_rna.bam"), loci,
  display = list(featureHeight = 3))
fig <- jbrowse_stack(genes, reads)
ggsave("stack.png", fig, width = 10, height = jbrowse_height(fig))
```

## Synteny

`jbrowse_synteny()`, `jbrowse_dotplot()` and `jbrowse_multiway()` read
alignments (PAF, indexed PIF, MUMmer delta, UCSC chain, MashMap) through
JBrowse's synteny displays.

- A file guessed from its extension names its sides `query` and `target`. A
  PAF names its sequences by assembly when they read `sample#haplotype#contig`
  ([the amylase locus](docs/gallery.md#amylase_synteny)), and a JBrowse
  adapter config names assemblies, which `assemblies = c(query, target)` puts
  on rows.
- `swap = TRUE` puts the file's target on the top row.
- `labels = c(query = "peach", target = "grape")` writes those names in the
  axis titles, whichever side each is on.
- `guides_jbrowse()` keys a continuous field such as
  `list(field = "identity")` as the display's `0%` to `100%` bar. A threshold
  field such as `mapq` keys its bins as chips.

### Two rows

`jbrowse_synteny()` returns the blocks two rows of regions share and the
ribbons the linear synteny display draws between them, in its colors and
opacity, with each CIGAR insertion, deletion and skip wider than a pixel
([HG002's chr8 inversion](docs/gallery.md#inversion_synteny),
[a CIGAR zoom](docs/gallery.md#amylase_zoom)).

- `view` takes the `LinearSyntenyView`'s settings (`color`, `opacity`,
  `cigarMode`, `minAlignmentLength`, `drawLocationMarkers`, ...).
- `autoplot(s, above = genes, below = genes2)` stacks tracks over either row
  ([genes on both rows](docs/gallery.md#synteny_genes)).
- `geom_jbrowse_ribbons()` draws `s$ribbons`, each foot's span on its row,
  with `thin` set where the band draws a one-pixel line. `s$alignments` holds
  the blocks.
- `geom_jbrowse_offscreen_mates()` draws `s$offscreen_strip`, the marks and
  contig names the band hangs off a row's edge for blocks whose mate the other
  row does not show (`s$offscreen`). `autoplot()` keys each contig from
  `s$offscreen_contigs`.

### Dotplots

`jbrowse_dotplot()` reads the same files through JBrowse's dotplot display,
which draws a line per block, or per CIGAR step where a block is wide enough
on the plot, and `view` takes the `DotplotView`'s settings
([HG002's chr8](docs/gallery.md#inversion_dotplot)). With no regions given,
each axis shows every contig of its assembly in natural order, so a
whole-genome dotplot is one call
([peach against grape](docs/gallery.md#whole_genome)).
`geom_jbrowse_dotplot()` draws its lines.

### Several genomes

`jbrowse_multiway()` reads an all-vs-all PAF, an indexed PIF or a
multi-genome adapter config through JBrowse's multi-way synteny display. It
returns an anchor lane over the locus given, a lane per other genome in the
alignments' own coordinates, and the ribbons between adjacent lanes
([TP53 in four primates](docs/gallery.md#multiway)).

- `lanes` picks and orders the genomes, or pins a lane to a contig.
- `display` takes the display's settings (`ribbonColor`, `drawCurves`,
  `bridgeSkippedLanes`, ...).
- `autoplot(m, above = genes)` stacks tracks over the anchor lane.

## From a JBrowse view

JBrowse's linear genome view is gaining an **Export R script** menu item, not
yet in a release. It writes one ggjbrowse reader per open track, then
`jbrowse_stack()` and `ggsave(height = jbrowse_height(fig))`.

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

- The variant display's row arrangements other than clustering (bands, a
  saved order) and its zoomed-out density fade
- The bisulfite layer
- Facets on a feature track
- The amino-acid overlay on a zoomed-in CDS
- A CRAN release, which needs the bundle's own dependencies credited (`cph`
  and a LICENSE.note) and its source named

## Docs

- [docs/gallery.md](docs/gallery.md): the R and `jb2export` recipe of each
  gallery figure
- [docs/api.md](docs/api.md): every reader, helper and geom, and what each
  returns
- [docs/how-does-this-work.md](docs/how-does-this-work.md): the engine,
  reading files and who decides what
- [docs/dataflow.md](docs/dataflow.md): a diagram of one call from path to plot
- [docs/faq.md](docs/faq.md): V8 versus QuickJS, caching, indexes, coordinates
- [agent-docs/engine-owns-display-decisions.md](agent-docs/engine-owns-display-decisions.md):
  the plan to move the rules R still copies into the engine
- [CONTRIBUTING.md](CONTRIBUTING.md): development, rebuilding the engine and CI

## See also

- [JBrowse 2](https://jbrowse.org)
- [JBrowseR](https://github.com/GMOD/JBrowseR) - JBrowse 2 as an R htmlwidget

Other R tools

- [ggbandage](https://github.com/GMOD/ggbandage) - pangenome graphs as ggplot2 layers
- [ggtubemap](https://github.com/GMOD/ggtubemap) - pangenome tubemaps as ggplot2 layers

## Academic use

Written with [NHGRI](http://genome.gov) funding as part of
[JBrowse](http://jbrowse.org). If you use this in a publication, please cite the
most recent JBrowse paper at [jbrowse.org](http://jbrowse.org).

## License

Apache License 2.0 © [Colin Diesh](https://github.com/cmdcolin)
