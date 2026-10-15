# ggjbrowse

> [!CAUTION]
> Very early results from work in progress!

Genome browser tracks as ggplot2 layers.

- **In:** a region of a BAM, CRAM, GFF3, GTF, BED, BigBed, BigWig or VCF, or
  two regions a PAF aligns.
- **Out:** the tables the [JBrowse 2](https://jbrowse.org) genome browser's
  own display draws from, run in V8: reads on pileup rows, gene glyphs,
  binned scores, genotypes and synteny ribbons, with the colors the display
  paints them.
- **Also:** reads CRAM, which Rsamtools cannot, and needs no Bioconductor
  packages.
- **Drawn with:** a geom per table, so a figure is an ordinary ggplot that
  `+` changes piece by piece. `autoplot()` composes the geoms into the
  browser's track.

```r
remotes::install_github("GMOD/ggjbrowse")
```

## Gallery

- Every figure is drawn twice from data installed with the package, by
  JBrowse's own renderer (`jb2export` from `@jbrowse/img`) on the left and by
  ggjbrowse on the right.
- A title links to the figure's R and `jb2export` recipes in
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
<tr><td width="50%"><img src="gallery/compact_reads-jbrowse.webp" alt="Compact reads over the same deletion in JBrowse. Sorted by the base at the center, at 3 pixels a row the pileup takes far less height, the deletion's dark bars banded through it." width="100%"></td><td width="50%"><img src="gallery/compact_reads-r.webp" alt="Compact reads over the same deletion drawn by ggjbrowse. Sorted by the base at the center, at 3 pixels a row the pileup takes far less height, the deletion's dark bars banded through it." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#base_zoom">Base-level zoom at the deletion's left edge</a></b> <code>chr10:122,835,300-122,835,420</code></td></tr>
<tr><td width="50%"><img src="gallery/base_zoom-jbrowse.webp" alt="Base-level zoom at the deletion's left edge in JBrowse. At 120 bp the pileup writes out each mismatch's base, and every read sorted to the top carries a G a few bases before its deletion." width="100%"></td><td width="50%"><img src="gallery/base_zoom-r.webp" alt="Base-level zoom at the deletion's left edge drawn by ggjbrowse. At 120 bp the pileup writes out each mismatch's base, and every read sorted to the top carries a G a few bases before its deletion." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#inversion">Split reads chained across an inversion, with their arcs</a></b> <code>chr1:197,786,900-197,789,700</code></td></tr>
<tr><td width="50%"><img src="gallery/inversion-jbrowse.webp" alt="Split reads chained across an inversion, with their arcs in JBrowse. The split reads, in their own section above the unsplit ones, are colored inverted between the two breakpoints, where inversion-type arcs join their segments." width="100%"></td><td width="50%"><img src="gallery/inversion-r.webp" alt="Split reads chained across an inversion, with their arcs drawn by ggjbrowse. The split reads, in their own section above the unsplit ones, are colored inverted between the two breakpoints, where inversion-type arcs join their segments." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#methylation">CpG methylation from MM/ML tags</a></b> <code>20:10,001-22,000</code></td></tr>
<tr><td width="50%"><img src="gallery/methylation-jbrowse.webp" alt="CpG methylation from MM/ML tags in JBrowse. Each read's 5mC calls in red where methylated and blue where not." width="100%"></td><td width="50%"><img src="gallery/methylation-r.webp" alt="CpG methylation from MM/ML tags drawn by ggjbrowse. Each read's 5mC calls in red where methylated and blue where not." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#bisulfite">CHH methylation from bisulfite reads</a></b> <code>NC_003070.9:4,398,001-4,412,000</code></td></tr>
<tr><td width="50%"><img src="gallery/bisulfite-jbrowse.webp" alt="CHH methylation from bisulfite reads in JBrowse. RefSeq genes and RepeatMasker over MethylDackel's three context bigWigs and each read's methylated CHH cytosines in red: CHH is dense over the LTR/Copia element on the right and rare over the gene on the left, which carries CpG methylation only." width="100%"></td><td width="50%"><img src="gallery/bisulfite-r.webp" alt="CHH methylation from bisulfite reads drawn by ggjbrowse. RefSeq genes and RepeatMasker over MethylDackel's three context bigWigs and each read's methylated CHH cytosines in red: CHH is dense over the LTR/Copia element on the right and rare over the gene on the left, which carries CpG methylation only." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#dnase_rows">DNase in twelve cell lines over the beta-globin locus</a></b> <code>chr11:5,270,001-5,320,000</code></td></tr>
<tr><td width="50%"><img src="gallery/dnase_rows-jbrowse.webp" alt="DNase in twelve cell lines over the beta-globin locus in JBrowse. Most lines share one peak near 5,281,000, and K562 alone shows several more to either side." width="100%"></td><td width="50%"><img src="gallery/dnase_rows-r.webp" alt="DNase in twelve cell lines over the beta-globin locus drawn by ggjbrowse. Most lines share one peak near 5,281,000, and K562 alone shows several more to either side." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#dnase_clustered">DNase in twelve cell lines as density, clustered by signal</a></b> <code>chr11:5,270,001-5,320,000</code></td></tr>
<tr><td width="50%"><img src="gallery/dnase_clustered-jbrowse.webp" alt="DNase in twelve cell lines as density, clustered by signal in JBrowse. K562 splits off alone from the other eleven lines." width="100%"></td><td width="50%"><img src="gallery/dnase_clustered-r.webp" alt="DNase in twelve cell lines as density, clustered by signal drawn by ggjbrowse. K562 splits off alone from the other eleven lines." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#conservation">phyloP conservation across GAPDH as bars, a line and density</a></b> <code>chr12:6,534,201-6,538,700</code></td></tr>
<tr><td width="50%"><img src="gallery/conservation-jbrowse.webp" alt="phyloP conservation across GAPDH as bars, a line and density in JBrowse. The exons stand out as runs of high scores in bars, in a mean line over its min-max band, and in a density strip." width="100%"></td><td width="50%"><img src="gallery/conservation-r.webp" alt="phyloP conservation across GAPDH as bars, a line and density drawn by ggjbrowse. The exons stand out as runs of high scores in bars, in a mean line over its min-max band, and in a density strip." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#two_regions">Genes, conservation and strand-colored reads over two regions end to end</a></b> <code>chr12:6,534,201-6,535,700 chr12:6,537,201-6,538,700</code></td></tr>
<tr><td width="50%"><img src="gallery/two_regions-jbrowse.webp" alt="Genes, conservation and strand-colored reads over two regions end to end in JBrowse. A grey divider separates the two loci, and each tick gives a position in its own region." width="100%"></td><td width="50%"><img src="gallery/two_regions-r.webp" alt="Genes, conservation and strand-colored reads over two regions end to end drawn by ggjbrowse. A grey divider separates the two loci, and each tick gives a position in its own region." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#amino_acids">Amino acids over GAPDH's first coding exons, read in frame across the intron</a></b> <code>chr12:6,534,811-6,534,880 chr12:6,536,476-6,536,545</code></td></tr>
<tr><td width="50%"><img src="gallery/amino_acids-jbrowse.webp" alt="Amino acids over GAPDH's first coding exons, read in frame across the intron in JBrowse. GAPDH's coding transcripts open MGKVKVGVNG, the last codon split by the intron and lettered red on each side, the frame running on in the next exon." width="100%"></td><td width="50%"><img src="gallery/amino_acids-r.webp" alt="Amino acids over GAPDH's first coding exons, read in frame across the intron drawn by ggjbrowse. GAPDH's coding transcripts open MGKVKVGVNG, the last codon split by the intron and lettered red on each side, the frame running on in the next exon." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#variants_by_population">100 people's genotypes across DCST2 (GRCh37), colored by superpopulation</a></b> <code>1:155,000,001-155,050,000</code></td></tr>
<tr><td width="50%"><img src="gallery/variants_by_population-jbrowse.webp" alt="100 people's genotypes across DCST2 (GRCh37), colored by superpopulation in JBrowse. 20 people from each superpopulation genotyped at 249 variants, a column per record, the East Asians homozygous for the alt allele at many of the first 40." width="100%"></td><td width="50%"><img src="gallery/variants_by_population-r.webp" alt="100 people's genotypes across DCST2 (GRCh37), colored by superpopulation drawn by ggjbrowse. 20 people from each superpopulation genotyped at 249 variants, a column per record, the East Asians homozygous for the alt allele at many of the first 40." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#variants_clustered">100 people's genotypes across DCST2 clustered, the tree beside the rows</a></b> <code>1:155,000,001-155,050,000</code></td></tr>
<tr><td width="50%"><img src="gallery/variants_clustered-jbrowse.webp" alt="100 people's genotypes across DCST2 clustered, the tree beside the rows in JBrowse. Most East Asians gather near the top, the Africans and Europeans further down." width="100%"></td><td width="50%"><img src="gallery/variants_clustered-r.webp" alt="100 people's genotypes across DCST2 clustered, the tree beside the rows drawn by ggjbrowse. Most East Asians gather near the top, the Africans and Europeans further down." width="100%"></td></tr>
<tr><td colspan="2"><b><a href="docs/gallery.md#variants_banded">100 people's genotypes across DCST2 in a band per superpopulation, each band clustered</a></b> <code>1:155,000,001-155,050,000</code></td></tr>
<tr><td width="50%"><img src="gallery/variants_banded-jbrowse.webp" alt="100 people's genotypes across DCST2 in a band per superpopulation, each band clustered in JBrowse. Each superpopulation clusters apart in its own band, and the records crowding a pixel fade by how many share it." width="100%"></td><td width="50%"><img src="gallery/variants_banded-r.webp" alt="100 people's genotypes across DCST2 in a band per superpopulation, each band clustered drawn by ggjbrowse. Each superpopulation clusters apart in its own band, and the records crowding a pixel fade by how many share it." width="100%"></td></tr>
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

`a` is a list of data frames, so the figure is also ordinary ggplot2:

```r
library(ggplot2)

ggplot(a) +
  geom_jbrowse_reads(aes(fill = mapq)) +
  geom_jbrowse_mismatches() +
  scale_fill_viridis_c() +
  scale_y_reverse()
```

- **Data:** every example runs as written on clips of public human data
  installed with the package. `scripts/extdata-*.sh` make the clips.
- **Regions:** a locus as JBrowse's location box writes it, 1-based
  (`"chr12:6,534,201-6,538,700"`), or `"chr12", 6534200, 6538700` in 0-based
  half-open coordinates. A vector of several loci draws them end to end.
- **Files:** a path or a URL. CRAM needs a `reference`.
- **Contig names:** where files spell a contig differently (`chr7` in a BAM,
  `7` in the FASTA), give each reader `assembly = jbrowse_assembly(fasta,
  aliases)`, with an alias table such as UCSC's `chromAlias.txt`. A region
  may then name the contig either way, and the reader names it in its tables
  as the region does.

## Tutorials

Whole analyses done in R, each on the data and question of a
[JBrowse tutorial](https://jbrowse.org/jb2/docs/tutorials/):

- [Methylation (bisulfite) in R](docs/tutorials/bisulfite.md): per-read CpG,
  CHG and CHH calls from an _Arabidopsis_ WGBS run, counted over a gene and a
  transposon and checked against MethylDackel
- [Methylation (long-read) in R](docs/tutorials/methylation.md): HG002
  nanopore 5mC split by haplotype at the _SNRPN_ imprinting center, checked
  against modkit and smoothed across the imprinted domain
- [Differential transcript usage in R](docs/tutorials/dtu.md): satuRn on
  ENCODE muscle and liver quantifications, written into GENCODE and drawn over
  each tissue's coverage, checked at the exon that tells _ATP5F1C_'s isoforms
  apart
- [Selection scans (Drosophila DGRP) in R](docs/tutorials/population_genomics.md):
  Fst, π and Tajima's D at the _Cyp6g1_ sweep and the `In(2L)t` inversion,
  each ranked against the rest of the genome
- [A loss-of-function allele across breeds (Dog10K) in R](docs/tutorials/dog10k_lof.md):
  the _CYP1A2_ stop codon derived from the reference, its frequency per breed
  with wolves and coyotes as controls, and the gene's copy number
- [RNA-seq splicing and strand in R](docs/tutorials/rnaseq.md): _ACTB_'s splice junctions from
  the reads, strand-specific coverage, IsoSeq isoforms, and a regtools junction
  table checked against the engine's own counts
- [Complex rearrangements and derivative alleles (COLO829) in R](docs/tutorials/cancer_sv.md):
  the chr3, chr10 and chr12 rearrangement chain walked from the calls, split
  reads counted in tumor and normal, and the der(3) allele against the truth
  set

## Docs

- [docs/alignments.md](docs/alignments.md): BAM and CRAM reads
  (`jbrowse_alignments()`): color, filter, sort, split and chain reads,
  methylation
- [docs/features.md](docs/features.md): GFF3, GTF, BED and BigBed
  (`jbrowse_features()`): gene glyphs, labels, facets
- [docs/wiggle.md](docs/wiggle.md): BigWig (`jbrowse_wiggle()`): bars, lines,
  density, rows per source, clustering
- [docs/variants.md](docs/variants.md): VCF (`jbrowse_variants()`): genotypes
  by sample, bands, clustering
- [docs/synteny.md](docs/synteny.md): PAF, PIF, delta, chain and MashMap
  (`jbrowse_synteny()`, `jbrowse_dotplot()`, `jbrowse_multiway()`): ribbons,
  dotplots, several genomes
- [docs/stack.md](docs/stack.md): several tracks (`jbrowse_stack()`):
  end-to-end regions, collapsed introns, figure size, the **Export R script**
  JBrowse menu item
- [docs/grammar.md](docs/grammar.md): the geoms, scales and facets that
  `autoplot()` composes
- [docs/api.md](docs/api.md): every reader, helper and geom, and what each
  returns
- [docs/gallery.md](docs/gallery.md): the R and `jb2export` code of each figure
- [docs/faq.md](docs/faq.md): how the engine works, what one call does, common
  errors
- [CONTRIBUTING.md](CONTRIBUTING.md): setup, rebuilding the engine, CI

## Not yet

- A CRAN submission: `R CMD check --as-cran` leaves the new-submission NOTE,
  and the 7.9 MB installed size, 4.3 MB of it `inst/extdata`, may need
  trimming or a word to CRAN

## See also

- [JBrowse 2](https://jbrowse.org)
- [JBrowseR](https://github.com/GMOD/JBrowseR) - JBrowse 2 as an R htmlwidget

Other R tools

- [ggbandage](https://github.com/GMOD/ggbandage) - pangenome graphs as ggplot2 layers
- [ggtubemap](https://github.com/GMOD/ggtubemap) - pangenome tubemaps as ggplot2 layers

## Academic use

Written with [NHGRI](https://www.genome.gov) funding as part of
[JBrowse](https://jbrowse.org). If you use this in a publication, please cite the
most recent JBrowse paper at [jbrowse.org](https://jbrowse.org).

## License

Apache License 2.0 © [Colin Diesh](https://github.com/cmdcolin)
