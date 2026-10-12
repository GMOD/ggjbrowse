#!/bin/sh
# inst/extdata's human BigWigs, hg38 chr12:6,532,001-6,541,000 around GAPDH:
#   gapdh_phylop100way.bw   UCSC's phyloP 100-way conservation, signed per base
#   gapdh_h3k27ac_*.bw      ENCODE Broad H3K27ac signal in seven cell lines
# Usage: scripts/extdata-bigwig.sh; needs bigWigToBedGraph and bedGraphToBigWig.
set -eu
out=$(cd "$(dirname "$0")" && pwd)/../inst/extdata
ucsc=https://hgdownload.soe.ucsc.edu
start=6532000
end=6541000
tmp=$(mktemp -d)
printf 'chr12\t133275309\n' > "$tmp/chrom.sizes"

clip() {
  bigWigToBedGraph -chrom=chr12 -start=$start -end=$end "$1" stdout |
    sort -k1,1 -k2,2n > "$tmp/clip.bedGraph"
  bedGraphToBigWig "$tmp/clip.bedGraph" "$tmp/chrom.sizes" "$out/$2"
}

clip "$ucsc/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw" gapdh_phylop100way.bw
for cell in Gm12878 H1hesc Hsmm Huvec K562 Nhek Nhlf; do
  clip "$ucsc/gbdb/hg38/bbi/wgEncodeReg/wgEncodeRegMarkH3k27ac/wgEncodeBroadHistone${cell}H3k27acStdSig.bigWig" \
    "gapdh_h3k27ac_$(echo "$cell" | tr '[:upper:]' '[:lower:]').bw"
done
rm -r "$tmp"
