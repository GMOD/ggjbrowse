#!/bin/sh
# inst/extdata's human BigWigs, clipped from UCSC's hg38 files:
#   gapdh_phylop100way.bw   phyloP 100-way conservation, signed per base,
#                           chr12:6,532,001-6,541,000 around GAPDH
#   hbb_dnase_*.bw          ENCODE DNase signal in twelve cell lines,
#                           chr11:5,270,001-5,320,000 over the beta-globin
#                           locus control region, strongest in K562
# Usage: scripts/extdata-bigwig.sh; needs bigWigToBedGraph and bedGraphToBigWig.
set -eu
out=$(cd "$(dirname "$0")" && pwd)/../inst/extdata
ucsc=https://hgdownload.soe.ucsc.edu
tmp=$(mktemp -d)
printf 'chr11\t135086622\nchr12\t133275309\n' > "$tmp/chrom.sizes"

clip() {
  bigWigToBedGraph -chrom="$3" -start="$4" -end="$5" "$1" stdout |
    sort -k1,1 -k2,2n > "$tmp/clip.bedGraph"
  bedGraphToBigWig "$tmp/clip.bedGraph" "$tmp/chrom.sizes" "$out/$2"
}

clip "$ucsc/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw" gapdh_phylop100way.bw \
  chr12 6532000 6541000
for cell in Gm12878 K562 Hepg2 Huvec Hmec Nhek Helas3 A549 Mcf7 Hsmm Hcm Hct116; do
  clip "$ucsc/gbdb/hg38/bbi/wgEncodeRegDnase/wgEncodeRegDnaseUw${cell}Signal.bw" \
    "hbb_dnase_$(echo "$cell" | tr '[:upper:]' '[:lower:]').bw" chr11 5270000 5320000
done
rm -r "$tmp"
