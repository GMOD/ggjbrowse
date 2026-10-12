#!/bin/sh
# inst/extdata's liftOver PIFs, UCSC's hg38 chains to other genomes cut to
# hg38 chr17:7,300,001-7,800,000 around TP53, in the order a multi-way view
# stacks them: chimp, gorilla, orangutan. Each keeps the rows of both tiers that
# overlap the window on hg38's side and, on the other genome's, the rows
# overlapping where those land.
# Usage: scripts/extdata-liftover.sh; needs tabix and bgzip.
set -eu
out=$(cd "$(dirname "$0")" && pwd)/../inst/extdata
base=https://jbrowse.org/ucsc/hg38/liftOver
start=7300000
end=7800000
tmp=$(mktemp -d)
cd "$tmp"

clip() {
  genome=$1
  pif=$base/hg38To$genome.over.pif.gz
  tabix -H "$pif" > "$genome.pif"
  for tier in T t; do
    tabix "$pif" "${tier}chr17:$start-$end" >> "$genome.rows"
  done
  awk -F'\t' -v OFS='\t' '{ print $6, $8, $9 }' "$genome.rows" |
    awk -F'\t' '{ lo[$1] = ($1 in lo && lo[$1] < $2) ? lo[$1] : $2
      hi[$1] = ($1 in hi && hi[$1] > $3) ? hi[$1] : $3 }
      END { for (c in lo) print c, lo[c], hi[c] }' > "$genome.partners"
  while read -r contig lo hi; do
    for tier in Q q; do
      tabix "$pif" "$tier$contig:$lo-$hi" |
        awk -F'\t' -v c="$contig" -v s=$start -v e=$end \
          '$6 == "chr17" && $9 > s && $8 < e' >> "$genome.rows" || true
    done
  done < "$genome.partners"
  sort -u "$genome.rows" | sort -k1,1 -k3,3n >> "$genome.pif"
  bgzip -c "$genome.pif" > "$out/hg38To$genome.over.pif.gz"
  tabix -f -C -0 -s1 -b3 -e4 -c '#' "$out/hg38To$genome.over.pif.gz"
}
for genome in PanTro6 GorGor6 PonAbe3; do clip "$genome"; done

cd "$out"
rm -r "$tmp"
