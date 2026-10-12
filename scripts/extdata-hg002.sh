#!/bin/sh
# inst/extdata's HG002 chr8 haplotype files, around the 8p23.1 inversion:
#   hg002_chr8.chain          the Q100 chain's chr8_MATERNAL to chr8_PATERNAL blocks
#   hg002_mat_genes.gff3.gz   JHU Liftoff v0.6 protein-coding genes on
#   hg002_pat_genes.gff3.gz   chr8_MATERNAL and chr8_PATERNAL, 5,250,001-14,250,000
# Usage: scripts/extdata-hg002.sh; needs curl, tabix and bgzip.
set -eu
out=$(cd "$(dirname "$0")" && pwd)/../inst/extdata
s3=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies
tmp=$(mktemp -d)
cd "$tmp"

curl -sL "$s3/changes/hg002v1.2_to_other_haplotype.chain.gz" | gzip -dc |
  awk 'BEGIN { RS = ""; ORS = "\n\n" } $3 == "chr8_MATERNAL" && $8 == "chr8_PATERNAL"' |
  cat > "$out/hg002_chr8.chain"

genes() {
  hap=$1; contig=$2; name=$3
  tabix "$s3/annotation/JHULiftoff/v0.6/hg002v1.1.$hap.loff.v0.6.gff.gz" \
    "$contig:5250001-14250000" > "$name.all"
  awk -F'\t' '
    function stem(attrs,   s) {
      if (!match(attrs, /(ID|Parent)=[^;]*/)) return ""
      s = substr(attrs, RSTART, RLENGTH)
      sub(/^[A-Za-z]+=/, "", s)
      sub(/\..*/, "", s)
      return s
    }
    NR == FNR { if ($3 == "gene" && $9 ~ /gene_biotype=protein_coding/) keep[stem($9)] = 1; next }
    stem($9) in keep' "$name.all" "$name.all" | sort -k1,1 -k4,4n > "$name.gff"
  { echo "##gff-version 3"; cat "$name.gff"; } | bgzip > "$out/$name.gff3.gz"
  tabix -f -p gff "$out/$name.gff3.gz"
}
genes MAT chr8_MATERNAL hg002_mat_genes
genes PAT chr8_PATERNAL hg002_pat_genes

cd "$out"
rm -r "$tmp"
