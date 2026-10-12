#!/bin/sh
# inst/extdata/1000g_dcst2.vcf.gz and its samples file: 100 people of the 1000
# Genomes phase 3 (GRCh37), 20 from each superpopulation, at
# 1:155,000,001-155,050,000 (DCST2, DCST1, ADAM15), from the slice and the
# population list in a jbrowse-components checkout.
# Usage: scripts/extdata-1000g.sh [jbrowse-components]; needs bcftools, tabix.
set -eu
jbc=${1:-../jbrowse-components}
vcf=$jbc/test_data/1000g_snpeff_chr1/ALL.chr1.phase3.155M.snpeff.vcf.gz
out=inst/extdata/1000g_dcst2
tmp=$(mktemp -d)
bcftools query -l "$vcf" > "$tmp/in_vcf"
Rscript scripts/extdata-1000g.R "$jbc/test_data/1000g_cnv/samples.tsv" \
  "$tmp/in_vcf" > "${out}_samples.tsv"
tail -n +2 "${out}_samples.tsv" | cut -f1 > "$tmp/picked"
bcftools view -s "$(paste -sd, "$tmp/picked")" -r 1:155000001-155050000 \
  --min-ac 1 -Ou "$vcf" | bcftools annotate -x INFO -Oz -o "$out.vcf.gz"
tabix -f -p vcf "$out.vcf.gz"
rm -r "$tmp"
