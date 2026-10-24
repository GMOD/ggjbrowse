#!/bin/sh
# inst/extdata/gapdh_pfam.bed.gz: GAPDH's two Pfam domains, as InterPro calls
# them on UniProt P04406, placed on every GENCODE transcript in
# gapdh_gencode.gff3.gz whose CDS translates to exactly that protein.
# Usage: scripts/extdata-gapdh-pfam.sh; needs curl, jq, gffread, bgzip, tabix
# and a `jbrowse` CLI with make-protein-features ($JBROWSE overrides it).
set -eu
out=$(pwd)/inst/extdata
jbrowse=${JBROWSE:-jbrowse}
tmp=$(mktemp -d)
cd "$tmp"

zcat "$out/gapdh_gencode.gff3.gz" > genes.gff3
zcat "$out/gapdh_hg38.fa.gz" > ref.fa
gffread -y proteins.fa -g ref.fa genes.gff3
protein=$(curl -fsS https://rest.uniprot.org/uniprotkb/P04406.fasta |
  tail -n +2 | tr -d '\n')
transcripts=$(awk -v want="$protein" '
  /^>/ { if (id != "" && seq == want) print id; id = substr($1, 2); seq = ""; next }
  { seq = seq $0 }
  END { sub(/\.$/, "", seq); if (seq == want) print id }' proteins.fa)

# One InterProScan TSV row per transcript and domain: protein, md5, length,
# analysis, accession, short name, start, end, E-value, status, date, IPR,
# description
curl -fsS https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/P04406 |
  jq -r '.results[] | .metadata as $m | .proteins[].entry_protein_locations[] |
    [$m.accession, $m.name, .fragments[0].start, .fragments[0].end, .score] |
    @tsv' > domains.tsv
for t in $transcripts; do
  awk -F'\t' -v OFS='\t' -v t="$t" '{
    short = $1 == "PF00044" ? "Gp_dh_N" : $1 == "PF02800" ? "Gp_dh_C" : $1
    print t, "-", "335", "Pfam", $1, short, $3, $4, $5, "T", "-", "-", $2
  }' domains.tsv
done > calls.tsv

$jbrowse make-protein-features genes.gff3 calls.tsv \
  --out "$out/gapdh_pfam.bed.gz"
cd "$out"
rm -r "$tmp"
