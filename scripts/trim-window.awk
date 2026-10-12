# Trims SAM records to the reference window [ws, we] (1-based, inclusive),
# turning what falls outside into hard clips, and keeps only the SA and HP tags.
# Usage: samtools view reads.bam | awk -v ws=N -v we=N -f trim-window.awk
BEGIN { FS = OFS = "\t" }

function emit(m) {
  cigar = cigar pending m
  pending = ""
}

/^@/ { print; next }

$6 != "*" {
  cig = $6
  n = 0
  while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
    len[n] = substr(cig, 1, RLENGTH - 1) + 0
    op[n] = substr(cig, RLENGTH, 1)
    n++
    cig = substr(cig, RLENGTH + 1)
  }
  p = $4; q = 0; qlen = length($10)
  cigar = ""; pending = ""; newpos = 0; qs = -1; qe = -1; seen_m = 0
  hard_lead = 0; hard_tail = 0; lead_s = 0; trail_s = 0
  for (i = 0; i < n; i++) {
    L = len[i]; o = op[i]
    if (o == "H") {
      if (i == 0) hard_lead = L; else hard_tail = L
    } else if (o == "S") {
      if (!seen_m) lead_s = L; else trail_s = L
      q += L
    } else if (o == "I") {
      if (qe >= 0 && p <= we) pending = pending L "I"
      q += L
    } else if (o == "D" || o == "N") {
      a = (p > ws) ? p : ws; b = (p + L - 1 < we) ? p + L - 1 : we
      if (qe >= 0 && b >= a) pending = pending (b - a + 1) o
      p += L
    } else {
      seen_m = 1
      a = (p > ws) ? p : ws; b = (p + L - 1 < we) ? p + L - 1 : we
      if (b >= a) {
        if (qs < 0) { qs = q + (a - p); newpos = a }
        emit((b - a + 1) o)
        qe = q + (b - p) + 1
      }
      p += L; q += L
    }
  }
  if (qs < 0) next
  refend = p - 1
  if (newpos == $4 && lead_s > 0) { qs = 0; cigar = lead_s "S" cigar }
  if (refend <= we && trail_s > 0) { qe = qlen; cigar = cigar trail_s "S" }
  if (qs > 0) cigar = (qs + hard_lead) "H" cigar
  else if (hard_lead > 0) cigar = hard_lead "H" cigar
  if (qlen - qe > 0) cigar = cigar (qlen - qe + hard_tail) "H"
  else if (hard_tail > 0) cigar = cigar hard_tail "H"
  $4 = newpos; $6 = cigar
  $10 = substr($10, qs + 1, qe - qs)
  if ($11 != "*") $11 = substr($11, qs + 1, qe - qs)
  line = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11
  for (j = 12; j <= NF; j++) if ($j ~ /^(SA|HP):/) line = line OFS $j
  print line
}
