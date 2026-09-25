#!/usr/bin/env bash

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------

PFAM_MAPPING="/path/to/Piriyapongsa_Pfam_name_mapping.tsv"
EGGNOG="/path/to/output.emapper.annotations"

OUTDIR="/path/to/output/TE_associated_Pfam"

mkdir -p "${OUTDIR}"

ALIASES="${OUTDIR}/Piriyapongsa_TE_Pfam_name_aliases.txt"
MATCHES="${OUTDIR}/Piriyapongsa_normalized_matches.txt"
GENES="${OUTDIR}/genes_with_Piriyapongsa_TE_Pfam_normalized.tsv"


# -----------------------------
# 1. Build Pfam name alias set
#
# Include both the historical Pfam names
# reported by Piriyapongsa et al. (2007)
# and current names resolved through the
# same Pfam accession.
# -----------------------------

tr -d '\r' < "${PFAM_MAPPING}" \
| awk -F'\t' '
NR > 1 {
    if ($2 != "")
        print $2

    if ($4 == "resolved" && $3 != "")
        print $3
}
' \
| sed '/^$/d' \
| sort -u \
> "${ALIASES}"


# -----------------------------
# 2. Identify TE-associated Pfam
# names present in eggNOG annotations
# -----------------------------

awk -F'\t' '
!/^#/ && $21 != "-" {
    n = split($21, a, ",")

    for (i = 1; i <= n; i++)
        print a[i]
}
' "${EGGNOG}" \
| sort -u \
> "${OUTDIR}/pfam_domains_unique.txt"

comm -12 \
    <(sort "${ALIASES}") \
    <(sort "${OUTDIR}/pfam_domains_unique.txt") \
    > "${MATCHES}"


# -----------------------------
# 3. Identify proteins assigned
# TE-associated Pfam annotations
#
# Column 21 of the eggNOG-mapper
# annotation file contains Pfam terms.
# -----------------------------

awk -F'\t' '
BEGIN {
    OFS = "\t"
}

NR == FNR {
    te[$1] = 1
    next
}

!/^#/ && $21 != "-" {
    n = split($21, a, ",")
    hits = ""

    for (i = 1; i <= n; i++) {

        if (a[i] in te) {

            if (hits == "")
                hits = a[i]
            else
                hits = hits "," a[i]
        }
    }

    if (hits != "")
        print $1, hits
}
' "${ALIASES}" \
  "${EGGNOG}" \
> "${GENES}"


# -----------------------------
# 4. Summary
# -----------------------------

echo "TE-associated Pfam aliases:"
wc -l "${ALIASES}"

echo
echo "TE-associated Pfam domains present in eggNOG annotations:"
wc -l "${MATCHES}"
cat "${MATCHES}"

echo
echo "Proteins with at least one matching TE-associated Pfam annotation:"
wc -l "${GENES}"

echo
echo "Pfam domain counts among matched proteins:"
cut -f2 "${GENES}" \
| tr ',' '\n' \
| sort \
| uniq -c \
| sort -nr
