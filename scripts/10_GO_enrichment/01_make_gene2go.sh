#!/bin/bash

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------

WORKDIR="/path/to/genome_comparison"
EGGNOG="/path/to/output.emapper.annotations"

OUTDIR="${WORKDIR}/go_enrichment"
OUTPUT="${OUTDIR}/gene2go.tsv"

mkdir -p "${OUTDIR}"

# Extract gene-to-GO mappings from eggNOG-mapper annotations.
#
# Column 10 of the eggNOG-mapper annotation file contains
# comma-separated GO terms.

awk -F'\t' '
BEGIN {
    OFS = "\t"
}

!/^#/ && $10 != "-" {

    split(
        $10,
        gos,
        ","
    )

    for (i in gos) {

        if (gos[i] ~ /^GO:/) {
            print $1, gos[i]
        }
    }
}
' "${EGGNOG}" > "${OUTPUT}"