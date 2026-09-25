#!/bin/bash

# Build a meryl k-mer database from PacBio HiFi reads for
# Merqury genome assembly evaluation.
# The k-mer size was determined using Merqury best_k.sh.

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
READS="/path/to/pacbio_hifi_reads.fasta"
KMER_SIZE="21"
OUT_DB="pacbio.meryl"

# -----------------------------
# Software requirements
# -----------------------------
# meryl

meryl count \
    k="${KMER_SIZE}" \
    "${READS}" \
    output "${OUT_DB}" \
    > meryl.log 2>&1
