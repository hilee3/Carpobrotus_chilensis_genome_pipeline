#!/bin/bash
#SBATCH --job-name=gffread
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=gffread_%j.out
#SBATCH --error=gffread_%j.err

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
BRAKER_DIR="/path/to/braker3_output"
GENOME_FA="/path/to/softmasked_genome.fa"

# Curated longest-isoform BRAKER3 annotation.
# One transcript model (g103.t2) containing an internal stop codon
# was removed prior to longest isoform selection.
GFF3="${BRAKER_DIR}/braker.g103t2_removed.longest_isoform.agat.gff3"

OUT_AA="${BRAKER_DIR}/braker.longest_isoform.aa"
OUT_CDS="${BRAKER_DIR}/braker.longest_isoform.cds.fa"

# -----------------------------
# Software requirements
# -----------------------------
# gffread must be available in PATH.

command -v gffread >/dev/null || {
    echo "[ERR] gffread not found in PATH"
    exit 1
}

[[ -r "$GENOME_FA" ]] || {
    echo "[ERR] Genome FASTA not found"
    exit 1
}

[[ -r "$GFF3" ]] || {
    echo "[ERR] GFF3 not found"
    exit 1
}

gffread \
    "$GFF3" \
    -g "$GENOME_FA" \
    -y "$OUT_AA" \
    -x "$OUT_CDS"

echo "[INFO] Protein FASTA: $OUT_AA"
echo "[INFO] CDS FASTA: $OUT_CDS"
