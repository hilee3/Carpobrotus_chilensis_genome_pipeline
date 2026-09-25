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
GENOME_FA="/path/to/assembly_split.softmasked.fa"

# Longest-isoform BRAKER3 annotation used for sequence extraction.
# Transcript g103.t2 was removed prior to longest-isoform selection
# because of an internal-stop issue.
GFF3="${BRAKER_DIR}/braker.g103t2_removed.longest_isoform.agat.gff3"

OUT_AA="${BRAKER_DIR}/braker.g103t2_removed.longest_isoform.gffread.aa"
OUT_CDS="${BRAKER_DIR}/braker.g103t2_removed.longest_isoform.gffread.cds.fa"

# -----------------------------
# Software requirements
# -----------------------------
# gffread must be available in PATH.

command -v gffread >/dev/null || {
    echo "[ERR] gffread not found in PATH"
    exit 1
}

[[ -r "$GENOME_FA" ]] || {
    echo "[ERR] Genome FASTA not found: $GENOME_FA"
    exit 1
}

[[ -r "$GFF3" ]] || {
    echo "[ERR] GFF3 not found: $GFF3"
    exit 1
}

gffread \
    "$GFF3" \
    -g "$GENOME_FA" \
    -y "$OUT_AA" \
    -x "$OUT_CDS"

echo "[INFO] Protein FASTA: $OUT_AA"
echo "[INFO] CDS FASTA: $OUT_CDS"