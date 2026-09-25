#!/usr/bin/env bash
#SBATCH --job-name=gffread_edulis
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=gffread_edulis_%j.out
#SBATCH --error=gffread_edulis_%j.err

set -euo pipefail

# Required software:
#   gffread

# -----------------------------
# Inputs
# -----------------------------

BRAKER_DIR="/path/to/edulis_braker3"
GENOME_FA="/path/to/edulis_nuclear_only.softmasked.fa"

GFF3="${BRAKER_DIR}/braker.longest_isoform.agat.gff3"

OUT_AA="${BRAKER_DIR}/braker.longest_isoform.gffread.aa"
OUT_CDS="${BRAKER_DIR}/braker.longest_isoform.gffread.cds.fa"

# -----------------------------
# Extract protein and CDS
# sequences
# -----------------------------

gffread \
    "${GFF3}" \
    -g "${GENOME_FA}" \
    -y "${OUT_AA}" \
    -x "${OUT_CDS}"

echo "[INFO] Protein FASTA: ${OUT_AA}"
echo "[INFO] CDS FASTA: ${OUT_CDS}"
