#!/bin/bash
#SBATCH --job-name=agat_stats_edulis
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=agat_stats_edulis_%j.out
#SBATCH --error=agat_stats_edulis_%j.err

set -euo pipefail

# Required software:
#   AGAT

# -----------------------------
# Inputs
# -----------------------------

# Complete BRAKER3 annotation, including all transcript isoforms.
GFF="/path/to/edulis_braker3/braker.gff3"
OUTDIR="/path/to/output/agat_stats"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# -----------------------------
# Annotation statistics
# -----------------------------

agat_sp_statistics.pl \
    --gff "${GFF}" \
    -o braker_annotation_stats.txt

echo "[INFO] Output: ${OUTDIR}/braker_annotation_stats.txt"
