#!/bin/bash
#SBATCH --job-name=btk_snail
#SBATCH --output=btk_snail_%j.out
#SBATCH --error=btk_snail_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail

# Required software:
#   BlobToolKit (blobtools)

# -----------------------------
# Inputs
# -----------------------------

WORKDIR="/path/to/blobtools"

DATASET="${WORKDIR}/chilensis_blobdir"
OUTDIR="${WORKDIR}/snail_plot"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"


# -----------------------------
# Generate snail plot
# -----------------------------

blobtools view \
    --plot \
    --view snail \
    --format svg \
    "${DATASET}"
