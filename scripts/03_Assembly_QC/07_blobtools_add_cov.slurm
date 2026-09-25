#!/bin/bash
#SBATCH --job-name=btk_cov
#SBATCH --output=btk_cov_%j.out
#SBATCH --error=btk_cov_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

# Required software:
#   BlobToolKit (blobtools)

# -----------------------------
# Inputs
# -----------------------------

WORKDIR="/path/to/blobtools"

BLOBDIR="${WORKDIR}/chilensis_blobdir"
BAM="${WORKDIR}/hifi_cov_split_sorted.bam"

cd "${WORKDIR}"


# -----------------------------
# Add HiFi read coverage
# -----------------------------

blobtools add \
    --cov "${BAM}" \
    "${BLOBDIR}"

blobtools validate \
    --basic \
    "${BLOBDIR}"
