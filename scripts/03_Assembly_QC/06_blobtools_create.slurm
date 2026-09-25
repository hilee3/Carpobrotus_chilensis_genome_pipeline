#!/bin/bash
#SBATCH --job-name=btk_create
#SBATCH --output=btk_create_%j.out
#SBATCH --error=btk_create_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

# Required software:
#   BlobToolKit (blobtools)

# -----------------------------
# Inputs
# -----------------------------

WORKDIR="/path/to/blobtools"
ASSEMBLY="/path/to/assembly_split.fa"

BLOBDIR="${WORKDIR}/chilensis_blobdir"

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"


# -----------------------------
# Create BlobToolKit dataset
# -----------------------------

blobtools create \
    --fasta "${ASSEMBLY}" \
    "${BLOBDIR}"

blobtools validate \
    --basic \
    "${BLOBDIR}"
